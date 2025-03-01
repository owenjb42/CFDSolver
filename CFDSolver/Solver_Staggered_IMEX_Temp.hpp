#pragma once

#include <vector>
#include <iostream>
#include <cmath>
#include <cmath>
#include <limits>
#include <thread>
#include <algorithm>

#include "PhysicalField.hpp"
#include "FlagField.hpp"
#include "Interface.hpp"
#include "BoundaryControl.hpp"

class SolverStaggeredIMEXTemp
{
public:
    SolverStaggeredIMEXTemp(Interface& interface) : interface(interface), nx(interface.nx), ny(interface.ny), dx(interface.dx), dy(interface.dy)
    {
        SetBlockedFaces(interface);
        SetBoundaries(interface);
        interface.SetData(*this);
    }
    bool is_solving{ true };
    Interface& interface;

    void ComputeDivergence() 
    {
        for (int i = 0; i < nx; ++i) 
        {
            for (int j = 0; j < ny; ++j) 
            {
                divergence(i, j) = ((u(i + 1, j) - u(i, j)) * dy + 
                                    (v(i, j + 1) - v(i, j)) * dx +
                                    p_scr(i,j));
            }
        }
    }

    void ComputeConvectiveTerms()// check density it should be added 
    {
        double flux = 0.0;
        double value = 0.0;
        double coeff = 0.0;

        for (int i = 1; i < nx; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                if (!u_face_flags(i, j, Flag::Open))
                    continue;

                if (u_face_flags(i - 1, j, Flag::Open))
                {
                    // Inline Left: u flux
                    flux = dy * ((u(i, j) + u(i - 1, j)) / 2.0);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? u(i - 1, j) : u(i, j);
                        coeff = fluid.density * flux / (dx * dy);
                        u_coeff(i, j) += coeff;
                        u_scr(i, j) += coeff * value;
                    }
                }

                if (u_face_flags(i + 1, j, Flag::Open))
                {
                    // Inline Right: u flux
                    flux = dy * -((u(i, j) + u(i + 1, j)) / 2.0);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? u(i + 1, j) : u(i, j);
                        coeff = fluid.density * flux / (dx * dy);
                        u_coeff(i, j) += coeff;
                        u_scr(i, j) += coeff * value;
                    }
                }

                if (j != 0 && v_face_flags(i, j, Flag::Open) && v_face_flags(i - 1, j, Flag::Open) && u_face_flags(i, j - 1, Flag::Open))
                {
                    // Bottom Left: u flux
                    flux = (dx / 2.0) * v(i - 1, j);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? u(i, j - 1) : u(i, j);
                        coeff = fluid.density * flux / (dx * dy);
                        u_coeff(i, j) += coeff;
                        u_scr(i, j) += coeff * value;
                    }

                    // Bottom Right: u flux
                    flux = (dx / 2.0) * v(i, j);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? u(i, j - 1) : u(i, j);
                        coeff = fluid.density * flux / (dx * dy);
                        u_coeff(i, j) += coeff;
                        u_scr(i, j) += coeff * value;
                    }
                }

                if (j != ny - 1 && v_face_flags(i, j + 1, Flag::Open) && v_face_flags(i - 1, j + 1, Flag::Open) && u_face_flags(i, j + 1, Flag::Open))
                {
                    // Top Left: u flux
                    flux = -(dx / 2.0) * v(i - 1, j + 1);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? u(i, j + 1) : u(i, j);
                        coeff = fluid.density * flux / (dx * dy);
                        u_coeff(i, j) += coeff;
                        u_scr(i, j) += coeff * value;
                    }

                    // Top Right: u flux
                    flux = -(dx / 2.0) * v(i, j + 1);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? u(i, j + 1) : u(i, j);
                        coeff = fluid.density * flux / (dx * dy);
                        u_coeff(i, j) += coeff;
                        u_scr(i, j) += coeff * value;
                    }
                }
            }
        }

        for (int i = 0; i < nx; ++i)
        {
            for (int j = 1; j < ny; ++j)
            {
                if (!v_face_flags(i, j, Flag::Open))
                    continue;

                if (v_face_flags(i, j - 1, Flag::Open))
                {
                    // Inline Left: v flux
                    flux = dy * ((v(i, j) + v(i, j - 1)) / 2.0);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? v(i, j - 1) : v(i, j);
                        coeff = fluid.density * flux / (dx * dy);
                        v_coeff(i, j) += coeff;
                        v_scr(i, j) += coeff * value;
                    }
                }

                if (v_face_flags(i, j + 1, Flag::Open))
                {
                    // Inline Right: v flux
                    flux = dy * -((v(i, j) + v(i, j + 1)) / 2.0);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? v(i, j + 1) : v(i, j);
                        coeff = fluid.density * flux / (dx * dy);
                        v_coeff(i, j) += coeff;
                        v_scr(i, j) += coeff * value;
                    }
                }

                if (i != 0 && u_face_flags(i + 1, j, Flag::Open) && u_face_flags(i + 1, j - 1, Flag::Open) && v_face_flags(i + 1, j, Flag::Open))
                {
                    // Bottom Left: v flux
                    flux = (dy / 2.0) * u(i, j - 1);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? v(i - 1, j) : v(i, j);
                        coeff = fluid.density * flux / (dx * dy);
                        v_coeff(i, j) += coeff;
                        v_scr(i, j) += coeff * value;
                    }

                    // Bottom Right: v flux
                    flux = (dy / 2.0) * u(i, j);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? v(i - 1, j) : v(i, j);
                        coeff = fluid.density * flux / (dx * dy);
                        v_coeff(i, j) += coeff;
                        v_scr(i, j) += coeff * value;
                    }
                }

                if (i != nx - 1 && u_face_flags(i, j, Flag::Open) && u_face_flags(i, j - 1, Flag::Open) && v_face_flags(i - 1, j, Flag::Open))
                {
                    // Top Left: v flux
                    flux = -(dy / 2.0) * u(i + 1, j - 1);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? v(i + 1, j) : v(i, j);
                        coeff = fluid.density * flux / (dx * dy);
                        v_coeff(i, j) += coeff;
                        v_scr(i, j) += coeff * value;
                    }

                    // Top Right: v flux
                    flux = -(dy / 2.0) * u(i + 1, j);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? v(i + 1, j) : v(i, j);
                        coeff = fluid.density * flux / (dx * dy);
                        v_coeff(i, j) += coeff;
                        v_scr(i, j) += coeff * value;
                    }
                }
            }
        }
    }

    void ComputeDiffusiveTerms()
    {
        double coeff = 0.0;
        double value = 0.0;
        for (int i = 1; i < nx; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                if (!u_face_flags(i, j, Flag::Open))
                    continue;

                if (u_face_flags(i - 1, j, Flag::Open))
                {
                    // Inline Left: u diff
                    coeff = fluid.kinematic_viscosity / (dx * dx);
                    value = u(i - 1, j);
                    u_coeff(i, j) += coeff;
                    u_scr(i, j) += coeff * value;
                }

                if (u_face_flags(i + 1, j, Flag::Open))
                {
                    // Inline Right: u diff
                    coeff = fluid.kinematic_viscosity / (dx * dx);
                    value = u(i + 1, j);
                    u_coeff(i, j) += coeff;
                    u_scr(i, j) += coeff * value;
                }

                if (j != 0 && v_face_flags(i, j, Flag::Open) && v_face_flags(i - 1, j, Flag::Open) && u_face_flags(i, j - 1, Flag::Open))
                {
                    // Bottom: u diff
                    coeff = fluid.kinematic_viscosity / (dy * dy);
                    value = u(i, j - 1);
                    u_coeff(i, j) += coeff;
                    u_scr(i, j) += coeff * value;
                }

                if (j != ny - 1 && v_face_flags(i, j + 1, Flag::Open) && v_face_flags(i - 1, j + 1, Flag::Open) && u_face_flags(i, j + 1, Flag::Open))
                {
                    // Top: u diff
                    coeff = fluid.kinematic_viscosity / (dy * dy);
                    value = u(i, j + 1);
                    u_coeff(i, j) += coeff;
                    u_scr(i, j) += coeff * value;
                }
            }
        }

        for (int i = 0; i < nx; ++i)
        {
            for (int j = 1; j < ny; ++j)
            {
                if (!v_face_flags(i, j, Flag::Open))
                    continue;

                if (v_face_flags(i, j - 1, Flag::Open))
                {
                    // Inline Left: v diff
                    coeff = fluid.kinematic_viscosity / (dy * dy);
                    value = v(i, j - 1);
                    v_coeff(i, j) += coeff;
                    v_scr(i, j) += coeff * value;
                }

                if (v_face_flags(i, j + 1, Flag::Open))
                {
                    // Inline Right: v diff
                    coeff = fluid.kinematic_viscosity / (dy * dy);
                    value = v(i, j + 1);
                    v_coeff(i, j) += coeff;
                    v_scr(i, j) += coeff * value;
                }

                if (i != 0 && u_face_flags(i, j, Flag::Open) && u_face_flags(i, j - 1, Flag::Open) && v_face_flags(i - 1, j, Flag::Open))
                {
                    // Bottom: v diff
                    coeff = fluid.kinematic_viscosity / (dx * dx);
                    value = v(i - 1, j);
                    v_coeff(i, j) += coeff;
                    v_scr(i, j) += coeff * value;
                }

                if (i != nx - 1 && u_face_flags(i + 1, j, Flag::Open) && u_face_flags(i + 1, j - 1, Flag::Open) && v_face_flags(i + 1, j, Flag::Open))
                {
                    // Top: v diff
                    coeff = fluid.kinematic_viscosity/ (dx * dx);
                    value = v(i + 1, j);
                    v_coeff(i, j) += coeff;
                    v_scr(i, j) += coeff * value;
                }
            }
        }
    }

    void ComputeConvectiveFluxesAndDiffusiveTermsForTemperature()
    {
        double value = 0.0;
        double coeff = 0.0;
        for (int i = 0; i < nx; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                if (u_face_flags(i, j, Flag::Open))
                {
                    value = t(i - 1, j);

                    // Left: Diff
                    coeff = (fluid.conductivity / fluid.cp) / (dx * dx);
                    t_coeff(i, j) += coeff;
                    t_scr(i, j) += coeff * value;

                    // Left: Cnv
                    coeff = fluid.density * dy * u(i, j) / (dy * dx);
                    if (coeff > 0.0)
                    {
                        t_coeff(i, j) += coeff;
                        t_scr(i, j) += coeff * value;
                    }
                }

                if (u_face_flags(i + 1, j, Flag::Open))
                {
                    value = t(i + 1, j);

                    // Right: Diff
                    coeff = (fluid.conductivity / fluid.cp) / (dx * dx);
                    t_coeff(i, j) += coeff;
                    t_scr(i, j) += coeff * value;

                    // Right: Cnv
                    coeff = -fluid.density * dy * u(i + 1, j) / (dy * dx);
                    if (coeff > 0.0)
                    {
                        t_coeff(i, j) += coeff;
                        t_scr(i, j) += coeff * value;
                    }
                }

                if (v_face_flags(i, j, Flag::Open))
                {
                    value = t(i, j - 1);

                    // Bottom: Diff
                    coeff = (fluid.conductivity / fluid.cp) / (dy * dy);
                    t_coeff(i, j) += coeff;
                    t_scr(i, j) += coeff * value;

                    // Bottom: Cnv
                    coeff = fluid.density * dy * v(i, j) / (dy * dx);
                    if (coeff > 0.0)
                    {
                        t_coeff(i, j) += coeff;
                        t_scr(i, j) += coeff * value;
                    }
                }

                if (v_face_flags(i, j + 1, Flag::Open))
                {
                    value = t(i, j + 1);

                    // Top: Diff
                    coeff = (fluid.conductivity / fluid.cp) / (dy * dy);
                    t_coeff(i, j) += coeff;
                    t_scr(i, j) += coeff * value;

                    // Top: Cnv
                    coeff = -fluid.density * dy * v(i, j + 1) / (dy * dx);
                    if (coeff > 0.0)
                    {
                        t_coeff(i, j) += coeff;
                        t_scr(i, j) += coeff * value;
                    }
                }
            }
        }
    }

    void CalculateCellVelocities()
    {
        for (int i = 0; i < nx; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                cell_u(i, j) = (u(i, j) + u(i + 1, j)) / 2.0;
                cell_v(i, j) = (v(i, j) + v(i, j + 1)) / 2.0;
            }
        }

        CorrectBoundaryCellVelocities();
    }

    void SolvePressure()
    {
        for (int itter = 0; itter < innerPressureItterations; ++itter)
        {
            std::swap(p_correction_old, p_correction);
            for (int i = 0; i < nx; ++i)
            {
                for (int j = 0; j < ny; ++j)
                {
                    p_correction(i, j) = 0.0;

                    double Ap = u_coeff(i, j) + u_coeff(i + 1, j) + v_coeff(i, j) + v_coeff(i, j + 1) + p_coef(i, j);

                    if (Ap == 0) { continue; }

                    if (u_face_flags(i, j, Flag::Open))
                        p_correction(i, j) += (u_coeff(i, j) / Ap) * p_correction_old(i - 1, j);
                    if (u_face_flags(i + 1, j, Flag::Open))
                        p_correction(i, j) += (u_coeff(i + 1, j) / Ap) * p_correction_old(i + 1, j);
                    if (v_face_flags(i, j, Flag::Open))
                        p_correction(i, j) += (v_coeff(i, j) / Ap) * p_correction_old(i, j - 1);
                    if (v_face_flags(i, j + 1, Flag::Open))
                        p_correction(i, j) += (v_coeff(i, j + 1) / Ap) * p_correction_old(i, j + 1);

                    p_correction(i, j) -= divergence(i, j);
                }
            }
        }

        for (int i = 0; i < nx; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                p(i, j) += relaxation * p_correction(i, j);
            }
        }
    }

    void ApplyPressureCorrectionToVelocity()
    {
        for (int i = 1; i < nx; ++i) 
        {
            for (int j = 0; j < ny; ++j) 
            {
                if (u_face_flags(i, j, Flag::Open))
                    u(i, j) += dt * (p_correction(i - 1, j) - p_correction(i, j)) / (dx * fluid.density);
            }
        }
        
        for (int i = 0; i < nx; ++i) 
        {
            for (int j = 1; j < ny; ++j) 
            {
                if (v_face_flags(i, j, Flag::Open))
                    v(i, j) += dt * (p_correction(i, j - 1) - p_correction(i, j)) / (dy * fluid.density);
            }
        }
    }

    void SolveVelocitiesForMomentumEquation() // check density it should be added 
    {
        for (int n = 0; n < innerVelocityItterations; ++n)
        {
            for (int i = 1; i < nx; ++i)
            {
                for (int j = 0; j < ny; ++j)
                {
                    if (u_face_flags(i, j, Flag::Open))
                        u(i, j) = (u(i, j) + dt * (u_scr(i, j) + (p(i - 1, j) - p(i, j)) / (dx * fluid.density))) / (1 + dt * u_coeff(i, j));
                }
            }
            for (int i = 0; i < nx; ++i)
            {
                for (int j = 1; j < ny; ++j)
                {
                    if (v_face_flags(i, j, Flag::Open))
                        v(i, j) = (v(i, j) + dt * (v_scr(i, j) + (p(i, j - 1) - p(i, j)) / (dy * fluid.density))) / (1 + dt * v_coeff(i, j));
                }
            }
        }
    }

    void SolveTemperatures()
    {
        for (int n = 0; n < innerTemperatureItterations; ++n)
        {
            for (int i = 0; i < nx; ++i)
            {
                for (int j = 0; j < ny; ++j)
                {
                    double denominator = (1 + dt * t_coeff(i, j) / fluid.density);
                    t(i, j) += relaxation * (t(i, j) * (1 - denominator) + dt * t_scr(i, j) / fluid.density) / denominator;
                }
            }
        }
    }

    // Main time stepping method
    void solve(int num_iterations) 
    {
        is_solving = true;

        for (int iter = 0; iter < num_iterations; ++iter) 
        {
            if (!is_solving)
                break;

            v_scr.reset();
            u_scr.reset();
            p_correction.reset();
            p_scr.reset();
            p_coef.reset();
            u_coeff.reset();
            v_coeff.reset();
            t_scr.reset();
            t_coeff.reset();
            t_previous = t;
            
            // Solve Temperature

            ApplyTemperatureBoundaryConditions();

            ComputeConvectiveFluxesAndDiffusiveTermsForTemperature();

            SolveTemperatures();

            // Solve Pressure/Velocity

            ApplyVelocityBoundaryConditions();

            ComputeConvectiveTerms();

            ComputeDiffusiveTerms();

            SolveVelocitiesForMomentumEquation();

            ApplyPressureBoundaryConditions();

            ComputeDivergence();

            SolvePressure();

            ApplyPressureCorrectionToVelocity();

            // Check Residuals
            if (iter % 10 == 0)
            {
                residual = p_correction / p;
                double maxPressureResidual = relaxation * (std::max(std::abs(*std::max_element(residual.begin(), residual.end())), std::abs(*std::min_element(residual.begin(), residual.end()))));

                residual = (t - t_previous) / t;
                double maxTemperatureResidual = std::max(std::abs(*std::max_element(residual.begin(), residual.end())), std::abs(*std::min_element(residual.begin(), residual.end())));

                double maxDivergence = std::max(std::abs(*std::ranges::min_element(divergence.begin(), divergence.end())), std::abs(*std::ranges::max_element(divergence.begin(), divergence.end())));

                printf("\rIteration: %d | Pressure Residual: %.3e | Max Abs Divergence %.3e | Temperature Residual: %.3e   ", iter, maxPressureResidual, maxDivergence, maxTemperatureResidual);

                if ((maxPressureResidual < residualLimit) && (maxDivergence < residualLimit) && (maxTemperatureResidual < residualLimit) && (iter > 2))
                    break;
            }

            CalculateCellVelocities();
            interface.SetData(*this);
        }

        is_solving = false;
    }

    void ApplyVelocityBoundaryConditions() 
    {
        for (auto& boundary : inlet_boundary_conditions)
            boundary.ApplyForVelocity(*this);
        
        for (auto& boundary : outlet_boundary_conditions)
            boundary.ApplyForVelocity(*this);

        double wall_velocity = 0.0;
        double coeff = 0.5 * fluid.kinematic_viscosity / (dy / 2);
        for (int i = 0; i < nx + 1; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                if (!u_face_flags(i, j, Flag::Open))
                {
                    if (i != 0)
                    {
                        if (v_face_flags(i - 1, j + 1, Flag::Open)) // Top Left
                        { 
                            v_coeff(i - 1, j + 1) += coeff;
                            v_scr(i - 1, j + 1) += wall_velocity * coeff;
                        }
                        if (v_face_flags(i - 1, j, Flag::Open)) // Bottom Left
                        {
                            v_coeff(i - 1, j) += coeff;
                            v_scr(i - 1, j) += wall_velocity * coeff;
                        }
                    }
                    if (i != nx)
                    {
                        if (v_face_flags(i, j + 1, Flag::Open)) // Top Right
                        {
                            v_coeff(i, j + 1) += coeff;
                            v_scr(i, j + 1) += wall_velocity * coeff;
                        }
                        if (v_face_flags(i, j, Flag::Open)) // Bottom Right
                        {
                            v_coeff(i, j) += coeff;
                            v_scr(i, j) += wall_velocity * coeff;
                        }
                    }
                }
            }
        }
        coeff = 0.5 * fluid.kinematic_viscosity / (dx / 2);
        for (int i = 0; i < nx; ++i)
        {
            for (int j = 0; j < ny + 1; ++j)
            {
                if (!v_face_flags(i, j, Flag::Open))
                {
                    if (j != 0)
                    {
                        if (u_face_flags(i, j - 1, Flag::Open)) // Bottom Left
                        {
                            u_coeff(i, j - 1) += coeff;
                            u_scr(i, j - 1) += wall_velocity * coeff;
                            if (u_face_flags(i + 1, j - 1, Flag::Open)) // Bottom Right
                            {
                                u_coeff(i + 1, j - 1) += coeff;
                                u_scr(i + 1, j - 1) += wall_velocity * coeff;
                            }
                        }
                    }
                    if (j != ny)
                    {
                        if (u_face_flags(i, j, Flag::Open)) // Top Left
                        {
                            u_coeff(i, j) += coeff;
                            u_scr(i, j) += wall_velocity * coeff;
                        }
                        if (u_face_flags(i + 1, j, Flag::Open)) // Top Right
                        {
                            u_coeff(i + 1, j) += coeff;
                            u_scr(i + 1, j) += wall_velocity * coeff;
                        }
                        
                    }
                }
            }
        }
    }

    void ApplyPressureBoundaryConditions()
    {
        for (auto& boundary : open_boundary_condition)
            boundary.ApplyForPressure(*this);
    }

    void ApplyTemperatureBoundaryConditions()
    {
        for (auto& boundary : inlet_boundary_conditions)
            boundary.ApplyForTemperature(*this);

        for (auto& boundary : open_boundary_condition)
            boundary.ApplyForTemperature(*this);
    }

    void CorrectBoundaryCellVelocities()
    {
        for (auto& boundary : inlet_boundary_conditions)
            boundary.CorrectBoundaryCellVelocities(*this);

        for (auto& boundary : outlet_boundary_conditions)
            boundary.CorrectBoundaryCellVelocities(*this);

        for (auto& boundary : open_boundary_condition)
            boundary.CorrectBoundaryCellVelocities(*this);
    }

    //////////////
    // Settings //
    //////////////

    struct FluidProperties
    {
        double density{ 1.0 };
        double kinematic_viscosity{ 0.01 };
        double conductivity{ 0.1 };
        double cp{ 1000.0 };
    };
    FluidProperties fluid;

    int nx, ny;
    double dx, dy;
    double dt{ 10e-3 };
    double relaxation{ 0.5 };

    int innerPressureItterations{ 40 };
    int innerVelocityItterations{ 2 };
    int innerTemperatureItterations{ 20 };
    double residualLimit{ 10e-7 };

    /////////////////
    // Solver Data //
    /////////////////

    // Cell Props

    PhysicalField p{ nx, ny };
    PhysicalField p_correction{ nx, ny };
    PhysicalField p_correction_old{ nx, ny };

    PhysicalField divergence{ nx, ny };
    PhysicalField p_scr{ nx, ny };
    PhysicalField p_coef{ nx, ny };

    PhysicalField t{ nx, ny };
    PhysicalField t_scr{ nx, ny };
    PhysicalField t_coeff{ nx, ny };

    PhysicalField residual{ nx, ny };
    PhysicalField t_previous{ nx, ny };

    // Face Props (staggered grid)

    PhysicalField u{ nx + 1, ny }, v{ nx, ny + 1 };

    PhysicalField u_scr{ nx + 1, ny }; 
    PhysicalField v_scr{ nx, ny + 1 };

    PhysicalField u_coeff{ nx + 1, ny };
    PhysicalField v_coeff{ nx, ny + 1 };

    FairMutex cell_data_mutex;

    PhysicalField cell_u{ nx, ny };
    PhysicalField cell_v{ nx, ny };

    // Face Flags (staggered grid): 1 - blocked, 2 - friction

    FlagField u_face_flags{ nx + 1, ny }; 
    FlagField v_face_flags{ nx, ny + 1 }; 

    ////////////////
    // Boundaries //
    ////////////////

    void SetBlockedFaces(Interface& interface)
    {
        u_face_flags.setFlag(Flag::Open);
        v_face_flags.setFlag(Flag::Open);

        // Boundary faces
        for (int j = 0; j < ny; ++j)
        {
            u_face_flags.clearFlag(0, j, Flag::Open);
            u_face_flags.clearFlag(nx, j, Flag::Open);
        }
        for (int i = 0; i < nx; ++i)
        {
            v_face_flags.clearFlag(i, 0, Flag::Open);
            v_face_flags.clearFlag(i, ny, Flag::Open);
        }

        // User set faces
        for (auto& face : interface.blocked_faces)
        {
            if (face.component_dir == 0) u_face_flags.clearFlag(face.i, face.j, Flag::Open);
            else v_face_flags.clearFlag(face.i, face.j, Flag::Open);
        }
    }

    void SetBoundaries(Interface& interface)
    {
        // User set 
        for (auto& boundary : interface.boundary_faces)
        {
            if (boundary.type == Boundary::FixedInflow)
            {
                auto& bc = inlet_boundary_conditions.emplace_back(boundary.boundary_dir, boundary.component_dir);
                bc.inlet_temperature = boundary.optional_temp;
                bc.inlet_velocity = boundary.optional_velocity;
                bc.boundary_faces.push_back({ boundary.i, boundary.j });
                boundary.component_dir == 0 ? u_face_flags.clearFlag(boundary.i, boundary.j, Flag::Open) : v_face_flags.clearFlag(boundary.i, boundary.j, Flag::Open);
            }
            if (boundary.type == Boundary::FixedOutflow)
            {
                auto& bc = outlet_boundary_conditions.emplace_back(boundary.boundary_dir, boundary.component_dir);
                bc.outlet_velocity = boundary.optional_velocity;
                bc.boundary_faces.push_back({ boundary.i, boundary.j });
                boundary.component_dir == 0 ? u_face_flags.clearFlag(boundary.i, boundary.j, Flag::Open) : v_face_flags.clearFlag(boundary.i, boundary.j, Flag::Open);
            }
            else if (boundary.type == Boundary::Open)
            {
                if ((boundary.component_dir == 0 && boundary.i != nx) || (boundary.component_dir == 1 && boundary.j != ny))
                {
                    auto& bc = open_boundary_condition.emplace_back(0, boundary.component_dir);
                    bc.temperature = boundary.optional_temp;
                    bc.pressure = boundary.optional_pressure;
                    bc.boundary_faces.push_back({ boundary.i, boundary.j });
                    boundary.component_dir == 0 ? u_face_flags.clearFlag(boundary.i, boundary.j, Flag::Open) : v_face_flags.clearFlag(boundary.i, boundary.j, Flag::Open);
                }
                if ((boundary.component_dir == 0 && boundary.i != 0) || (boundary.component_dir == 1 && boundary.j != 0))
                {
                    auto& bc = open_boundary_condition.emplace_back(1, boundary.component_dir);
                    bc.temperature = boundary.optional_temp;
                    bc.pressure = boundary.optional_pressure;
                    bc.boundary_faces.push_back({ boundary.i, boundary.j });
                    boundary.component_dir == 0 ? u_face_flags.clearFlag(boundary.i, boundary.j, Flag::Open) : v_face_flags.clearFlag(boundary.i, boundary.j, Flag::Open);
                }
            }
        }
    }

    std::vector<FixedInletBoundaryCondition> inlet_boundary_conditions;
    std::vector<FixedOutletBoundaryCondition> outlet_boundary_conditions;
    std::vector<OpenBoundaryCondition> open_boundary_condition;
    std::vector<FrictionBoundaryCondition> friction_boundary_condition;

    // TODO
    // - blocked cells / solid cells
    // - turbulance
    // - gravity force

    // Coding ToDo
    // could the api for faces be better?

};