#pragma once

#include <vector>
#include <iostream>
#include <cmath>
#include <cmath>
#include <limits>
#include <thread>
#include <algorithm>
#include "FairMutex.hpp"

#include "PhysicalField.hpp"
#include "Plotter.hpp"

class SolverStaggeredExplicit
{
public:
    SolverStaggeredExplicit(int nx, int ny, double dx, double dy) : nx(nx), ny(ny), dx(dx), dy(dy)
    {
        if (dt > 0.5 * std::pow(std::min(dx, dy), 2) / fluid.kinematic_viscosity)
            std::cout << "Solution Will Be Unstable\n";

        SetBlockedFaces();
    }

    void ComputeDivergence() 
    {
        for (int i = 0; i < nx; ++i) 
        {
            for (int j = 0; j < ny; ++j) 
            {
                divergence(i, j) = ((u(i + 1, j) - u(i, j)) * dy + 
                                    (v(i, j + 1) - v(i, j)) * dx);
            }
        }
    }

    void ComputeConvectiveFluxes()
    {
        double flux = 0.0;
        double value = 0.0;
        double coeff = 0.0;

        for (int i = 1; i < nx; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                if (u_face_flags(i, j) == 1)
                    continue;

                if (u_face_flags(i - 1, j) != 1)
                {
                    // Inline Left: u flux
                    flux = dy * ((u(i, j) + u(i - 1, j)) / 2.0);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? u(i - 1, j) : u(i, j);
                        value = value - u(i, j);
                        coeff = flux / (dx * dy);
                        u_coeff(i, j) += coeff;
                        u_scr(i, j) += coeff * value;
                    }
                }

                if (u_face_flags(i + 1, j) != 1)
                {
                    // Inline Right: u flux
                    flux = dy * -((u(i, j) + u(i + 1, j)) / 2.0);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? u(i + 1, j) : u(i, j);
                        value = value - u(i, j);
                        coeff = flux / (dx * dy);
                        u_coeff(i, j) += coeff;
                        u_scr(i, j) += coeff * value;
                    }
                }

                if (j != 0 && v_face_flags(i, j) != 1 && v_face_flags(i - 1, j) != 1 && u_face_flags(i, j - 1) != 1)
                {
                    // Bottom Left: u flux
                    flux = (dx / 2.0) * v(i - 1, j);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? u(i, j - 1) : u(i, j);
                        value = value - u(i, j);
                        coeff = flux / (dx * dy);
                        u_coeff(i, j) += coeff;
                        u_scr(i, j) += coeff * value;
                    }

                    // Bottom Right: u flux
                    flux = (dx / 2.0) * v(i, j);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? u(i, j - 1) : u(i, j);
                        value = value - u(i, j);
                        coeff = flux / (dx * dy);
                        u_coeff(i, j) += coeff;
                        u_scr(i, j) += coeff * value;
                    }
                }

                if (j != ny - 1 && v_face_flags(i, j + 1) != 1 && v_face_flags(i - 1, j + 1) != 1 && u_face_flags(i, j + 1) != 1)
                {
                    // Top Left: u flux
                    flux = -(dx / 2.0) * v(i - 1, j + 1);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? u(i, j + 1) : u(i, j);
                        value = value - u(i, j);
                        coeff = flux / (dx * dy);
                        u_coeff(i, j) += coeff;
                        u_scr(i, j) += coeff * value;
                    }

                    // Top Right: u flux
                    flux = -(dx / 2.0) * v(i, j + 1);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? u(i, j + 1) : u(i, j);
                        value = value - u(i, j);
                        coeff = flux / (dx * dy);
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
                if (v_face_flags(i, j) == 1)
                    continue;

                if (v_face_flags(i, j - 1) != 1)
                {
                    // Inline Left: v flux
                    flux = dy * ((v(i, j) + v(i, j - 1)) / 2.0);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? v(i, j - 1) : v(i, j);
                        value = value - v(i, j);
                        coeff = flux / (dx * dy);
                        v_coeff(i, j) += coeff;
                        v_scr(i, j) += coeff * value;
                    }
                }

                if (v_face_flags(i, j + 1) != 1)
                {
                    // Inline Right: v flux
                    flux = dy * -((v(i, j) + v(i, j + 1)) / 2.0);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? v(i, j + 1) : v(i, j);
                        value = value - v(i, j);
                        coeff = flux / (dx * dy);
                        v_coeff(i, j) += coeff;
                        v_scr(i, j) += coeff * value;
                    }
                }

                if (i != 0 && u_face_flags(i + 1, j) + u_face_flags(i + 1, j - 1) + v_face_flags(i + 1, j))
                {
                    // Bottom Left: v flux
                    flux = (dy / 2.0) * u(i, j - 1);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? v(i - 1, j) : v(i, j);
                        value = value - v(i, j);
                        coeff = flux / (dx * dy);
                        v_coeff(i, j) += coeff;
                        v_scr(i, j) += coeff * value;
                    }

                    // Bottom Right: v flux
                    flux = (dy / 2.0) * u(i, j);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? v(i - 1, j) : v(i, j);
                        value = value - v(i, j);
                        coeff = flux / (dx * dy);
                        v_coeff(i, j) += coeff;
                        v_scr(i, j) += coeff * value;
                    }
                }

                if (i != nx - 1 && u_face_flags(i, j) + u_face_flags(i, j - 1) + v_face_flags(i - 1, j))
                {
                    // Top Left: v flux
                    flux = -(dy / 2.0) * u(i + 1, j - 1);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? v(i + 1, j) : v(i, j);
                        value = value - v(i, j);
                        coeff = flux / (dx * dy);
                        v_coeff(i, j) += coeff;
                        v_scr(i, j) += coeff * value;
                    }

                    // Top Right: v flux
                    flux = -(dy / 2.0) * u(i + 1, j);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? v(i + 1, j) : v(i, j);
                        value = value - v(i, j);
                        coeff = flux / (dx * dy);
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
                if (u_face_flags(i, j) == 1)
                    continue;

                if (u_face_flags(i - 1, j) != 1)
                {
                    // Inline Left: u diff
                    coeff = fluid.kinematic_viscosity / (dx * dx);
                    value = u(i - 1, j);
                    value = (value - u(i, j));
                    u_coeff(i, j) += coeff;
                    u_scr(i, j) += coeff * value;
                }

                if (u_face_flags(i + 1, j) != 1)
                {
                    // Inline Right: u diff
                    coeff = fluid.kinematic_viscosity / (dx * dx);
                    value = u(i + 1, j);
                    value = (value - u(i, j));
                    u_coeff(i, j) += coeff;
                    u_scr(i, j) += coeff * value;
                }

                if (j != 0 && v_face_flags(i, j) != 1 && v_face_flags(i - 1, j) != 1 && u_face_flags(i, j - 1) != 1)
                {
                    // Bottom: u diff
                    coeff = fluid.kinematic_viscosity / (dy * dy);
                    value = u(i, j - 1);
                    value = (value - u(i, j));
                    u_coeff(i, j) += coeff;
                    u_scr(i, j) += coeff * value;
                }

                if (j != ny - 1 && v_face_flags(i, j + 1) != 1 && v_face_flags(i - 1, j + 1) != 1 && u_face_flags(i, j + 1) != 1)
                {
                    // Top: u diff
                    coeff = fluid.kinematic_viscosity / (dy * dy);
                    value = u(i, j + 1);
                    value = (value - u(i, j));
                    u_coeff(i, j) += coeff;
                    u_scr(i, j) += coeff * value;
                }
            }
        }

        for (int i = 0; i < nx; ++i)
        {
            for (int j = 1; j < ny; ++j)
            {
                if (v_face_flags(i, j) == 1)
                    continue;

                if (v_face_flags(i, j - 1) != 1)
                {
                    // Inline Left: v diff
                    coeff = fluid.kinematic_viscosity / (dy * dy);
                    value = v(i, j - 1);
                    value = (value - v(i, j));
                    v_coeff(i, j) += coeff;
                    v_scr(i, j) += coeff * value;
                }

                if (v_face_flags(i, j + 1) != 1)
                {
                    // Inline Right: v diff
                    coeff = fluid.kinematic_viscosity / (dy * dy);
                    value = v(i, j + 1);
                    value = (value - v(i, j));
                    v_coeff(i, j) += coeff;
                    v_scr(i, j) += coeff * value;
                }

                if (i != 0 && u_face_flags(i + 1, j) + u_face_flags(i + 1, j - 1) + v_face_flags(i + 1, j))
                {
                    // Bottom: v diff
                    coeff = fluid.kinematic_viscosity / (dx * dx);
                    value = v(i - 1, j);
                    value = (value - v(i, j));
                    v_coeff(i, j) += coeff;
                    v_scr(i, j) += coeff * value;
                }

                if (i != nx - 1 && u_face_flags(i, j) + u_face_flags(i, j - 1) + v_face_flags(i - 1, j))
                {
                    // Top: v diff
                    coeff = fluid.kinematic_viscosity/ (dx * dx);
                    value = v(i + 1, j);
                    value = (value - v(i, j));
                    v_coeff(i, j) += coeff;
                    v_scr(i, j) += coeff * value;
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

                    double Ap = u_coeff(i, j) + u_coeff(i + 1, j) + v_coeff(i, j) + v_coeff(i, j + 1);

                    if (Ap == 0)
                        continue;

                    if (u_face_flags(i,j) != 1)
                        p_correction(i, j) += (u_coeff(i, j) / Ap) * p_correction_old(i - 1, j);
                    if (u_face_flags(i + 1, j) != 1)
                        p_correction(i, j) += (u_coeff(i + 1, j) / Ap) * p_correction_old(i + 1, j);
                    if (v_face_flags(i, j) != 1)
                        p_correction(i, j) += (v_coeff(i, j) / Ap) * p_correction_old(i, j - 1);
                    if (v_face_flags(i, j + 1) != 1)
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
                if (u_face_flags(i, j) != 1)
                    u(i, j) += dt * (p_correction(i - 1, j) - p_correction(i, j)) / (dx * fluid.density);
            }
        }
        
        for (int i = 0; i < nx; ++i) 
        {
            for (int j = 1; j < ny; ++j) 
            {
                if (v_face_flags(i, j) != 1)
                    v(i, j) += dt * (p_correction(i, j - 1) - p_correction(i, j)) / (dy * fluid.density);
            }
        }
    }

    void SolveVelocitiesForMomentumEquation()
    {
        for (int i = 1; i < nx; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                if (u_face_flags(i, j) != 1)
                    u(i, j) += dt * (u_scr(i, j) - u(i, j) * u_coeff(i, j) + (p(i - 1, j) - p(i, j)) / (dy * fluid.density));
            }
        }
        for (int i = 0; i < nx; ++i)
        {
            for (int j = 1; j < ny; ++j)
            {
                if (v_face_flags(i, j) != 1)
                    v(i, j) += dt * (v_scr(i, j) - v(i, j) * v_coeff(i, j) + (p(i, j - 1) - p(i, j)) / (dx * fluid.density));
            }
        }
    }

    // Main time stepping method
    void solve(int num_iterations) 
    {
        for (int iter = 0; iter < num_iterations; ++iter) 
        {
            v_scr.reset();
            u_scr.reset();
            p_correction.reset();
            u_coeff.reset();
            v_coeff.reset();
            
            ApplyVelocityBoundaryConditions();

            ComputeConvectiveFluxes();

            ComputeDiffusiveTerms();

            SolveVelocitiesForMomentumEquation();

            ComputeDivergence();

            SolvePressure();

            ApplyPressureCorrectionToVelocity();

            // Check Residuals
            if (iter % 5 == 0)
            {
                double maxResidual{ 0.0 };
                for (int i = 0; i < p.values.size(); ++i)
                {
                    double residual = std::abs((relaxation * p_correction.values[i]) / (p.values[i] + 1e-20));
                    maxResidual = std::max(maxResidual, residual);
                }
                double maxDivergence = std::max(std::abs(*std::ranges::min_element(divergence.values)), std::abs(*std::ranges::max_element(divergence.values)));

                printf("\rIteration: %d | Pressure Residual: %.3e | Max Abs Divergence %.3e   ", iter, maxResidual, maxDivergence);
                if ((maxResidual < residualLimit) && (maxDivergence < divergenceLimit) && (iter > 2))
                    break;
            }
        }

        CalculateCellVelocities();

        auto plotter = CFDVisualizer(nx, ny, (float)dx, (float)dy, cell_u, cell_v, u, v, p, divergence, cell_data_mutex);
        plotter.Render();
    }

    void ApplyVelocityBoundaryConditions() 
    {
        double coeff = 0.0;

        // Inflow
        double inflowVelocity = 1.0;
        double flux = inflowVelocity * dy;
        coeff = flux / (dx * dy);
        for (int j = 0; j < ny/4; ++j)
        {
            u(0, j) = inflowVelocity;
            u_coeff(1, j) += coeff;
            u_scr(1, j) += coeff * (inflowVelocity - u(1, j));
            u_scr(1, j) += coeff * inflowVelocity;
            
            u(nx, j) = inflowVelocity;
        }

        // Friction - ToDo apply from each blocked face take half for completeness
        coeff = fluid.kinematic_viscosity / (dy / 2);
        for (int i = 1; i < nx; ++i)
        {
            // Top boundary
            u_coeff(i, 0) += coeff;
            u_scr(i, 0) += (0.0 - u(i, 0)) * coeff;
            u_scr(i, 0) += 0.0 * coeff;
        
            // Bottom boundary
            u_coeff(i, ny - 1) += coeff;
            u_scr(i, ny - 1) += (0.0 - u(i, ny - 1)) * coeff;
            u_scr(i, ny - 1) += 0.0 * coeff;
        }
        coeff = fluid.kinematic_viscosity / (dx / 2);
        for (int j = 1; j < ny; ++j)
        {
            // Left boundary
            v_coeff(0, j) += coeff;
            v_scr(0, j) += (0.0 - v(0, j)) * coeff;
            v_scr(0, j) += 0.0 * coeff;
        
            // Right boundary
            v_coeff(nx - 1, j) += coeff;
            v_scr(nx - 1, j) += (0.0 - v(nx - 1, j)) * coeff;
            v_scr(nx - 1, j) += 0.0 * coeff;
        }
    }

    void SetBlockedFaces()
    {
        // Boundary faces
        for (int j = 0; j < ny; ++j)
        {
            u_face_flags(0, j) = 1;
            u_face_flags(nx, j) = 1;
        }
        for (int i = 0; i < nx; ++i)
        {
            v_face_flags(i, 0) = 1;
            v_face_flags(i, ny) = 1;
        }

        // Inner blocked faces
        for (int j = 0; j < nx / 2; ++j)
        {
            u_face_flags(nx / 3, j) = 1;
        }
    }

    struct FluidProperties
    {
        double density{ 1.0 };
        double kinematic_viscosity{ 0.1 };
    };
    FluidProperties fluid;

    int nx, ny;
    double dx, dy;
    double dt{ 10e-3 };
    double relaxation{ 0.5 };

    int innerPressureItterations{ 40 };
    double residualLimit{ 10e-8 };
    double divergenceLimit{ 10e-10 };

    // Cell prop

    PhysicalField p{ nx, ny };
    PhysicalField p_correction{ nx, ny };
    PhysicalField p_correction_old{ nx, ny };

    PhysicalField divergence{ nx, ny };

    PhysicalField t{ nx, ny };

    // Face Prop

    PhysicalField u{ nx + 1, ny }, v{ nx, ny + 1 }; // staggered grid

    PhysicalField u_scr{ nx + 1, ny }; // staggered grid
    PhysicalField v_scr{ nx, ny + 1 }; // staggered grid

    PhysicalField u_coeff{ nx + 1, ny }; // staggered grid
    PhysicalField v_coeff{ nx, ny + 1 }; // staggered grid

    PhysicalField u_face_flags{ nx + 1, ny }; // staggered grid
    PhysicalField v_face_flags{ nx, ny + 1 }; // staggered grid

    FairMutex cell_data_mutex;

    PhysicalField cell_u{ nx, ny };
    PhysicalField cell_v{ nx, ny };

    // TODO
    // - blocked cells / solid cells
    // - turbulance
    // - temperature
    // - gravity force

    // Coding ToDo
    // template PhysicalField better to make a bool field
    // could the api for faces be better?

};