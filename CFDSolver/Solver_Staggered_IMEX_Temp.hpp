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

class SolverStaggeredIMEXTemp
{
public:
    SolverStaggeredIMEXTemp(int nx, int ny, double dx, double dy) : nx(nx), ny(ny), dx(dx), dy(dy)
    {
        //if (dt > 0.5 * std::pow(std::min(dx, dy), 2) / fluid.kinematic_viscosity)
        //    std::cout << "Solution Will Be Unstable\n";

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

    void ComputeConvectiveFluxes()// check density it should be added 
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
                        coeff = flux / (dx * dy);
                        u_coeff(i, j) += coeff;
                        u_scr(i, j) += coeff * value;
                    }

                    // Bottom Right: u flux
                    flux = (dx / 2.0) * v(i, j);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? u(i, j - 1) : u(i, j);
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
                        coeff = flux / (dx * dy);
                        u_coeff(i, j) += coeff;
                        u_scr(i, j) += coeff * value;
                    }

                    // Top Right: u flux
                    flux = -(dx / 2.0) * v(i, j + 1);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? u(i, j + 1) : u(i, j);
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
                        coeff = flux / (dx * dy);
                        v_coeff(i, j) += coeff;
                        v_scr(i, j) += coeff * value;
                    }

                    // Bottom Right: v flux
                    flux = (dy / 2.0) * u(i, j);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? v(i - 1, j) : v(i, j);
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
                        coeff = flux / (dx * dy);
                        v_coeff(i, j) += coeff;
                        v_scr(i, j) += coeff * value;
                    }

                    // Top Right: v flux
                    flux = -(dy / 2.0) * u(i + 1, j);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? v(i + 1, j) : v(i, j);
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
                    u_coeff(i, j) += coeff;
                    u_scr(i, j) += coeff * value;
                }

                if (u_face_flags(i + 1, j) != 1)
                {
                    // Inline Right: u diff
                    coeff = fluid.kinematic_viscosity / (dx * dx);
                    value = u(i + 1, j);
                    u_coeff(i, j) += coeff;
                    u_scr(i, j) += coeff * value;
                }

                if (j != 0 && v_face_flags(i, j) != 1 && v_face_flags(i - 1, j) != 1 && u_face_flags(i, j - 1) != 1)
                {
                    // Bottom: u diff
                    coeff = fluid.kinematic_viscosity / (dy * dy);
                    value = u(i, j - 1);
                    u_coeff(i, j) += coeff;
                    u_scr(i, j) += coeff * value;
                }

                if (j != ny - 1 && v_face_flags(i, j + 1) != 1 && v_face_flags(i - 1, j + 1) != 1 && u_face_flags(i, j + 1) != 1)
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
                if (v_face_flags(i, j) == 1)
                    continue;

                if (v_face_flags(i, j - 1) != 1)
                {
                    // Inline Left: v diff
                    coeff = fluid.kinematic_viscosity / (dy * dy);
                    value = v(i, j - 1);
                    v_coeff(i, j) += coeff;
                    v_scr(i, j) += coeff * value;
                }

                if (v_face_flags(i, j + 1) != 1)
                {
                    // Inline Right: v diff
                    coeff = fluid.kinematic_viscosity / (dy * dy);
                    value = v(i, j + 1);
                    v_coeff(i, j) += coeff;
                    v_scr(i, j) += coeff * value;
                }

                if (i != 0 && u_face_flags(i + 1, j) + u_face_flags(i + 1, j - 1) + v_face_flags(i + 1, j))
                {
                    // Bottom: v diff
                    coeff = fluid.kinematic_viscosity / (dx * dx);
                    value = v(i - 1, j);
                    v_coeff(i, j) += coeff;
                    v_scr(i, j) += coeff * value;
                }

                if (i != nx - 1 && u_face_flags(i, j) + u_face_flags(i, j - 1) + v_face_flags(i - 1, j))
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

    void ComputeConvectiveFluxesAndDiffusiveTerms()
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
                        coeff = flux / (dx * dy);
                        u_coeff(i, j) += coeff;
                        u_scr(i, j) += coeff * value;
                    }

                    // Inline Left: u diff
                    coeff = fluid.kinematic_viscosity / (dx * dx);
                    value = u(i - 1, j);
                    u_coeff(i, j) += coeff;
                    u_scr(i, j) += coeff * value;
                }

                if (u_face_flags(i + 1, j) != 1)
                {
                    // Inline Right: u flux
                    flux = dy * -((u(i, j) + u(i + 1, j)) / 2.0);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? u(i + 1, j) : u(i, j);
                        coeff = flux / (dx * dy);
                        u_coeff(i, j) += coeff;
                        u_scr(i, j) += coeff * value;
                    }

                    // Inline Right: u diff
                    coeff = fluid.kinematic_viscosity / (dx * dx);
                    value = u(i + 1, j);
                    u_coeff(i, j) += coeff;
                    u_scr(i, j) += coeff * value;
                }

                if (j != 0 && v_face_flags(i, j) != 1 && v_face_flags(i - 1, j) != 1 && u_face_flags(i, j - 1) != 1)
                {
                    // Bottom Left: u flux
                    flux = (dx / 2.0) * v(i - 1, j);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? u(i, j - 1) : u(i, j);
                        coeff = flux / (dx * dy);
                        u_coeff(i, j) += coeff;
                        u_scr(i, j) += coeff * value;
                    }

                    // Bottom Right: u flux
                    flux = (dx / 2.0) * v(i, j);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? u(i, j - 1) : u(i, j);
                        coeff = flux / (dx * dy);
                        u_coeff(i, j) += coeff;
                        u_scr(i, j) += coeff * value;
                    }

                    // Bottom: u diff
                    coeff = fluid.kinematic_viscosity / (dy * dy);
                    value = u(i, j - 1);
                    u_coeff(i, j) += coeff;
                    u_scr(i, j) += coeff * value;
                }

                if (j != ny - 1 && v_face_flags(i, j + 1) != 1 && v_face_flags(i - 1, j + 1) != 1 && u_face_flags(i, j + 1) != 1)
                {
                    // Top Left: u flux
                    flux = -(dx / 2.0) * v(i - 1, j + 1);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? u(i, j + 1) : u(i, j);
                        coeff = flux / (dx * dy);
                        u_coeff(i, j) += coeff;
                        u_scr(i, j) += coeff * value;
                    }

                    // Top Right: u flux
                    flux = -(dx / 2.0) * v(i, j + 1);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? u(i, j + 1) : u(i, j);
                        coeff = flux / (dx * dy);
                        u_coeff(i, j) += coeff;
                        u_scr(i, j) += coeff * value;
                    }

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
                if (v_face_flags(i, j) == 1)
                    continue;

                if (v_face_flags(i, j - 1) != 1)
                {
                    // Inline Left: v flux
                    flux = dy * ((v(i, j) + v(i, j - 1)) / 2.0);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? v(i, j - 1) : v(i, j);
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
                        coeff = flux / (dx * dy);
                        v_coeff(i, j) += coeff;
                        v_scr(i, j) += coeff * value;
                    }

                    // Bottom Right: v flux
                    flux = (dy / 2.0) * u(i, j);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? v(i - 1, j) : v(i, j);
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
                        coeff = flux / (dx * dy);
                        v_coeff(i, j) += coeff;
                        v_scr(i, j) += coeff * value;
                    }

                    // Top Right: v flux
                    flux = -(dy / 2.0) * u(i + 1, j);
                    if (flux > 0.0)
                    {
                        value = flux >= 0.0 ? v(i + 1, j) : v(i, j);
                        coeff = flux / (dx * dy);
                        v_coeff(i, j) += coeff;
                        v_scr(i, j) += coeff * value;
                    }
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
                if (u_face_flags(i, j) != 1)
                {
                    value = t(i - 1, j);// -t(i, j);

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

                if (u_face_flags(i + 1, j) != 1)
                {
                    value = t(i + 1, j);// -t(i, j);

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

                if (v_face_flags(i, j) != 1)
                {
                    value = t(i, j - 1);// -t(i, j);

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

                if (v_face_flags(i, j + 1) != 1)
                {
                    value = t(i, j + 1);// -t(i, j);

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

                    if (u_face_flags(i, j) != 1)
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

    void SolveVelocitiesForMomentumEquation() // check density it should be added 
    {
        for (int n = 0; n < innerVelocityItterations; ++n)
        {
            for (int i = 1; i < nx; ++i)
            {
                for (int j = 0; j < ny; ++j)
                {
                    if (u_face_flags(i, j) != 1)
                        u(i, j) = (u(i, j) + dt * (u_scr(i, j) + (p(i - 1, j) - p(i, j)) / (dy * fluid.density))) / (1 + dt * u_coeff(i, j));
                }
            }
            for (int i = 0; i < nx; ++i)
            {
                for (int j = 1; j < ny; ++j)
                {
                    if (v_face_flags(i, j) != 1)
                        v(i, j) = (v(i, j) + dt * (v_scr(i, j) + (p(i, j - 1) - p(i, j)) / (dx * fluid.density))) / (1 + dt * v_coeff(i, j));
                }
            }
        }
    }

    void SolveTemperatures()
    {
        for (int n = 0; n < 5; ++n)
        {
            for (int i = 0; i < nx; ++i)
            {
                for (int j = 0; j < ny; ++j)
                {
                    t(i, j) = (t(i, j) + dt * t_scr(i, j) / fluid.density) / (1 + dt * t_coeff(i, j) / fluid.density);
                }
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
            t_scr.reset();
            t_coeff.reset();
            
            // Solve Temperature

            ApplyTemperatureBoundaryConditions();

            ComputeConvectiveFluxesAndDiffusiveTermsForTemperature();

            SolveTemperatures();

            // Solve Pressure/Velocity

            ApplyVelocityBoundaryConditions();

            ComputeConvectiveFluxes();

            ComputeDiffusiveTerms();

            SolveVelocitiesForMomentumEquation();

            ComputeDivergence();

            ApplyPressureBoundaryConditions();

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

        auto plotter = CFDVisualizer(nx, ny, (float)dx, (float)dy, cell_u, cell_v, u, v, p, t, divergence, cell_data_mutex);
        plotter.Render();
    }

    void ApplyVelocityBoundaryConditions() 
    {
        double coeff = 0.0;

        // Inflow top left
        double inflowVelocity = 1.0;
        double flux = inflowVelocity * dy;
        coeff = flux / (dx * dy);
        for (int j = 0; j < ny/4; ++j)
        {
            u(0, j) = inflowVelocity;
            u_coeff(1, j) += coeff;
            u_scr(1, j) += coeff * inflowVelocity;

            //u(nx, j) = inflowVelocity;
        }

        // Inflow top
        inflowVelocity = 1.0;
        flux = inflowVelocity * dx;
        coeff = flux / (dx * dy);
        for (int i = ny / 3; i < ny / 2; ++i)
        {
            v(i, 0) = inflowVelocity;
            v_coeff(i, 1) += coeff;
            v_scr(i, 1) += coeff * inflowVelocity;
        }
        
        //for (int j = 2 * ny / 4; j < 3 * ny / 4; ++j)
        //{
        //    u(0, j) = inflowVelocity;
        //    u_coeff(1, j) += coeff;
        //    u_scr(1, j) += coeff * inflowVelocity;
        //}

        //v(25, 0) = inflowVelocity;
        //v_coeff(25, 1) = coeff;
        //v_scr(25, 1) = coeff * inflowVelocity;

        // Friction - ToDo apply from each blocked face take half for completeness - TODO blocked faces
        //coeff = fluid.kinematic_viscosity / (dy / 2);
        //for (int i = 1; i < nx; ++i)
        //{
        //    // Top boundary
        //    u_coeff(i, 0) += coeff;
        //    u_scr(i, 0) += 0.0 * coeff;
        //
        //    // Bottom boundary
        //    u_coeff(i, ny - 1) += coeff;
        //    u_scr(i, ny - 1) += 0.0 * coeff;
        //}
        //coeff = fluid.kinematic_viscosity / (dx / 2);
        //for (int j = 1; j < ny; ++j)
        //{
        //    // Left boundary
        //    v_coeff(0, j) += coeff;
        //    v_scr(0, j) += 0.0 * coeff;
        //
        //    // Right boundary
        //    v_coeff(nx - 1, j) += coeff;
        //    v_scr(nx - 1, j) += 0.0 * coeff;
        //}
    }

    void ApplyPressureBoundaryConditions()
    {
        // Top right outflow
        for (int j = 0; j < ny / 4; ++j)
        {
            divergence(nx - 1, j) = 0.0;
        }
    }

    void ApplyTemperatureBoundaryConditions()
    {
        // Inflow
        for (int j = 0; j < ny / 4; ++j)
        {
            double inflowTemp = 20.0;

            double flux = u(0,j) * fluid.density * dy;
            double coeff = flux / (dx * dy);

            t_coeff(0, j) += coeff;
            t_scr(0, j) += coeff * inflowTemp;
        }
        
        for (int i = ny / 3; i < ny / 2; ++i)
        {
            double inflowTemp = 40.0;

            double flux = v(i, 0) * fluid.density * dy;
            double coeff = flux / (dx * dy);

            t_coeff(i, 0) += coeff;
            t_scr(i, 0) += coeff * inflowTemp;
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
        for (int j = 0; j < ny / 2; ++j)
        {
            u_face_flags(nx / 2, j) = 1;
        }
        for (int i = nx / 4; i < nx / 2; ++i)
        {
            v_face_flags(i, ny / 2) = 1;
        }
    }

    struct FluidProperties
    {
        double density{ 1.0 };
        double kinematic_viscosity{ 0.01 };
        double conductivity{ 0.01 };
        double cp{ 1000.0 };
    };
    FluidProperties fluid;

    int nx, ny;
    double dx, dy;
    double dt{ 10e-4 };
    double relaxation{ 0.5 };

    int innerPressureItterations{ 40 };
    int innerVelocityItterations{ 2 };
    double residualLimit{ 10e-5 };
    double divergenceLimit{ 10e-5 };

    // Cell prop

    PhysicalField p{ nx, ny };
    PhysicalField p_correction{ nx, ny };
    PhysicalField p_correction_old{ nx, ny };

    PhysicalField divergence{ nx, ny };

    PhysicalField t{ nx, ny };
    PhysicalField t_scr{ nx, ny };
    PhysicalField t_coeff{ nx, ny };

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