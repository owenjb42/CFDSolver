#pragma once

#include <vector>
#include <iostream>
#include <cmath>
#include <cmath>
#include <limits>
#include <thread>
#include "FairMutex.hpp"

#include "PhysicalField.hpp"
#include "Plotter.hpp"

class SolverStaggered
{
public:
    SolverStaggered(int nx, int ny, double dx, double dy) : nx(nx), ny(ny), dx(dx), dy(dy)
    {
        if (dt > 0.5 * std::pow(std::min(dx, dy), 2) / fluid.kinematic_viscosity)
            std::cout << "Unstable\n";
    }

    void ComputeDivergence() 
    {
        for (int i = 0; i < nx; ++i) 
        {
            for (int j = 0; j < ny; ++j) 
            {
                divergence(i, j) = ((u(i + 1, j) - u(i, j)) * dy + 
                                    (v(i, j + 1) - v(i, j)) * dx)
                                    / (dx * dy);
            }
        }
    }

    void ComputeConvectiveFluxes()
    {
        double flux = 0.0;
        double value = 0.0;

        for (int i = 1; i < nx; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                // Inline Left: u flux
                flux = dy * ((u(i, j) + u(i - 1, j)) / 2.0);
                value = flux >= 0.0 ? u(i - 1, j) : u(i, j);
                u_flux(i, j) += flux * value / (dx * dy);

                // Inline Right: u flux
                flux = dy * -((u(i, j) + u(i + 1, j)) / 2.0);
                value = flux >= 0.0 ? u(i + 1, j) : u(i, j);
                u_flux(i, j) += flux * value / (dx * dy);

                if (j != 0)
                {
                    // Bottom Left: u flux
                    flux = (dx / 2.0) * v(i - 1, j);
                    value = flux >= 0.0 ? u(i, j - 1) : u(i, j);
                    u_flux(i, j) += flux * value / (dx * dy);

                    // Bottom Right: u flux
                    flux = (dx / 2.0) * v(i, j);
                    value = flux >= 0.0 ? u(i, j - 1) : u(i, j);
                    u_flux(i, j) += flux * value / (dx * dy);
                }

                if (j != ny - 1)
                {
                    // Top Left: u flux
                    flux = -(dx / 2.0) * v(i - 1, j + 1);
                    value = flux >= 0.0 ? u(i, j + 1) : u(i, j);
                    u_flux(i, j) += flux * value / (dx * dy);

                    // Top Right: u flux
                    flux = -(dx / 2.0) * v(i, j + 1);
                    value = flux >= 0.0 ? u(i, j + 1) : u(i, j);
                    u_flux(i, j) += flux * value / (dx * dy);
                }
            }
        }

        for (int i = 0; i < nx; ++i)
        {
            for (int j = 1; j < ny; ++j)
            {
                // Inline Left: v flux
                flux = dy * ((v(i, j) + v(i, j - 1)) / 2.0);
                value = flux >= 0.0 ? v(i, j - 1) : v(i, j);
                v_flux(i, j) += flux * value / (dx * dy);

                // Inline Right: v flux
                flux = dy * -((v(i, j) + v(i, j + 1)) / 2.0);
                value = flux >= 0.0 ? v(i, j + 1) : v(i, j);
                v_flux(i, j) += flux * value / (dx * dy);

                if (i != 0)
                {
                    // Bottom Left: v flux
                    flux = (dy / 2.0) * u(i, j - 1);
                    value = flux >= 0.0 ? v(i - 1, j) : v(i, j);
                    v_flux(i, j) += flux * value / (dx * dy);

                    // Bottom Right: v flux
                    flux = (dy / 2.0) * u(i, j);
                    value = flux >= 0.0 ? v(i - 1, j) : v(i, j);
                    v_flux(i, j) += flux * value / (dx * dy);
                }

                if (i != nx - 1)
                {
                    // Top Left: v flux
                    flux = -(dy / 2.0) * u(i + 1, j - 1);
                    value = flux >= 0.0 ? v(i + 1, j) : v(i, j);
                    v_flux(i, j) += flux * value / (dx * dy);

                    // Top Right: v flux
                    flux = -(dy / 2.0) * u(i + 1, j);
                    value = flux >= 0.0 ? v(i + 1, j) : v(i, j);
                    v_flux(i, j) += flux * value / (dx * dy);
                }
            }
        }
    }

    // Compute diffusive terms (on faces)
    void ComputeDiffusiveTerms() 
    {
        for (int i = 1; i < nx; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                // Inline Left: u diff
                u_dif(i, j) += fluid.kinematic_viscosity * (
                    dy * (u(i - 1, j) - u(i, j)) / dx
                    ) / (dx * dy);

                // Inline Right: u diff
                u_dif(i, j) += fluid.kinematic_viscosity * (
                    dy * (u(i + 1, j) - u(i, j)) / dx
                    ) / (dx * dy);

                if (j != 0)
                {
                    // Bottom: u diff
                    u_dif(i, j) += fluid.kinematic_viscosity * (
                        dx * (u(i, j - 1) - u(i, j)) / dy
                        ) / (dx * dy);
                }

                if (j != ny - 1)
                {
                    // Top: u diff
                    u_dif(i, j) += fluid.kinematic_viscosity * (
                        dx * (u(i, j + 1) - u(i, j)) / dy
                        ) / (dx * dy);
                }
            }
        }

        for (int i = 0; i < nx; ++i)
        {
            for (int j = 1; j < ny; ++j)
            {
                // Inline Left: v diff
                v_dif(i, j) += fluid.kinematic_viscosity * (
                    dx * (v(i, j - 1) - v(i, j)) / dy
                    ) / (dx * dy);

                // Inline Right: v diff
                v_dif(i, j) += fluid.kinematic_viscosity * (
                    dx * (v(i, j + 1) - v(i, j)) / dy
                    ) / (dx * dy);

                if (i != 0)
                {
                    // Bottom: v diff
                    v_dif(i, j) += fluid.kinematic_viscosity * (
                        dy * (v(i - 1, j) - v(i, j)) / dx
                        ) / (dx * dy);
                }

                if (i != nx - 1)
                {
                    // Top: v diff
                    v_dif(i, j) += fluid.kinematic_viscosity * (
                        dy * (v(i + 1, j) - v(i, j)) / dx
                        ) / (dx * dy);
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
            p_old = p;
            for (int i = 1; i < nx - 1; ++i)
            {
                for (int j = 1; j < ny - 1; ++j)
                {
                    p(i, j) = 0.25 * (p_old(i - 1, j) + p_old(i, j - 1) + p_old(i + 1, j) + p_old(i, j + 1) - (divergence(i, j) * dx * dy) / (fluid.density * dt));
                }
            }
            ApplyPressureBoundaryConditions();
        }

        for (int i = 1; i < nx - 1; ++i)
        {
            for (int j = 1; j < ny - 1; ++j)
            {
                p(i, j) = p(i, j) * (1 - relaxation) + p_old(i, j) * relaxation;
            }
        }
    }

    void ApplyPressureCorrection()
    {
        for (int i = 1; i < nx; ++i) 
        {
            for (int j = 0; j < ny; ++j) 
            {
                u(i, j) -= (dt / dx) * (p(i, j) - p(i - 1, j)) / fluid.density;
            }
        }
        
        for (int i = 0; i < nx; ++i) 
        {
            for (int j = 1; j < ny; ++j) 
            {
                v(i, j) -= (dt / dx) * (p(i, j) - p(i, j - 1)) / fluid.density;
            }
        }
    }

    void SolveVelocitiesForMomentumEquation()
    {
        for (int i = 0; i < nx; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                u(i, j) = u_old(i, j) + dt * (u_flux(i, j) + u_dif(i, j) + u_scr(i,j));

                v(i, j) = v_old(i, j) + dt * (v_flux(i, j) + v_dif(i, j) + v_scr(i, j));
            }
        }
    }

    // Main time stepping method
    void solve(int num_iterations) 
    {
        //auto plotter = CFDVisualizer(nx, ny, (float)dx, (float)dy, cell_u, cell_v, p, cell_data_mutex);
        //auto plotterThread = std::jthread(&CFDVisualizer::Render, plotter);

        for (int iter = 0; iter < num_iterations; ++iter) 
        {
            // Store old fields and reset fluxes
            u_old = u;//swap
            v_old = v;//swap
            p_previous = p;//swap
            v_flux.reset();
            u_flux.reset();
            v_dif.reset();
            u_dif.reset();

            ComputeConvectiveFluxes();

            ComputeDiffusiveTerms();

            SolveVelocitiesForMomentumEquation();

            ComputeDivergence();

            SolvePressure();

            ApplyPressureCorrection();

            ApplyVelocityBoundaryConditions();

            // Check Residuals
            double maxResidual{ 0.0 };
            for (int i = 0; i < p.values.size(); ++i)
            {
                double residual = (p.values[i] - p_previous.values[i]) / (p.values[i] + 1e-20);
                maxResidual = std::max(maxResidual, residual);
            }
            printf("\rIteration: %d | Pressure Residual: %f", iter, maxResidual);
            if (maxResidual < residualLimit && iter > 2) break;
        }

        CalculateCellVelocities();

        auto plotter = CFDVisualizer(nx, ny, (float)dx, (float)dy, cell_u, cell_v, p, cell_data_mutex);
        plotter.Render();
    }

    // Simplified boundary conditions
    void ApplyVelocityBoundaryConditions() 
    {
        for (int i = 1; i < nx; ++i)
        {
            // Top boundary
            u_scr(i, 0) = 0.5 * (1.0 - u(i,0));
        }
    }

    void ApplyPressureBoundaryConditions()
    {
        for (int j = 1; j < ny - 1; ++j)
        {
            // Left boundary
            p(0, j) = p(1, j);

            // Right boundary
            p(nx - 1, j) = p(nx - 2, j);
        }

        for (int i = 1; i < nx - 1; ++i)
        {
            // Top boundary
            p(i, 0) = p(i, 1);

            // Bottom boundary
            p(i, ny - 1) = p(i, ny - 2);
        }

        p(0, 0) = (p(0, 1) + p(1, 0)) / 2.0;
        p(nx - 1, 0) = (p(nx - 1, 1) + p(nx - 2, 0)) / 2.0;
        p(0, ny - 1) = (p(0, ny - 2) + p(1, ny - 1)) / 2.0;
        p(nx - 1, ny - 1) = (p(nx - 1, ny - 2) + p(nx - 2, ny - 1)) / 2.0;
    }

    struct FluidProperties
    {
        double density{ 1.0 };
        double kinematic_viscosity{ 0.1 };
    };

    FluidProperties fluid;

    int nx, ny;
    double dx, dy;
    double dt{ 0.0001 };

    double relaxation{ 0.5 };

    int innerPressureItterations{ 10 };
    double residualLimit{ 0.0001 };

    PhysicalField p{ nx, ny };
    PhysicalField p_old{ nx, ny };
    PhysicalField p_previous{ nx, ny };

    PhysicalField cell_u{ nx, ny };
    PhysicalField cell_v{ nx, ny };

    PhysicalField u{ nx + 1, ny }, v{ nx, ny + 1 }; // staggered grid
    PhysicalField u_old{ nx + 1, ny }, v_old{ nx, ny + 1 }; // staggered grid

    PhysicalField divergence{ nx, ny };

    PhysicalField u_flux{ nx + 1, ny }; // staggered grid
    PhysicalField v_flux{ nx, ny + 1 }; // staggered grid

    PhysicalField u_dif{ nx + 1, ny }; // staggered grid
    PhysicalField v_dif{ nx, ny + 1 }; // staggered grid

    PhysicalField u_scr{ nx + 1, ny }; // staggered grid
    PhysicalField v_scr{ nx, ny + 1 }; // staggered grid

    FairMutex cell_data_mutex;
};