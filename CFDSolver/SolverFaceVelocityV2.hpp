#pragma once

#include <vector>
#include <iostream>
#include <cmath>
#include <cmath>
#include <limits>

#include "PhysicalField.hpp"
#include "Plotter.hpp"

class SolverFaceVelocity_V2
{
public:
    SolverFaceVelocity_V2(int nx, int ny, double dx, double dy) : nx(nx), ny(ny), dx(dx), dy(dy)
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
                divergence(i, j) = (u(i + 1,j) - u(i, j)) / dx + (v(i,j + 1) - v(i, j)) / dy;
            }
        }
    }

    void ComputeConvectiveFluxes()
    {
        for (int i = 1; i < nx; ++i)
        {
            for (int j = 1; j < ny - 1; ++j)
            {
                // U-velocity flux (central difference)
                u_flux(i, j) = -u(i, j) * (u(i, j) - u(i - 1, j)) / dx - v(i, j) * (u(i, j) - u(i,j - 1)) / dy;
            }
        }
        for (int i = 1; i < nx - 1; ++i)
        {
            for (int j = 1; j < ny; ++j)
            {
                // V-velocity flux (central difference)
                v_flux(i, j) = -u(i, j) * (v(i, j) - v(i - 1, j)) / dx - v(i, j) * (v(i, j) - v(i, j - 1)) / dy;
            }
        }
    }

    // Compute diffusive terms (on faces)
    void ComputeDiffusiveTerms() 
    {
        for (int i = 1; i < nx - 1; ++i)
        {
            for (int j = 1; j < ny - 1; ++j)
            {
                // U-velocity diffusion
                u_dif(i, j) = fluid.kinematic_viscosity * (
                    (u(i + 1, j) - 2 * u(i, j) + u(i - 1, j)) / (dx * dx) +
                    (u(i, j + 1) - 2 * u(i, j) + u(i, j - 1)) / (dy * dy));
    
                // V-velocity diffusion
                v_dif(i, j) = fluid.kinematic_viscosity * (
                    (v(i + 1, j) - 2 * v(i, j) + v(i - 1, j)) / (dx * dx) +
                    (v(i, j + 1) - 2 * v(i, j) + v(i, j - 1)) / (dy * dy));
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
                cell_v(i, j) = (v(i, j) + v(i + 1, j)) / 2.0;
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
    }

    void ApplyPressureCorrection()
    {
        for (int i = 1; i < nx; ++i) 
        {
            for (int j = 0; j < ny; ++j) 
            {
                u(i, j) -= dt / dx * (p(i, j) - p(i - 1, j)) / fluid.density;
            }
        }
        
        for (int i = 0; i < nx; ++i) 
        {
            for (int j = 1; j < ny; ++j) 
            {
                v(i, j) -= dt / dx * (p(i, j) - p(i, j - 1)) / fluid.density;
            }
        }
    }

    void SolveVelocitiesForMomentumEquation()
    {
        for (int i = 1; i < nx - 1; ++i)
        {
            for (int j = 1; j < ny - 1; ++j)
            {
                u(i, j) = u_old(i, j) + dt * (u_flux(i, j) + u_dif(i, j));

                v(i, j) = v_old(i, j) + dt * (v_flux(i, j) + v_dif(i, j));
            }
        }
    }

    // Main time stepping method
    void solve(int num_iterations) 
    {
        for (int iter = 0; iter < num_iterations; ++iter) 
        {
            // Store old fields
            u_old = u;
            v_old = v;
            p_previous = p;

            ApplyVelocityBoundaryConditions();

            ComputeDivergence();

            SolvePressure();

            ApplyPressureCorrection();

            ComputeConvectiveFluxes();

            ComputeDiffusiveTerms();

            SolveVelocitiesForMomentumEquation();

            // Check Residuals
            double maxResidual{ 0.0 };
            for (int i = 0; i < p.values.size(); ++i)
            {
                double residual = (p.values[i] - p_previous.values[i]) / (p.values[i] + 1e-20);
                maxResidual = std::max(maxResidual, residual);
            }
            printf("\rIteration: %d | Pressure Residual: %f", iter, maxResidual);
            if (maxResidual < residualLimit && iter > 1) break;
        }

        CalculateCellVelocities();

        // Plot Results
        auto plotter = CFDVisualizer(nx, ny, (float)dx, (float)dy, cell_u, cell_v, p);
        plotter.Render();
    }

    // Simplified boundary conditions
    void ApplyVelocityBoundaryConditions() 
    {
        for (int j = 0; j < ny + 1; ++j)
        {
            // Left boundary
            v(0, j) = 0.0;

            // Right boundary
            v(nx - 1, j) = 0.0;
        }
        for (int j = 0; j < ny; ++j)
        {
            // Left boundary
            u(0, j) = 0.0;

            // Right boundary
            u(nx - 1, j) = 0.0;
        }

        for (int i = 0; i < nx + 1; ++i)
        {
            // Top boundary
            u(i, 0) = 1.0;

            // Bottom boundary
            u(i, ny - 1) = 0.0;
        }

        for (int i = 0; i < nx; ++i)
        {
            // Top boundary
            v(i, 0) = 0.0;

            // Bottom boundary
            v(i, ny - 1) = 0.0;
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

        for (int i = 0; i < nx; ++i)
        {
            // Top boundary
            p(i, 0) = p(i, 1);

            // Bottom boundary
            p(i, ny - 1) = p(i, ny - 2);
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
    double dt{ 0.00001 };

    int innerPressureItterations{ 1000 };
    double residualLimit{ 0.00001 };

    PhysicalField p{ nx, ny };
    PhysicalField p_old{ nx, ny };
    PhysicalField p_previous{ nx, ny };

    PhysicalField cell_u{ nx, ny };
    PhysicalField cell_v{ nx, ny };

    PhysicalField u{ nx + 1, ny }, v{ nx, ny + 1 };
    PhysicalField u_old{ nx + 1, ny }, v_old{ nx, ny + 1 };

    PhysicalField divergence{ nx, ny };

    PhysicalField u_flux{ nx + 1, ny };
    PhysicalField v_flux{ nx, ny + 1 };

    PhysicalField u_dif{ nx + 1, ny };
    PhysicalField v_dif{ nx, ny + 1 };
};