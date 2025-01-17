#pragma once

#include <vector>
#include <iostream>
#include <cmath>
#include <cmath>
#include <limits>

#include "PhysicalField.hpp"
#include "Plotter.hpp"

class Solver 
{
public:
    Solver(int nx, int ny, double dx, double dy) : nx(nx), ny(ny), dx(dx), dy(dy)
    {
        if (dt > 0.5 * std::pow(std::min(dx, dy), 2) / fluid.kinematic_viscosity)
            std::cout << "Unstable\n";
    }

    // Calculate convective terms
    void ComputeConvectiveFluxes() 
    {
        for (int i = 1; i < nx - 1; ++i)
        {
            for (int j = 1; j < ny - 1; ++j) 
            {
                // U-velocity flux (central difference)
                u_flux(i, j) = (u(i + 1, j) - u(i - 1, j)) * u(i, j) / (2 * dx) +
                               (u(i, j + 1) - u(i, j - 1)) * v(i, j) / (2 * dy);

                // V-velocity flux (central difference)
                v_flux(i, j) = (v(i + 1, j) - v(i - 1, j)) * u(i, j) / (2 * dx) +
                               (v(i, j + 1) - v(i, j - 1)) * v(i, j) / (2 * dy);
            }
        }
    }

    // Compute diffusive terms
    void ComputeDiffusiveTerms() 
    {
        for (int i = 1; i < nx - 1; ++i)
        {
            for (int j = 1; j < ny - 1; ++j)
            {
                // U-velocity diffusion
                u_dif(i, j) = fluid.kinematic_viscosity * (
                    (u(i + 1, j) - 2 * u(i, j) + u(i - 1, j)) / (dx * dx) +
                    (u(i, j + 1) - 2 * u(i, j) + u(i, j - 1)) / (dy * dy)
                    );

                // V-velocity diffusion
                v_dif(i, j) = fluid.kinematic_viscosity * (
                    (v(i + 1, j) - 2 * v(i, j) + v(i - 1, j)) / (dx * dx) +
                    (v(i, j + 1) - 2 * v(i, j) + v(i, j - 1)) / (dy * dy)
                    );
            }
        }
    }

    // Pressure correction step (SIMPLE algorithm)
    void PressureCorrection()
    {
        for (int i = 1; i < nx - 1; ++i)
        {
            for (int j = 1; j < ny - 1; ++j)
            {
                p_equation_rhs(i, j) = (fluid.density / dt) * 
                    ((u(i + 1, j) - u(i - 1, j)) / (2.0f * dx) +
                     (v(i, j + 1) - v(i, j - 1)) / (2.0f * dy));

                //p_equation_rhs(i, j) = (fluid.density * dx / 16) * (
                //    (2 / dt) * (u(i, j) - u(j, i - 1) + v(j + 1, i) - v(j - 1, i))
                //    - (2 / dx) * (u(j + 1, i) - u(j - 1, i)) * (v(j, i + 1) - v(j, i - 1))
                //    - (u(j, i + 1) - u(j, i - 1)) * (u(j, i + 1) - u(j, i - 1)) / dx
                //    - (v(j + 1, i) - v(j - 1, i)) * (v(j + 1, i) - v(j - 1, i)) / dx
                //    );
            }
        }

        for (int itter = 0; itter < innerPressureItterations; ++itter)
        {
            p_old = p;
            for (int i = 1; i < nx - 1; ++i)
            {
                for (int j = 1; j < ny - 1; ++j)
                {
                    //p(i, j) = (p_old(i + 1, j) + p_old(i - 1, j) + (p_old(i, j + 1) + p_old(i, j - 1))) / 4
                    //    - (fluid.density * dx / 16) * (
                    //        (2 / dt) * (u(i, j) - u(j, i - 1) + v(j + 1, i) - v(j - 1, i))
                    //        - (2 / dx) * (u(j + 1, i) - u(j - 1, i)) * (v(j, i + 1) - v(j, i - 1))
                    //        - (u(j, i + 1) - u(j, i - 1)) * (u(j, i + 1) - u(j, i - 1)) / dx
                    //        - (v(j + 1, i) - v(j - 1, i)) * (v(j + 1, i) - v(j - 1, i)) / dx
                    //        );

                    p(i, j) = 0.25 * (p_old(i - 1, j) + p_old(i, j - 1) + p_old(i + 1, j) + p_old(i, j + 1) - dx * dy * p_equation_rhs(i, j));
                }
            }
            ApplyPressureBoundaryConditions();
        }

        // Update pressure and correct velocities 
        for (int j = 1; j < ny - 1; ++j)
        {
            for (int i = 1; i < nx - 1; ++i)
            {
                // Correct velocities using the pressure gradient
                u(i, j) -= dt * (p(i, j) - p(i - 1, j)) / (2.0f * fluid.density * dx);
                v(i, j) -= dt * (p(i, j + 1) - p(i, j - 1)) / (2.0f * fluid.density * dy);
            }
        }
    }

    // Main time stepping method
    void solve(int num_iterations) 
    {
        for (int iter = 0; iter < num_iterations; ++iter) 
        {
            // Store old velocities
            u_old = u;
            v_old = v;

            ApplyVelocityBoundaryConditions();

            ComputeConvectiveFluxes();

            ComputeDiffusiveTerms();

            // Update velocities (simplified momentum equation)
            for (int i = 1; i < nx - 1; ++i)
            {
                for (int j = 1; j < ny - 1; ++j) 
                {                
                    u(i, j) = u_old(i, j) - dt * u_flux(i, j) + dt * u_dif(i, j);

                    v(i, j) = v_old(i, j) - dt * v_flux(i, j) + dt * v_dif(i, j);
                }
            }

            PressureCorrection();

            // Check Residuals
            double maxResidual{ 0.0 };
            for (int i = 0; i < u.values.size(); ++i)
            {
                double residual = (u.values[i] - u_old.values[i]) / (u.values[i] + 1e-20);
                maxResidual = std::max(maxResidual, residual);
            }
            printf("\rIteration: %d | X Velocity Residual: %f", iter, maxResidual);
            if (maxResidual < residualLimit && iter != 0) break;
        }

        //generateVectorFieldWithGrid(u, v);
        //auto plotter = CFDVisualizer(nx, ny, (float)dx, (float)dy, u, v, p);
        //plotter.Render();
    }

    // Simplified boundary conditions
    void ApplyVelocityBoundaryConditions() 
    {
        for (int j = 0; j < ny; ++j)
        {
            // Left boundary
            u(0, j) = -u(1, j);
            v(0, j) = 0.0;

            // Right boundary
            u(nx - 1, j) = -u(nx - 2, j);
            v(nx - 1, j) = 0.0;
        }

        for (int i = 0; i < nx; ++i)
        {
            // Top boundary
            u(i, 0) = 1.0;
            v(i, 0) = -v(i, 1);

            // Bottom boundary
            u(i, ny - 1) = 0.0;
            v(i, ny - 1) = -v(i, ny - 2);
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
            p(i, ny - 1) = 0.0;// p(i, ny - 2);
        }
    }

    struct FluidProperties
    {
        double density{ 1.0 };
        double kinematic_viscosity{ 0.1 /*0.000015*/ };
    };

    FluidProperties fluid;

    int nx, ny;
    double dx, dy;
    double dt{ 0.0000001 };

    int innerPressureItterations{ 20 };
    double residualLimit{ 0.0001 };

    PhysicalField p{ nx, ny };
    PhysicalField p_old{ nx, ny };
    PhysicalField p_equation_rhs{ nx, ny };

    PhysicalField u{ nx, ny }, v{ nx, ny };
    PhysicalField u_old{ nx, ny }, v_old{ nx, ny };

    PhysicalField u_flux{ nx, ny };
    PhysicalField v_flux{ nx, ny };

    PhysicalField u_dif{ nx, ny };
    PhysicalField v_dif{ nx, ny };

    PhysicalField mass_flux{ nx, ny };
};