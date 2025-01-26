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
            std::cout << "Solution Will Be Unstable\n";
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
                // Inline Left: u flux
                flux = dy * ((u(i, j) + u(i - 1, j)) / 2.0);
                value = flux >= 0.0 ? u(i - 1, j) : u(i, j);
                coeff = flux / (dx * dy);
                u_flux(i, j) += coeff * value;
                u_coeff(i, j) += coeff;

                // Inline Right: u flux
                flux = dy * -((u(i, j) + u(i + 1, j)) / 2.0);
                value = flux >= 0.0 ? u(i + 1, j) : u(i, j);
                coeff = flux / (dx * dy);
                u_flux(i, j) += coeff * value;
                u_coeff(i, j) += coeff;

                if (j != 0)
                {
                    // Bottom Left: u flux
                    flux = (dx / 2.0) * v(i - 1, j);
                    value = flux >= 0.0 ? u(i, j - 1) : u(i, j);
                    coeff = flux / (dx * dy);
                    u_flux(i, j) += coeff * value;
                    u_coeff(i, j) += coeff;

                    // Bottom Right: u flux
                    flux = (dx / 2.0) * v(i, j);
                    value = flux >= 0.0 ? u(i, j - 1) : u(i, j);
                    coeff = flux / (dx * dy);
                    u_flux(i, j) += coeff * value;
                    u_coeff(i, j) += coeff;
                }

                if (j != ny - 1)
                {
                    // Top Left: u flux
                    flux = -(dx / 2.0) * v(i - 1, j + 1);
                    value = flux >= 0.0 ? u(i, j + 1) : u(i, j);
                    coeff = flux / (dx * dy);
                    u_flux(i, j) += coeff * value;
                    u_coeff(i, j) += coeff;

                    // Top Right: u flux
                    flux = -(dx / 2.0) * v(i, j + 1);
                    value = flux >= 0.0 ? u(i, j + 1) : u(i, j);
                    coeff = flux / (dx * dy);
                    u_flux(i, j) += coeff * value;
                    u_coeff(i, j) += coeff;
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
                coeff = flux / (dx * dy);
                v_flux(i, j) += coeff * value;
                v_coeff(i, j) += coeff;

                // Inline Right: v flux
                flux = dy * -((v(i, j) + v(i, j + 1)) / 2.0);
                value = flux >= 0.0 ? v(i, j + 1) : v(i, j);
                coeff = flux / (dx * dy);
                v_flux(i, j) += coeff * value;
                v_coeff(i, j) += coeff;

                if (i != 0)
                {
                    // Bottom Left: v flux
                    flux = (dy / 2.0) * u(i, j - 1);
                    value = flux >= 0.0 ? v(i - 1, j) : v(i, j);
                    coeff = flux / (dx * dy);
                    v_flux(i, j) += coeff * value;
                    v_coeff(i, j) += coeff;

                    // Bottom Right: v flux
                    flux = (dy / 2.0) * u(i, j);
                    value = flux >= 0.0 ? v(i - 1, j) : v(i, j);
                    coeff = flux / (dx * dy);
                    v_flux(i, j) += coeff * value;
                    v_coeff(i, j) += coeff;
                }

                if (i != nx - 1)
                {
                    // Top Left: v flux
                    flux = -(dy / 2.0) * u(i + 1, j - 1);
                    value = flux >= 0.0 ? v(i + 1, j) : v(i, j);
                    coeff = flux / (dx * dy);
                    v_flux(i, j) += coeff * value;
                    v_coeff(i, j) += coeff;

                    // Top Right: v flux
                    flux = -(dy / 2.0) * u(i + 1, j);
                    value = flux >= 0.0 ? v(i + 1, j) : v(i, j);
                    coeff = flux / (dx * dy);
                    v_flux(i, j) += coeff * value;
                    v_coeff(i, j) += coeff;
                }
            }
        }
    }

    void ComputeDiffusiveTerms() 
    {
        double coeff = 0.0;
        for (int i = 1; i < nx; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                // Inline Left: u diff
                coeff = fluid.kinematic_viscosity * (dy / dx) / (dx * dy);
                u_dif(i, j) += coeff * (u(i - 1, j) - u(i, j));
                u_coeff(i, j) += coeff;

                // Inline Right: u diff
                coeff = fluid.kinematic_viscosity * (dy / dx) / (dx * dy);
                u_dif(i, j) += coeff * (u(i + 1, j) - u(i, j));
                u_coeff(i, j) += coeff;

                if (j != 0)
                {
                    // Bottom: u diff
                    coeff = fluid.kinematic_viscosity * (dx / dy) / (dx * dy);
                    u_dif(i, j) += coeff * (u(i, j - 1) - u(i, j));
                    u_coeff(i, j) += coeff;
                }

                if (j != ny - 1)
                {
                    // Top: u diff
                    coeff = fluid.kinematic_viscosity * (dx / dy) / (dx * dy);
                    u_dif(i, j) += coeff * (u(i, j + 1) - u(i, j));
                    u_coeff(i, j) += coeff;
                }
            }
        }

        for (int i = 0; i < nx; ++i)
        {
            for (int j = 1; j < ny; ++j)
            {
                // Inline Left: v diff
                coeff = fluid.kinematic_viscosity * (dx / dy) / (dx * dy);
                v_dif(i, j) += coeff * (v(i, j - 1) - v(i, j));
                v_coeff(i, j) += coeff;

                // Inline Right: v diff
                coeff = fluid.kinematic_viscosity * (dx / dy) / (dx * dy);
                v_dif(i, j) += coeff * (v(i, j + 1) - v(i, j));
                v_coeff(i, j) += coeff;

                if (i != 0)
                {
                    // Bottom: v diff
                    coeff = fluid.kinematic_viscosity * (dy / dx) / (dx * dy);
                    v_dif(i, j) += coeff * (v(i - 1, j) - v(i, j));
                    v_coeff(i, j) += coeff;
                }

                if (i != nx - 1)
                {
                    // Top: v diff
                    coeff = fluid.kinematic_viscosity * (dy / dx) / (dx * dy);
                    v_dif(i, j) += coeff * (v(i + 1, j) - v(i, j));
                    v_coeff(i, j) += coeff;
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
            p_correction_old = p_correction;
            p_correction.reset();
            for (int i = 0; i < nx; ++i)
            {
                for (int j = 0; j < ny; ++j)
                {
                    double Ap = u_coeff(i, j) + u_coeff(i + 1, j) + v_coeff(i, j) + v_coeff(i, j + 1);

                    if (i != 0)
                        p_correction(i, j) += (u_coeff(i, j) / Ap) * p_correction_old(i - 1, j);
                    if (i != nx - 1)
                        p_correction(i, j) += (u_coeff(i + 1, j) / Ap) * p_correction_old(i + 1, j);
                    if (j != 0)
                        p_correction(i, j) += (v_coeff(i, j) / Ap) * p_correction_old(i, j - 1);
                    if (j != ny - 1)
                        p_correction(i, j) += (v_coeff(i, j + 1) / Ap) * p_correction_old(i, j + 1);

                    p_correction(i, j) -= divergence(i, j);
                }
            }
        }

        for (int i = 0; i < nx; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                p(i, j) += p_correction(i, j);
            }
        }

        //ApplyPressureBoundaryConditions();
    }

    void ApplyPressureCorrection()
    {
        for (int i = 1; i < nx; ++i) 
        {
            for (int j = 0; j < ny; ++j) 
            {
                u(i, j) += dt * (p(i - 1, j) - p(i, j)) / (dx * fluid.density); //p_correction
            }
        }
        
        for (int i = 0; i < nx; ++i) 
        {
            for (int j = 1; j < ny; ++j) 
            {
                v(i, j) += dt * (p(i, j - 1) - p(i, j)) / (dy * fluid.density); //p_correction
            }
        }
    }

    void SolveVelocitiesForMomentumEquation()
    {
        for (int i = 1; i < nx; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                u(i, j) = u_old(i, j) + dt * (u_flux(i, j) + u_dif(i, j) + u_scr(i, j)/* + (p(i - 1, j) - p(i, j)) / dy*/);
            }
        }
        for (int i = 0; i < nx; ++i)
        {
            for (int j = 1; j < ny; ++j)
            {

                v(i, j) = v_old(i, j) + dt * (v_flux(i, j) + v_dif(i, j) + v_scr(i, j)/* + (p(i, j - 1) - p(i, j)) / dx*/);
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
            p_correction.reset();
            u_coeff.reset();
            v_coeff.reset();
            
            ApplyVelocityBoundaryConditions();

            ComputeConvectiveFluxes();

            ComputeDiffusiveTerms();

            SolveVelocitiesForMomentumEquation();

            ComputeDivergence();

            SolvePressure();

            ApplyPressureCorrection();

            // Check Residuals
            if (iter % 5 == 0)
            {
                double maxResidual{ 0.0 };
                for (int i = 0; i < p.values.size(); ++i)
                {
                    double residual = (p.values[i] - p_previous.values[i]) / (p.values[i] + 1e-20);
                    maxResidual = std::max(maxResidual, residual);
                }
                printf("\rIteration: %d | Pressure Residual: %f  ", iter, maxResidual);
                if (maxResidual < residualLimit && iter > 2) break;
            }
        }

        CalculateCellVelocities();

        auto plotter = CFDVisualizer(nx, ny, (float)dx, (float)dy, cell_u, cell_v, u, v, p, cell_data_mutex);
        plotter.Render();
    }

    // Simplified boundary conditions
    void ApplyVelocityBoundaryConditions() 
    {
        for (int i = 1; i < nx; ++i)
        {
            // Top boundary
            u_scr(i, 0) = (1.0 - u(i, 0));

            // Bottom boundary
            u_scr(i, ny - 1) = (0.0 - u(i, ny - 1));
        }
        for (int j = 1; j < ny; ++j)
        {
            // Left boundary
            v_scr(0, j) = (0.0 - v(0, j));

            // Right boundary
            v_scr(nx - 1, j) = (0.0 - v(nx - 1, j));
        }

        //u_scr(0, ny / 2) = u_scr(nx, ny / 2) = 1.0;
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
    double dt{ 0.00001 };

    //double relaxation{ 0.0 };

    int innerPressureItterations{ 20 };
    double residualLimit{ 0.000001 };

    PhysicalField p{ nx, ny };
    PhysicalField p_previous{ nx, ny };
    PhysicalField p_correction{ nx, ny };
    PhysicalField p_correction_old{ nx, ny };

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

    PhysicalField u_coeff{ nx + 1, ny }; // staggered grid
    PhysicalField v_coeff{ nx, ny + 1 }; // staggered grid

    FairMutex cell_data_mutex;
};