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
                //if (i != 0)
                {
                    // Inline Left: u flux
                    flux = dy * ((u(i, j) + u(i - 1, j)) / 2.0);
                    //if (massFlux > 0.0)
                    {
                        value = flux >= 0.0 ? u(i - 1, j) : u(i, j);
                        coeff = flux / (dx * dy);
                        u_coeff(i, j) += coeff;
                        u_flux(i, j) += coeff * value;
                    }
                }

                //if (i != nx - 1)
                {
                    // Inline Right: u flux
                    flux = dy * -((u(i, j) + u(i + 1, j)) / 2.0);
                    //if (massFlux > 0.0)
                    {
                        value = flux >= 0.0 ? u(i + 1, j) : u(i, j);
                        coeff = flux / (dx * dy);
                        u_coeff(i, j) += coeff;
                        u_flux(i, j) += coeff * value;
                    }
                }

                if (j != 0)
                {
                    // Bottom Left: u flux
                    flux = (dx / 2.0) * v(i - 1, j);
                    //if (massFlux > 0.0)
                    {
                        value = flux >= 0.0 ? u(i, j - 1) : u(i, j);
                        coeff = flux / (dx * dy);
                        u_coeff(i, j) += coeff;
                        u_flux(i, j) += coeff * value;
                    }

                    // Bottom Right: u flux
                    flux = (dx / 2.0) * v(i, j);
                    //if (massFlux > 0.0)
                    {
                        value = flux >= 0.0 ? u(i, j - 1) : u(i, j);
                        coeff = flux / (dx * dy);
                        u_coeff(i, j) += coeff;
                        u_flux(i, j) += coeff * value;
                    }
                }

                if (j != ny - 1)
                {
                    // Top Left: u flux
                    flux = -(dx / 2.0) * v(i - 1, j + 1);
                    //if (massFlux > 0.0)
                    {
                        value = flux >= 0.0 ? u(i, j + 1) : u(i, j);
                        coeff = flux / (dx * dy);
                        u_coeff(i, j) += coeff;
                        u_flux(i, j) += coeff * value;
                    }

                    // Top Right: u flux
                    flux = -(dx / 2.0) * v(i, j + 1);
                    //if (massFlux > 0.0)
                    {
                        value = flux >= 0.0 ? u(i, j + 1) : u(i, j);
                        coeff = flux / (dx * dy);
                        u_coeff(i, j) += coeff;
                        u_flux(i, j) += coeff * value;
                    }
                }
            }
        }

        for (int i = 0; i < nx; ++i)
        {
            for (int j = 1; j < ny; ++j)
            {
                //if (j != 0)
                {
                    // Inline Left: v flux
                    flux = dy * ((v(i, j) + v(i, j - 1)) / 2.0);
                    //if (massFlux > 0.0)
                    {
                        value = flux >= 0.0 ? v(i, j - 1) : v(i, j);
                        coeff = flux / (dx * dy);
                        v_coeff(i, j) += coeff;
                        v_flux(i, j) += coeff * value;
                    }
                }

                //if (j != ny - 1)
                {
                    // Inline Right: v flux
                    flux = dy * -((v(i, j) + v(i, j + 1)) / 2.0);
                    //if (massFlux > 0.0)
                    {
                        value = flux >= 0.0 ? v(i, j + 1) : v(i, j);
                        coeff = flux / (dx * dy);
                        v_coeff(i, j) += coeff;
                        v_flux(i, j) += coeff * value;
                    }
                }

                if (i != 0)
                {
                    // Bottom Left: v flux
                    flux = (dy / 2.0) * u(i, j - 1);
                    //if (massFlux > 0.0)
                    {
                        value = flux >= 0.0 ? v(i - 1, j) : v(i, j);
                        coeff = flux / (dx * dy);
                        v_coeff(i, j) += coeff;
                        v_flux(i, j) += coeff * value;
                    }

                    // Bottom Right: v flux
                    flux = (dy / 2.0) * u(i, j);
                    //if (massFlux > 0.0)
                    {
                        value = flux >= 0.0 ? v(i - 1, j) : v(i, j);
                        coeff = flux / (dx * dy);
                        v_coeff(i, j) += coeff;
                        v_flux(i, j) += coeff * value;
                    }
                }

                if (i != nx - 1)
                {
                    // Top Left: v flux
                    flux = -(dy / 2.0) * u(i + 1, j - 1);
                    //if (massFlux > 0.0)
                    {
                        value = flux >= 0.0 ? v(i + 1, j) : v(i, j);
                        coeff = flux / (dx * dy);
                        v_coeff(i, j) += coeff;
                        v_flux(i, j) += coeff * value;
                    }

                    // Top Right: v flux
                    flux = -(dy / 2.0) * u(i + 1, j);
                    //if (massFlux > 0.0)
                    {
                        value = flux >= 0.0 ? v(i + 1, j) : v(i, j);
                        coeff = flux / (dx * dy);
                        v_coeff(i, j) += coeff;
                        v_flux(i, j) += coeff * value;
                    }
                }
            }
        }
    }

    void ComputeDiffusiveTerms()  // are we adding boundary faces - maybe best to add the blocked face check now
    {
        double coeff = 0.0;
        double value = 0.0;
        for (int i = 1; i < nx; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                // Inline Left: u diff
                coeff = fluid.kinematic_viscosity / (dx * dx);
                value = (u(i - 1, j) - u(i, j));
                u_coeff(i, j) += coeff;
                u_dif(i, j) += coeff * value;

                // Inline Right: u diff
                coeff = fluid.kinematic_viscosity / (dx * dx);
                value = (u(i + 1, j) - u(i, j));
                u_coeff(i, j) += coeff;
                u_dif(i, j) += coeff * value;

                if (j != 0)
                {
                    // Bottom: u diff
                    coeff = fluid.kinematic_viscosity / (dy * dy);
                    value = (u(i, j - 1) - u(i, j));
                    u_coeff(i, j) += coeff;
                    u_dif(i, j) += coeff * value;
                }

                if (j != ny - 1)
                {
                    // Top: u diff
                    coeff = fluid.kinematic_viscosity / (dy * dy);
                    value = (u(i, j + 1) - u(i, j));
                    u_coeff(i, j) += coeff;
                    u_dif(i, j) += coeff * value;
                }
            }
        }

        for (int i = 0; i < nx; ++i)
        {
            for (int j = 1; j < ny; ++j)
            {
                // Inline Left: v diff
                coeff = fluid.kinematic_viscosity / (dy * dy);
                value = (v(i, j - 1) - v(i, j));
                v_coeff(i, j) += coeff;
                v_dif(i, j) += coeff * value;

                // Inline Right: v diff
                coeff = fluid.kinematic_viscosity / (dy * dy);
                value = (v(i, j + 1) - v(i, j));
                v_coeff(i, j) += coeff;
                v_dif(i, j) += coeff * value;

                if (i != 0)
                {
                    // Bottom: v diff
                    coeff = fluid.kinematic_viscosity / (dx * dx);
                    value = (v(i - 1, j) - v(i, j));
                    v_coeff(i, j) += coeff;
                    v_dif(i, j) += coeff * value;
                }

                if (i != nx - 1)
                {
                    // Top: v diff
                    coeff = fluid.kinematic_viscosity/ (dx * dx);
                    value = (v(i + 1, j) - v(i, j));
                    v_coeff(i, j) += coeff;
                    v_dif(i, j) += coeff * value;
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
                u(i, j) += dt * (u_flux(i, j) + u_dif(i, j) + u_scr(i, j)/* + (p(i - 1, j) - p(i, j)) / dy*/);
            }
        }
        for (int i = 0; i < nx; ++i)
        {
            for (int j = 1; j < ny; ++j)
            {

                v(i, j) += dt * (v_flux(i, j) + v_dif(i, j) + v_scr(i, j)/* + (p(i, j - 1) - p(i, j)) / dx*/);
            }
        }
    }

    // Main time stepping method
    void solve(int num_iterations) 
    {
        for (int iter = 0; iter < num_iterations; ++iter) 
        {
            p_previous = p;
            v_flux.reset();
            u_flux.reset();
            v_dif.reset();
            u_dif.reset();
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
                    double residual = std::abs((p.values[i] - p_previous.values[i]) / (p.values[i] + 1e-20));
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

    // Simplified boundary conditions
    void ApplyVelocityBoundaryConditions() 
    {
        double coeff = 0.0;

        // Inflow
        //double inflowVelocity = 1.0;
        //double massFlux = inflowVelocity * fluid.density * dy;
        //coeff = massFlux / (dx * dy);
        //for (int j = 0; j < ny; ++j)
        //{
        //
        //    //u_coeff(1, j) += coeff;
        //    //u_scr(1, j) += coeff * inflowVelocity;
        //    
        //    //u_coeff(nx - 1, j) += coeff;
        //    //u_scr(nx - 1, j) += coeff * inflowVelocity;
        //
        //    //u(1, j) = u(nx - 1, j) = 1.0;
        //}

        // Friction
        coeff = fluid.kinematic_viscosity / (dy / 2);
        for (int i = 1; i < nx; ++i)
        {
            // Top boundary
            u_coeff(i, 0) += coeff;
            u_scr(i, 0) += (1.0 - u(i, 0)) * coeff;
        
            // Bottom boundary
            u_coeff(i, ny - 1) += coeff;
            u_scr(i, ny - 1) += (0.0 - u(i, ny - 1)) * coeff;
        }
        coeff = fluid.kinematic_viscosity / (dx / 2);
        for (int j = 1; j < ny; ++j)
        {
            // Left boundary
            v_coeff(0, j) += coeff;
            v_scr(0, j) += (0.0 - v(0, j)) * coeff;
        
            // Right boundary
            v_coeff(nx - 1, j) += coeff;
            v_scr(nx - 1, j) += (0.0 - v(nx - 1, j)) * coeff;
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
    double dt{ 10e-6 };
    double relaxation{ 1.0 };

    int innerPressureItterations{ 20 };
    double residualLimit{ 10e-7 };
    double divergenceLimit{ 10e-10 };

    PhysicalField p{ nx, ny };
    PhysicalField p_scr{ nx, ny };
    PhysicalField p_previous{ nx, ny };
    PhysicalField p_correction{ nx, ny };
    PhysicalField p_correction_old{ nx, ny };

    PhysicalField cell_u{ nx, ny };
    PhysicalField cell_v{ nx, ny };

    PhysicalField u{ nx + 1, ny }, v{ nx, ny + 1 }; // staggered grid

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

    // TODO
    // - blocked cells / faces
    // - inlet velocity boundary (what are the pressure conditions)
    // - friction
    // - turbulance
    // - plot divergence as it should be zero
};