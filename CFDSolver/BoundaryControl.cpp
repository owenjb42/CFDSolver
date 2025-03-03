#include "BoundaryControl.hpp"
#include "Solver_Staggered_IMEX_Temp.hpp"

/////////////////
// Fixed Inlet //
/////////////////

void FixedInletBoundaryCondition::ApplyForVelocity(SolverStaggeredIMEXTemp& solver)
{
    int offset = direction == 0 ? 1 : -1;
    int cell_offset = direction == 0 ? 0 : -1;
    if (component == 0) // u
    {
        double flux = inlet_velocity * solver.dy;
        double coeff = solver.fluid.density * flux / (solver.dx * solver.dy);
        for (auto [i, j] : boundary_faces)
        {
            solver.p_scr(i + cell_offset, j) -= inlet_velocity * solver.dy;
            solver.u_coeff(i + offset, j) += coeff;
            solver.u_scr(i + offset, j) += coeff * offset * inlet_velocity;
        }
    }
    else // v
    {
        double flux = inlet_velocity * solver.dx;
        double coeff = solver.fluid.density * flux / (solver.dx * solver.dy);
        for (auto [i, j] : boundary_faces)
        {
            solver.p_scr(i, j + cell_offset) -= inlet_velocity * solver.dx;
            solver.v_coeff(i, j + offset) += coeff;
            solver.v_scr(i, j + offset) += coeff * offset * inlet_velocity;
        }
    }
}

void FixedInletBoundaryCondition::ApplyForTemperature(SolverStaggeredIMEXTemp& solver)
{
    int offset = direction == 0 ? 0 : -1;
    if (component == 0) // u
    {
        double massflux = inlet_velocity * solver.fluid.density * solver.dy;
        double coeff = massflux / (solver.dx * solver.dy);

        for (auto [i, j] : boundary_faces)
        {
            solver.t_coeff(i + offset, j) += coeff;
            solver.t_scr(i + offset, j) += coeff * inlet_temperature;
        }
    }
    else // v
    {
        double massflux = inlet_velocity * solver.fluid.density * solver.dx;
        double coeff = massflux / (solver.dx * solver.dy);

        for (auto [i, j] : boundary_faces)
        {
            solver.t_coeff(i, j + offset) += coeff;
            solver.t_scr(i, j + offset) += coeff * inlet_temperature;
        }
    }
}

void FixedInletBoundaryCondition::CorrectBoundaryCellVelocities(SolverStaggeredIMEXTemp& solver)
{
    int offset = direction == 0 ? 1 : -1;
    int cell_offset = direction == 0 ? 0 : -1;
    if (component == 0) // u
    {
        for (auto [i, j] : boundary_faces)
        {
            solver.cell_u(i + cell_offset, j) += offset * inlet_velocity / 2.0;
        }
    }
    else // v
    {
        int cell_offset = direction == 0 ? 0 : -1;
        for (auto [i, j] : boundary_faces)
        {
            solver.cell_v(i, j + cell_offset) += offset * inlet_velocity / 2.0;
        }
    }
}

//////////////////
// Fixed Outlet //
//////////////////

void FixedOutletBoundaryCondition::ApplyForVelocity(SolverStaggeredIMEXTemp& solver)
{
    if (component == 0) // u
    {
        for (auto [i, j] : boundary_faces)
        {
            int offset = direction == 0 ? -1 : 1;
            int cell_offset = direction == 0 ? 0 : -1;
            solver.p_scr(i + cell_offset, j) += outlet_velocity * solver.dy;
        }
    }
    else // v
    {
        for (auto [i, j] : boundary_faces)
        {
            int offset = direction == 0 ? -1 : 1;
            int cell_offset = direction == 0 ? 0 : -1;
            solver.p_scr(i, j + cell_offset) += outlet_velocity * solver.dx;
        }
    }
}
void FixedOutletBoundaryCondition::CorrectBoundaryCellVelocities(SolverStaggeredIMEXTemp& solver)
{
    int offset = direction == 0 ? -1 : 1;
    int cell_offset = direction == 0 ? 0 : -1;
    if (component == 0) // u
    {
        for (auto [i, j] : boundary_faces)
        {
            solver.cell_u(i + cell_offset, j) += offset * outlet_velocity / 2.0;
        }
    }
    else // v
    {
        int cell_offset = direction == 0 ? 0 : -1;
        for (auto [i, j] : boundary_faces)
        {
            solver.cell_v(i, j + cell_offset) += offset * outlet_velocity / 2.0;
        }
    }
}

///////////////////
// Open Boundary //
///////////////////
void OpenBoundaryCondition::ApplyForPressure(SolverStaggeredIMEXTemp& solver)
{
    double fix_coef = 10.0;

    int offset = direction == 0 ? 1 : -1;
    int cell_offset = direction == 0 ? 0 : -1;
    if (component == 0) // u
    {
        for (auto [i, j] : boundary_faces)
        {
            double coeff = solver.dy * fix_coef;
            double dp = (pressure - solver.p(i + cell_offset, j));
            double u_b = -dp / solver.fluid.density;
            solver.p_scr(i + cell_offset, j) += coeff * u_b;
            solver.p_coef(i + cell_offset, j) += coeff;
        }
    }
    else // v
    {
        for (auto [i, j] : boundary_faces)
        {
            double coeff = solver.dx * fix_coef;
            double dp = (pressure - solver.p(i, j + cell_offset));
            double u_b = -dp / solver.fluid.density;
            solver.p_scr(i, j + cell_offset) += coeff * u_b;
            solver.p_coef(i, j + cell_offset) += coeff;
        }
    }
}
void OpenBoundaryCondition::ApplyForTemperature(SolverStaggeredIMEXTemp& solver)
{
    double fix_coef = 10.0;

    int offset = direction == 0 ? 0 : -1;
    int cell_offset = direction == 0 ? 0 : -1;
    if (component == 0) // u
    {
        for (auto [i, j] : boundary_faces)
        {
            double massflux = fix_coef * (pressure - solver.p(i + cell_offset, j)) * solver.dy;
            if (massflux > 0.0)
            {
                double coeff = massflux / (solver.dx * solver.dy);
                solver.t_coeff(i + offset, j) += coeff;
                solver.t_scr(i + offset, j) += coeff * temperature;
            }
        }
    }
    else // v
    {
        for (auto [i, j] : boundary_faces)
        {
            double massflux = fix_coef * (pressure - solver.p(i, j + cell_offset)) * solver.dx;
            if (massflux > 0.0)
            {
                double coeff = massflux / (solver.dx * solver.dy);
                solver.t_coeff(i, j + offset) += coeff;
                solver.t_scr(i, j + offset) += coeff * temperature;
            }
        }
    }
}
void OpenBoundaryCondition::CorrectBoundaryCellVelocities(SolverStaggeredIMEXTemp& solver)
{
    double fix_coef = 10.0;

    int offset = direction == 0 ? 1 : -1;
    int cell_offset = direction == 0 ? 0 : -1;
    if (component == 0) // u
    {
        for (auto [i, j] : boundary_faces)
        {
            double dp = (pressure - solver.p(i + cell_offset, j));
            double u_b = fix_coef * dp / solver.fluid.density;
            solver.cell_u(i + cell_offset, j) += offset * u_b / 2.0;
        }
    }
    else // v
    {
        for (auto [i, j] : boundary_faces)
        {
            double dp = (pressure - solver.p(i, j + cell_offset));
            double u_b = fix_coef * dp / solver.fluid.density;
            solver.cell_v(i, j + cell_offset) += offset * u_b / 2.0;
        }
    }
}

///////////////////////
// Friction Boundary //
///////////////////////
void FrictionBoundaryCondition::ApplyFriction(SolverStaggeredIMEXTemp& solver)
{
    double coeff = 0.5 * solver.fluid.kinematic_viscosity / (solver.dy / 2);
    for (int i = 0; i < solver.nx + 1; ++i)
    {
        for (int j = 0; j < solver.ny; ++j)
        {
            if (!solver.u_face_flags(i, j, Flag::Open))
            {
                if (i != 0)
                {
                    if (solver.v_face_flags(i - 1, j + 1, Flag::Open)) // Top Left
                    {
                        solver.v_coeff(i - 1, j + 1) += coeff;
                        solver.v_scr(i - 1, j + 1) += wall_velocity * coeff;
                    }
                    if (solver.v_face_flags(i - 1, j, Flag::Open)) // Bottom Left
                    {
                        solver.v_coeff(i - 1, j) += coeff;
                        solver.v_scr(i - 1, j) += wall_velocity * coeff;
                    }
                }
                if (i != solver.nx)
                {
                    if (solver.v_face_flags(i, j + 1, Flag::Open)) // Top Right
                    {
                        solver.v_coeff(i, j + 1) += coeff;
                        solver.v_scr(i, j + 1) += wall_velocity * coeff;
                    }
                    if (solver.v_face_flags(i, j, Flag::Open)) // Bottom Right
                    {
                        solver.v_coeff(i, j) += coeff;
                        solver.v_scr(i, j) += wall_velocity * coeff;
                    }
                }
            }
        }
    }
    coeff = 0.5 * solver.fluid.kinematic_viscosity / (solver.dx / 2);
    for (int i = 0; i < solver.nx; ++i)
    {
        for (int j = 0; j < solver.ny + 1; ++j)
        {
            if (!solver.v_face_flags(i, j, Flag::Open))
            {
                if (j != 0)
                {
                    if (solver.u_face_flags(i, j - 1, Flag::Open)) // Bottom Left
                    {
                        solver.u_coeff(i, j - 1) += coeff;
                        solver.u_scr(i, j - 1) += wall_velocity * coeff;
                    }
                    if (solver.u_face_flags(i + 1, j - 1, Flag::Open)) // Bottom Right
                    {
                        solver.u_coeff(i + 1, j - 1) += coeff;
                        solver.u_scr(i + 1, j - 1) += wall_velocity * coeff;
                    }
                }
                if (j != solver.ny)
                {
                    if (solver.u_face_flags(i, j, Flag::Open)) // Top Left
                    {
                        solver.u_coeff(i, j) += coeff;
                        solver.u_scr(i, j) += wall_velocity * coeff;
                    }
                    if (solver.u_face_flags(i + 1, j, Flag::Open)) // Top Right
                    {
                        solver.u_coeff(i + 1, j) += coeff;
                        solver.u_scr(i + 1, j) += wall_velocity * coeff;
                    }

                }
            }
        }
    }
}
