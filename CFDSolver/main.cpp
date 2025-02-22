#include "Solver.hpp"
#include "Interface.hpp"
#include "Solver_Staggered_IMEX.hpp"
#include "Solver_Staggered_IMEX_Temp.hpp"

int main()
{
    Interface interface;

    while (!WindowShouldClose())
    {
        // Define grid

        interface.RenderModel();

        SolverStaggeredIMEXTemp s(interface);

        // Set Boundaries

        auto& b1 = s.inlet_boundary_conditions.emplace_back(0, 0);
        b1.inlet_temperature = 20;
        b1.inlet_velocity = 1.0;
        for (int j = 0; j < s.ny / 4; ++j)
            b1.boundary_faces.push_back({ 0, j });

        auto& b2 = s.inlet_boundary_conditions.emplace_back(0, 1);
        b2.inlet_temperature = 40;
        b2.inlet_velocity = 1.0;
        for (int i = s.ny / 3; i < s.ny / 2; ++i)
            b2.boundary_faces.push_back({ i, 0 });

        auto& b3 = s.outlet_boundary_condition.emplace_back(1, 0);
        for (int j = 0; j < s.ny / 4; ++j)
            b3.boundary_faces.push_back({ s.nx, j });

        // Solve for max n itterations

        auto solve_thread = std::jthread(&SolverStaggeredIMEXTemp::solve, &s, 10000);

        interface.RenderResults();
        if (s.is_solving)
        {
            s.is_solving = false;
            interface.RenderResults();
        }
    }
	
	return 0;
}