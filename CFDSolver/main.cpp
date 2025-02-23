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

        // Solve for max n itterations

        auto solve_thread = std::jthread(&SolverStaggeredIMEXTemp::solve, &s, 100000);

        interface.RenderResults();
        if (s.is_solving)
        {
            s.is_solving = false;
            interface.RenderResults();
        }
    }
	
	return 0;
}