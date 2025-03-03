#include "Solver.hpp"
#include "Interface.hpp"
#include "Solver_Staggered_IMEX.hpp"
#include "Solver_Staggered_IMEX_Temp.hpp"

int main()
{
    Interface interface;

    while (!WindowShouldClose())
    {
        // Define Model

        interface.RenderModel();

        SolverStaggeredIMEXTemp s(interface);

        // Solve for max n itterations

        auto solve_thread = std::jthread(&SolverStaggeredIMEXTemp::solve, &s, 1000000);

        interface.RenderResults();
        if (s.is_solving) // If it is still solving stop and post process
        {
            s.is_solving = false;
            solve_thread.join();
            interface.RenderResults();
        }
    }
	
	return 0;
}