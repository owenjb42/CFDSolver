#include "Solver.hpp"
#include "Solver_Staggered_IMEX.hpp"
#include "Solver_Staggered_IMEX_Temp.hpp"

int main()
{
	SolverStaggeredIMEXTemp s(50, 50, 0.005, 0.005);
	s.solve(20000);
	
	return 0;
}