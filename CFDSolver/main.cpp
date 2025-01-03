#include "Solver.hpp"
#include "SolverFaceVelocity.hpp"
#include "SolverFaceVelocityV2.hpp"

int main()
{
	Solver s(5, 5, 0.01, 0.01);
	s.solve(20000);
	
	//SolverFaceVelocity s(50, 50, 0.01, 0.01);
	//s.solve(20000);

	return 0;
}