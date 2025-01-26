#include "Solver.hpp"
#include "SolverFaceVelocity.hpp"
#include "SolverFaceVelocityV2.hpp"

int main()
{
	SolverStaggered s(100, 100, 0.005, 0.005);
	s.solve(10000);
	
	return 0;
}