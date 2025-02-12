#include "Solver.hpp"
#include "SolverFaceVelocity.hpp"
#include "SolverFaceVelocityV2.hpp"

int main()
{
	SolverStaggered s(20, 20, 0.005, 0.005);
	s.solve(100000);
	
	return 0;
}