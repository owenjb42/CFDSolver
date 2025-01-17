#include "Solver.hpp"
#include "SolverFaceVelocity.hpp"
#include "SolverFaceVelocityV2.hpp"

int main()
{
	SolverStaggered s(100, 100, 0.01, 0.01);
	s.solve(4000);
	
	return 0;
}