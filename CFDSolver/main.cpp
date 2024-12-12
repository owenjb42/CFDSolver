#include "Solver.hpp"

int main()
{
	Solver s(50, 50, 1.0, 1.0);
	s.solve(1000);

	return 0;
}