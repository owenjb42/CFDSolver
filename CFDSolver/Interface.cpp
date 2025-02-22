#include "Interface.hpp"
#include "Solver_Staggered_IMEX_Temp.hpp"

void Interface::SetData(SolverStaggeredIMEXTemp& solver)
{
    std::lock_guard lk(mutex);

    u = solver.cell_u;
    v = solver.cell_v;
    p = solver.p;
    t = solver.t;
    divergence = solver.divergence;
    u_face = solver.u;
    v_face = solver.v;

    needs_update = true;
}