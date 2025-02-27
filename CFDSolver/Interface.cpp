#include "Interface.hpp"
#include "Solver_Staggered_IMEX_Temp.hpp"

void Interface::SetData(SolverStaggeredIMEXTemp& solver)
{
    std::lock_guard lk(mutex);

    u_buffer = solver.cell_u;
    v_buffer = solver.cell_v;
    p_buffer = solver.p;
    t_buffer = solver.t;
    divergence_buffer = solver.divergence;
    u_face_buffer = solver.u;
    v_face_buffer = solver.v;

    needs_update = true;
    is_solving = true;
}

void Interface::GetDataFromBuffer()
{
    std::lock_guard lk(mutex);
    if (needs_update)
    {
        u = u_buffer;
        v = v_buffer;
        p = p_buffer;
        t = t_buffer;
        divergence = divergence_buffer;
        u_face = u_face_buffer;
        v_face = v_face_buffer;

        needs_update = false;
        recalculate_auxiliary_data = true;
    }
}