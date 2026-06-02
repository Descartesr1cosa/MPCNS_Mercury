#include "Z0_Solver.h"

#include "Z0_Control.h"

void Z0_Solver::Setup_()
{
    if (setup_done_)
        return;
    control_ = Z0::LoadControl(*param_);
    setup_done_ = true;
}
