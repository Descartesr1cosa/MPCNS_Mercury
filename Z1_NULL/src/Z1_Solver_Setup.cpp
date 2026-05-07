#include "Z1_Solver.h"

#include "Z1_Control.h"

void Z1_Solver::Setup_()
{
    if (setup_done_)
        return;
    control_ = Z1::LoadControl(*param_);
    setup_done_ = true;
}
