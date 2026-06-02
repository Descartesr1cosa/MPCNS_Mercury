#include "Z0_Solver.h"

#include "Z0_Initial.h"

void Z0_Solver::Initialize()
{
    Setup_();
    if (initialized_)
        return;
    Z0::InitializeFields(*field_, *param_, *boundary_);
    initialized_ = true;
}

void Z0_Solver::Advance()
{
    Setup_();
    Initialize();

    while (!control_.should_stop())
    {
        PrepareStep_();
        ComputeTimeStep_();
        ComputeRHS_();
        UpdateFields_();
        ApplyBoundaryAndSync_();
        Diagnostics_();
        Output_();
        AdvanceRunState_();
    }
}
