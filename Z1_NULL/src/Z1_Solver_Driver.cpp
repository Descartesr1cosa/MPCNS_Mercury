#include "Z1_Solver.h"

#include "Z1_Initial.h"

void Z1_Solver::Initialize()
{
    Setup_();
    if (initialized_)
        return;
    Z1::InitializeFields(*field_, *param_, *boundary_);
    initialized_ = true;
}

void Z1_Solver::Advance()
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
