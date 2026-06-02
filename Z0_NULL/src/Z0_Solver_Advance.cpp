#include "Z0_Solver.h"

#include "Z0_Boundary.h"
#include "Z0_Diagnostics.h"

#include "0_basic/MPI_WRAPPER.h"
#include <iostream>

void Z0_Solver::PrepareStep_()
{
    boundary_->PrepareStepBoundary();
}

void Z0_Solver::UpdateFields_()
{
    // no-op template hook.
    // New physics modules should update fields here.
}

void Z0_Solver::ApplyBoundaryAndSync_()
{
    boundary_->ApplyAllPhysicalBoundaries();
    boundary_->SyncGroup("null_phi");
}

void Z0_Solver::Diagnostics_()
{
    int myid = 0;
    PARALLEL::mpi_rank(&myid);
    if (myid == 0)
        Z0::PrintStepDiagnostics(*field_, control_.step, control_.time, control_.dt, std::cout);
}

void Z0_Solver::Output_()
{
    // no-op template hook for output.
}

void Z0_Solver::AdvanceRunState_()
{
    control_.advance();
}
