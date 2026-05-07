#include "Z1_Solver.h"

#include "Z1_Boundary.h"
#include "Z1_Diagnostics.h"

#include "0_basic/MPI_WRAPPER.h"

#include <iostream>

void Z1_Solver::PrepareStep_()
{
    boundary_->PrepareStepBoundary();
}

void Z1_Solver::UpdateFields_()
{
    // no-op template hook.
    // New physics modules should update fields here, for example:
    // U^{n+1} = U^n + dt * RHS.
}

void Z1_Solver::ApplyBoundaryAndSync_()
{
    boundary_->ApplyAllPhysicalBoundaries();
    boundary_->SyncGroup("U");
    boundary_->SyncGroup("V");
    boundary_->SyncGroup("rho");
    boundary_->SyncGroup("Bface");
    boundary_->SyncGroup("Eedge");
}

void Z1_Solver::Diagnostics_()
{
    int myid = 0;
    PARALLEL::mpi_rank(&myid);
    if (myid == 0)
        Z1::PrintStepDiagnostics(*field_, control_.step, control_.time, control_.dt, std::cout);
}

void Z1_Solver::Output_()
{
    // no-op template hook.
    // A real solver can configure IOModule here for restart/plot output.
}

void Z1_Solver::AdvanceRunState_()
{
    control_.advance();
}
