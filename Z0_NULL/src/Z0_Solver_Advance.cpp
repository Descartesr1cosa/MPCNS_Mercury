#include "Z0_Solver.h"

#include "Z0_Boundary.h"
#include "Z0_Diagnostics.h"

#include "0_basic/MPI_WRAPPER.h"
#include <cstdlib>
#include <iostream>
#include <string>

namespace
{
    bool env_enabled(const char *name)
    {
        const char *value = std::getenv(name);
        if (!value)
            return false;
        const std::string s(value);
        return s == "1" || s == "true" || s == "TRUE" || s == "on" || s == "ON";
    }

    bool skip_solver_sync()
    {
        return env_enabled("Z0_OWNER_ONLY") || env_enabled("Z0_TEST_DRIVEN_SYNC");
    }
}

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
    if (!skip_solver_sync())
        boundary_->SyncAllRegistered();
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
    io_.WriteTecplotBinFile(control_.step, control_.time);
}

void Z0_Solver::AdvanceRunState_()
{
    control_.advance();
}
