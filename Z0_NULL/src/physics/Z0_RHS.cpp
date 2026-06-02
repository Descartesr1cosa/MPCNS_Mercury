#include "Z0_Solver.h"

void Z0_Solver::ComputeRHS_()
{
    ComputeFlux_();
    ComputeSource_();
    // no-op template hook.
    // New physics modules should assemble geometry operators, fluxes, sources,
    // and RHS fields here.
}
