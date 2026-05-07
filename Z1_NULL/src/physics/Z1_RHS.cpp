#include "Z1_Solver.h"

void Z1_Solver::ComputeRHS_()
{
    ComputeFlux_();
    ComputeSource_();
    // no-op template hook.
    // New physics modules should assemble geometry operators, fluxes, sources,
    // and RHS fields here.
}
