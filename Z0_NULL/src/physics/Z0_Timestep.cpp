#include "Z0_Solver.h"

#include "Z0_Const.h"

void Z0_Solver::ComputeTimeStep_()
{
    // Fixed placeholder dt. A real solver should compute CFL/physics limits here.
    control_.dt = Z0::DefaultDt;
}
