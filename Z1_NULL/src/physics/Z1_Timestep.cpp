#include "Z1_Solver.h"

#include "Z1_Const.h"

void Z1_Solver::ComputeTimeStep_()
{
    // Fixed placeholder dt. A real solver should compute CFL/physics limits here.
    control_.dt = Z1::DefaultDt;
}
