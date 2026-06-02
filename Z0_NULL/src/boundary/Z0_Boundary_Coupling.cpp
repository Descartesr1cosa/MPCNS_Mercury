#include "Z0_Boundary.h"

void Z0_Boundary::ApplyCoupling()
{
    // no-op template hook.
    // A real multi-physics solver should call coupling_trans_* and then consume
    // coupling buffers through its boundary/coupling handlers here.
}
