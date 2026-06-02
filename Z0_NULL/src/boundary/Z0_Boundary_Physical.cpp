#include "Z0_Boundary.h"

void Z0_Boundary::ApplyPhysicalBoundary(const std::string &group)
{
    // no-op template hook.
    // New physics modules should dispatch physical boundary handlers by group here.
    (void)group;
}

void Z0_Boundary::ApplyAllPhysicalBoundaries()
{
    ApplyPhysicalBoundary("null_phi");
}
