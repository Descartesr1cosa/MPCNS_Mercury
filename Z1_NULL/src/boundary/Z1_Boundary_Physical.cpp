#include "Z1_Boundary.h"

void Z1_Boundary::ApplyPhysicalBoundary(const std::string &group)
{
    // no-op template hook.
    // New physics modules should dispatch physical boundary handlers by group here.
    (void)group;
}

void Z1_Boundary::ApplyAllPhysicalBoundaries()
{
    ApplyPhysicalBoundary("U");
    ApplyPhysicalBoundary("V");
    ApplyPhysicalBoundary("rho");
    ApplyPhysicalBoundary("Bface");
    ApplyPhysicalBoundary("Eedge");
}
