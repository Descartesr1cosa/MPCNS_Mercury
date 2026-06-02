#include "Z0_Boundary.h"

Z0_Boundary::Z0_Boundary(Grid *grid,
                         Field *field,
                         Halo *halo,
                         TOPO::Topology *topology,
                         TOPO::Topology *topology_equiv,
                         Param *param)
    : grid_(grid),
      field_(field),
      halo_(halo),
      topology_(topology),
      topology_equiv_(topology_equiv),
      param_(param)
{
}

void Z0_Boundary::PrepareStepBoundary()
{
    ApplyAllPhysicalBoundaries();
    ApplyCoupling();
}

void Z0_Boundary::FinishStepBoundary()
{
    SyncAllRegistered();
}
