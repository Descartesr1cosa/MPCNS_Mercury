#include "Z1_Boundary.h"

Z1_Boundary::Z1_Boundary(Grid *grid,
                         Field *field,
                         Halo *halo,
                         TOPO::Topology *topology,
                         TOPO::TopologyEquiv *topology_equiv,
                         Param *param)
    : grid_(grid),
      field_(field),
      halo_(halo),
      topology_(topology),
      topology_equiv_(topology_equiv),
      param_(param)
{
}

void Z1_Boundary::PrepareStepBoundary()
{
    ApplyAllPhysicalBoundaries();
    ApplyCoupling();
}

void Z1_Boundary::FinishStepBoundary()
{
    SyncAllRegistered();
}
