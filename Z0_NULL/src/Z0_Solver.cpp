#include "Z0_Solver.h"

Z0_Solver::Z0_Solver(Grid *grid,
                     Field *field,
                     Halo *halo,
                     TOPO::Topology *topology,
                     TOPO::Topology *topology_equiv,
                     Z0_Boundary *boundary,
                     Param *param)
    : grid_(grid),
      field_(field),
      halo_(halo),
      topology_(topology),
      topology_equiv_(topology_equiv),
      boundary_(boundary),
      param_(param)
{
}
