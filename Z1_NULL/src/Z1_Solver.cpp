#include "Z1_Solver.h"

Z1_Solver::Z1_Solver(Grid *grid,
                     Field *field,
                     Halo *halo,
                     TOPO::Topology *topology,
                     TOPO::Topology *topology_equiv,
                     Z1_Boundary *boundary,
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
