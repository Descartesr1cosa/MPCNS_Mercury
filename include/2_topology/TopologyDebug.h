#pragma once

#include "2_topology/2_MPCNS_Topology.h"

namespace TOPO_DEBUG
{
    void dump_topology_summary(const TOPO::Topology &topo, int my_rank);

    void validate_topology_or_abort(TOPO::Topology &topo,
                                    Grid &grid,
                                    int my_rank,
                                    int dimension);
}