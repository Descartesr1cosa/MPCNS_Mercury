#pragma once

#include "2_topology/TopologyTypes.h"

class Grid;

namespace TOPO_DEBUG
{
    void dump_topology_summary(const TOPO::Topology &topo, int my_rank);

    void validate_topology_or_abort(TOPO::Topology &topo,
                                    Grid &grid,
                                    int my_rank,
                                    int dimension);
}