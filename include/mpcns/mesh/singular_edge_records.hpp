#pragma once

#include "2_topology/Topology.h"

namespace mpcns::mesh
{

struct SingularEdgeAlias
{
    TOPO::EntityKey edge;
    int orientation = +1;
    bool owner = false;
};

struct WeightedIncidentEntity
{
    TOPO::EntityKey entity;
    double measure = 0.0;
    double weight = 0.0;
    TOPO::EntityKey source_alias{
        TOPO::EntityDim::Edge, 0, 0, 0, 0, 0, TOPO::EntityAxis::Xi};
    int source_orientation = +1;
};

} // namespace mpcns::mesh
