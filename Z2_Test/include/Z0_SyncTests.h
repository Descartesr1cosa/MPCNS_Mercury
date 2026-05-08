#pragma once

#include "Z0_Diagnostics.h"

#include <iosfwd>

class Field;
class Grid;
class Halo;

namespace TOPO
{
    struct Topology;
    struct TopologyEquiv;
}

namespace Z0
{
    TestResult test_field_extents(Field &fields, Grid &grid, int dimension, int my_rank, std::ostream &os);
    TestResult test_component_halo(Field &fields, Halo &halo, const TOPO::Topology &topology, int my_rank, std::ostream &os);
    TestResult test_edge_1form_triplet_halo(Field &fields, Halo &halo, const TOPO::Topology &topology, int my_rank, std::ostream &os);
    TestResult test_face_2form_triplet_halo(Field &fields, Halo &halo, const TOPO::Topology &topology, int my_rank, std::ostream &os);
    TestResult test_owner_alias_sync(Field &fields, Halo &halo, const TOPO::TopologyEquiv &equiv, int my_rank, std::ostream &os);
    TestResult test_sync_group_order(Field &fields, const TOPO::Topology &topology, int my_rank, std::ostream &os);
}
