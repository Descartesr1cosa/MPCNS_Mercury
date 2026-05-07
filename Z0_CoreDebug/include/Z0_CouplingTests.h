#pragma once

#include "Z0_Diagnostics.h"

#include <iosfwd>

class Field;
class Halo;

namespace TOPO
{
    struct Topology;
}

namespace Z0
{
    TestResult test_coupling_cell_scalar(Field &fields, Halo &halo, const TOPO::Topology &topology, int dimension, int my_rank, std::ostream &os);
    TestResult test_coupling_edge_1form(Field &fields, Halo &halo, const TOPO::Topology &topology, int dimension, int my_rank, std::ostream &os);
    TestResult test_coupling_face_2form(Field &fields, Halo &halo, const TOPO::Topology &topology, int dimension, int my_rank, std::ostream &os);
}
