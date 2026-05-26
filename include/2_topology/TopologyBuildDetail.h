#pragma once

#include <string>

#include "2_topology/Topology.h"

class Grid;

namespace TOPO::detail
{
    // Internal construction stages. Public callers should use build_topology().
    void build_edge_patches(Grid &grid, Topology &topology, int dimension);
    void build_vertex_patches(Grid &grid, Topology &topology, int dimension);
    void append_coupling_faces_as_physical_patches(
        Grid &grid,
        Topology &topology,
        int dimension,
        const std::string &prefix);
    void build_equivalence(
        Topology &topology,
        Grid &grid,
        int my_rank,
        int dimension);
} // namespace TOPO::detail
