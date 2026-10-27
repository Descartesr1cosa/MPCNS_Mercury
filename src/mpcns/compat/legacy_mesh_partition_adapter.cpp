#include "mpcns/compat/legacy_mesh_partition_adapter.hpp"

#include "1_grid/1_MPCNS_Grid.h"
#include "2_topology/Topology.h"

namespace
{
using mpcns::core::LogicalBox;
using mpcns::core::LogicalIndex;
using mpcns::mesh::BlockAddress;
using mpcns::mesh::IndexTransform;
using mpcns::mesh::Patch;
using mpcns::mesh::PatchKind;

LogicalBox logical_box(const Box3 &box)
{
    return {{box.lo.i, box.lo.j, box.lo.k}, {box.hi.i, box.hi.j, box.hi.k}};
}
IndexTransform transform(const TOPO::IndexTransform &legacy)
{
    IndexTransform value;
    for (int axis = 0; axis < 3; ++axis)
    {
        value.permutation[axis] = legacy.perm[axis];
        value.sign[axis] = legacy.sign[axis];
    }
    value.offset = {legacy.offset.i, legacy.offset.j, legacy.offset.k};
    return value;
}
PatchKind kind(TOPO::PatchKind legacy)
{
    return legacy == TOPO::PatchKind::Inner ? PatchKind::Inner :
           legacy == TOPO::PatchKind::Parallel ? PatchKind::Parallel : PatchKind::Physical;
}
Patch face(const TOPO::InterfacePatch &legacy)
{
    Patch value;
    value.kind = kind(legacy.kind);
    value.this_block = {legacy.this_rank, legacy.this_block};
    value.neighbor_block = {legacy.nb_rank, legacy.nb_block};
    value.this_block_name = legacy.this_block_name;
    value.neighbor_block_name = legacy.nb_block_name;
    value.this_node_box = logical_box(legacy.this_box_node);
    value.neighbor_node_box = logical_box(legacy.nb_box_node);
    value.transform = transform(legacy.trans);
    value.has_neighbor = true;
    value.is_coupling = legacy.is_coupling;
    value.direction = legacy.direction;
    value.neighbor_direction = legacy.nb_direction;
    return value;
}
template <class LegacyPatch>
Patch corner(const LegacyPatch &legacy, PatchKind patch_kind, int directions)
{
    Patch value;
    value.kind = patch_kind;
    value.this_block = {legacy.this_rank, legacy.this_block};
    if (patch_kind != PatchKind::Physical)
    {
        value.neighbor_block = {legacy.nb_rank, legacy.nb_block};
        value.neighbor_block_name = legacy.nb_block_name;
        value.neighbor_node_box = logical_box(legacy.nb_box_node);
        value.transform = transform(legacy.trans);
        value.has_neighbor = true;
    }
    value.this_block_name = legacy.this_block_name;
    value.this_node_box = logical_box(legacy.this_box_node);
    value.is_coupling = legacy.is_coupling;
    value.directions[0] = legacy.dir1;
    value.directions[1] = legacy.dir2;
    if constexpr (requires { legacy.dir3; }) { if (directions == 3) value.directions[2] = legacy.dir3; }
    return value;
}
}
namespace mpcns::compat
{
mesh::MeshPartition LegacyMeshPartitionAdapter::build(const Grid &grid, const TOPO::Topology &topology,
                                                      int rank, int size)
{
    // Array1D predates const indexing.  This scoped cast is a read-only bridge;
    // no Grid member is written and the Mesh API itself remains legacy-free.
    Grid &legacy_grid = const_cast<Grid &>(grid);
    std::vector<mesh::LocalBlock> blocks;
    blocks.reserve(grid.nblock);
    for (int ib = 0; ib < grid.nblock; ++ib)
    {
        const Block &block = legacy_grid.grids(ib);
        blocks.push_back({{rank, ib}, block.block_name,
                          {{0, 0, 0}, {block.mx + 1, block.my + 1, block.mz + 1}}});
    }
    std::vector<Patch> faces, edges, vertices;
    for (const auto &p : topology.inner_patches) faces.push_back(face(p));
    for (const auto &p : topology.parallel_patches) faces.push_back(face(p));
    for (const auto &p : topology.physical_patches)
    {
        Patch value;
        value.kind = PatchKind::Physical;
        value.this_block = {p.this_rank, p.this_block};
        value.this_block_name = p.this_block_name;
        value.this_node_box = logical_box(p.this_box_node);
        value.direction = p.direction;
        value.boundary_id = p.bc_id;
        value.boundary_name = p.bc_name;
        faces.push_back(std::move(value));
    }
    for (const auto &p : topology.inner_edge_patches) edges.push_back(corner(p, PatchKind::Inner, 2));
    for (const auto &p : topology.parallel_edge_patches) edges.push_back(corner(p, PatchKind::Parallel, 2));
    for (const auto &p : topology.physical_edge_patches) edges.push_back(corner(p, PatchKind::Physical, 2));
    for (const auto &p : topology.inner_vertex_patches) vertices.push_back(corner(p, PatchKind::Inner, 3));
    for (const auto &p : topology.parallel_vertex_patches) vertices.push_back(corner(p, PatchKind::Parallel, 3));
    for (const auto &p : topology.physical_vertex_patches) vertices.push_back(corner(p, PatchKind::Physical, 3));
    return mesh::MeshPartition(rank, size, std::move(blocks), std::move(faces),
                               std::move(edges), std::move(vertices));
}
} // namespace mpcns::compat
