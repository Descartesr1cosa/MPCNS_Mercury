#include "Z0_Tests.h"
#include "Z0_TestCommon.h"

#include "0_basic/MPI_WRAPPER.h"
#include "1_grid/1_MPCNS_Grid.h"
#include "2_topology/Topology.h"
#include "mpcns/compat/legacy_mesh_partition_adapter.hpp"
#include "mpcns/mesh/mesh_partition.hpp"

#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace
{
using mpcns::mesh::MeshPartition;
using mpcns::mesh::Patch;
using mpcns::mesh::PatchKind;

struct Summary
{
    long long checked = 0;
    long long failed = 0;
};

const char *kind_name(PatchKind kind)
{
    switch (kind)
    {
    case PatchKind::Inner: return "Inner";
    case PatchKind::Parallel: return "Parallel";
    case PatchKind::Physical: return "Physical";
    }
    return "Unknown";
}

PatchKind mesh_kind(TOPO::PatchKind kind)
{
    switch (kind)
    {
    case TOPO::PatchKind::Inner: return PatchKind::Inner;
    case TOPO::PatchKind::Parallel: return PatchKind::Parallel;
    case TOPO::PatchKind::Physical: return PatchKind::Physical;
    }
    return PatchKind::Physical;
}

template <class Value>
std::string value_string(const Value &value)
{
    std::ostringstream out;
    out << value;
    return out.str();
}

template <>
std::string value_string<bool>(const bool &value)
{
    return value ? "true" : "false";
}

template <>
std::string value_string<PatchKind>(const PatchKind &value)
{
    return kind_name(value);
}

struct Context
{
    int rank;
    const char *type;
    PatchKind kind;
    std::size_t index;
    int this_block;
    int neighbor_block;
};

template <class Expected, class Actual>
void check_equal(Summary &summary, const Context &ctx, const char *field,
                 const Expected &expected, const Actual &actual)
{
    ++summary.checked;
    if (expected == actual)
        return;
    ++summary.failed;
    std::cerr << "[MeshPartition][rank " << ctx.rank << "][" << ctx.type
              << "][" << kind_name(ctx.kind) << "][patch " << ctx.index
              << "][this " << ctx.this_block << "][neighbor " << ctx.neighbor_block
              << "] field=" << field << " expected=" << value_string(expected)
              << " actual=" << value_string(actual) << "\n";
}

void check_box(Summary &summary, const Context &ctx, const char *prefix,
               const Box3 &expected, const mpcns::core::LogicalBox &actual)
{
    const int e[6] = {expected.lo.i, expected.lo.j, expected.lo.k,
                      expected.hi.i, expected.hi.j, expected.hi.k};
    const int a[6] = {actual.lo.i, actual.lo.j, actual.lo.k,
                      actual.hi.i, actual.hi.j, actual.hi.k};
    const char *suffix[6] = {".lo.i", ".lo.j", ".lo.k", ".hi.i", ".hi.j", ".hi.k"};
    for (int i = 0; i < 6; ++i)
        check_equal(summary, ctx, (std::string(prefix) + suffix[i]).c_str(), e[i], a[i]);
}

void check_transform(Summary &summary, const Context &ctx,
                     const TOPO::IndexTransform &expected, const Patch &actual)
{
    check_equal(summary, ctx, "transform.has_value", true, actual.transform.has_value());
    if (!actual.transform)
        return;
    for (int axis = 0; axis < 3; ++axis)
    {
        check_equal(summary, ctx, ("transform.permutation[" + std::to_string(axis) + "]").c_str(),
                    expected.perm[axis], actual.transform->permutation[axis]);
        check_equal(summary, ctx, ("transform.sign[" + std::to_string(axis) + "]").c_str(),
                    expected.sign[axis], actual.transform->sign[axis]);
    }
    check_equal(summary, ctx, "transform.offset.i", expected.offset.i, actual.transform->offset.i);
    check_equal(summary, ctx, "transform.offset.j", expected.offset.j, actual.transform->offset.j);
    check_equal(summary, ctx, "transform.offset.k", expected.offset.k, actual.transform->offset.k);
}

template <class LegacyPatch>
void check_neighbor_common(Summary &summary, int rank, const char *type,
                           std::size_t index, const LegacyPatch &expected,
                           const Patch &actual)
{
    const Context ctx{rank, type, mesh_kind(expected.kind), index,
                      expected.this_block, expected.nb_block};
    check_equal(summary, ctx, "kind", mesh_kind(expected.kind), actual.kind);
    check_equal(summary, ctx, "this.rank", expected.this_rank, actual.this_block.rank);
    check_equal(summary, ctx, "this.local_block", expected.this_block, actual.this_block.local_block);
    check_equal(summary, ctx, "neighbor.rank", expected.nb_rank, actual.neighbor_block.rank);
    check_equal(summary, ctx, "neighbor.local_block", expected.nb_block, actual.neighbor_block.local_block);
    check_equal(summary, ctx, "this_block_name", expected.this_block_name, actual.this_block_name);
    check_equal(summary, ctx, "neighbor_block_name", expected.nb_block_name, actual.neighbor_block_name);
    check_equal(summary, ctx, "has_neighbor", true, actual.has_neighbor);
    check_equal(summary, ctx, "is_coupling", expected.is_coupling, actual.is_coupling);
    check_box(summary, ctx, "this_node_box", expected.this_box_node, actual.this_node_box);
    check_box(summary, ctx, "neighbor_node_box", expected.nb_box_node, actual.neighbor_node_box);
    check_transform(summary, ctx, expected.trans, actual);
}

void check_face(Summary &summary, int rank, std::size_t index,
                const TOPO::InterfacePatch &expected, const Patch &actual)
{
    check_neighbor_common(summary, rank, "face", index, expected, actual);
    const Context ctx{rank, "face", mesh_kind(expected.kind), index,
                      expected.this_block, expected.nb_block};
    check_equal(summary, ctx, "direction", expected.direction, actual.direction);
    check_equal(summary, ctx, "neighbor_direction", expected.nb_direction,
                actual.neighbor_direction);
}

template <class LegacyPatch>
void check_corner(Summary &summary, int rank, const char *type, std::size_t index,
                  const LegacyPatch &expected, const Patch &actual, int dimensions)
{
    check_neighbor_common(summary, rank, type, index, expected, actual);
    const Context ctx{rank, type, mesh_kind(expected.kind), index,
                      expected.this_block, expected.nb_block};
    check_equal(summary, ctx, "dir1", expected.dir1, actual.directions[0]);
    check_equal(summary, ctx, "dir2", expected.dir2, actual.directions[1]);
    if (dimensions == 3)
    {
        if constexpr (requires { expected.dir3; })
            check_equal(summary, ctx, "dir3", expected.dir3, actual.directions[2]);
    }
}

template <class LegacyPatch>
void check_physical_common(Summary &summary, int rank, const char *type,
                           std::size_t index, const LegacyPatch &expected,
                           const Patch &actual)
{
    const Context ctx{rank, type, PatchKind::Physical, index, expected.this_block, -1};
    check_equal(summary, ctx, "kind", PatchKind::Physical, actual.kind);
    check_equal(summary, ctx, "this.rank", expected.this_rank, actual.this_block.rank);
    check_equal(summary, ctx, "this.local_block", expected.this_block, actual.this_block.local_block);
    check_equal(summary, ctx, "this_block_name", expected.this_block_name, actual.this_block_name);
    check_box(summary, ctx, "this_node_box", expected.this_box_node, actual.this_node_box);
    check_equal(summary, ctx, "has_neighbor", false, actual.has_neighbor);
    check_equal(summary, ctx, "neighbor.valid", false, actual.neighbor_block.valid());
    check_equal(summary, ctx, "neighbor_block_name.empty", true, actual.neighbor_block_name.empty());
    check_equal(summary, ctx, "transform.has_value", false, actual.transform.has_value());
}

void check_physical_face(Summary &summary, int rank, std::size_t index,
                         const TOPO::PhysicalPatch &expected, const Patch &actual)
{
    check_physical_common(summary, rank, "face", index, expected, actual);
    const Context ctx{rank, "face", PatchKind::Physical, index, expected.this_block, -1};
    check_equal(summary, ctx, "direction", expected.direction, actual.direction);
    check_equal(summary, ctx, "boundary_id", expected.bc_id, actual.boundary_id);
    check_equal(summary, ctx, "boundary_name", expected.bc_name, actual.boundary_name);
}

template <class LegacyPatch>
void check_physical_corner(Summary &summary, int rank, const char *type,
                           std::size_t index, const LegacyPatch &expected,
                           const Patch &actual, int dimensions)
{
    check_physical_common(summary, rank, type, index, expected, actual);
    const Context ctx{rank, type, PatchKind::Physical, index, expected.this_block, -1};
    check_equal(summary, ctx, "is_coupling", expected.is_coupling, actual.is_coupling);
    check_equal(summary, ctx, "dir1", expected.dir1, actual.directions[0]);
    check_equal(summary, ctx, "dir2", expected.dir2, actual.directions[1]);
    if (dimensions == 3)
    {
        if constexpr (requires { expected.dir3; })
            check_equal(summary, ctx, "dir3", expected.dir3, actual.directions[2]);
    }
}

bool rejected(const std::function<void()> &operation)
{
    try { operation(); }
    catch (const std::invalid_argument &) { return true; }
    return false;
}

Patch neighbor_patch(PatchKind kind, int this_rank, int this_block,
                     int neighbor_rank, int neighbor_block)
{
    Patch patch;
    patch.kind = kind;
    patch.this_block = {this_rank, this_block};
    patch.neighbor_block = {neighbor_rank, neighbor_block};
    patch.has_neighbor = true;
    patch.transform = mpcns::mesh::IndexTransform{{0, 1, 2}, {1, 1, 1}, {0, 0, 0}};
    return patch;
}

Summary check_constructor_contracts()
{
    using mpcns::mesh::LocalBlock;
    Summary summary;
    const std::vector<LocalBlock> one{{{0, 0}, "a", {{0, 0, 0}, {2, 2, 2}}}};
    const std::vector<LocalBlock> two{{{0, 0}, "a", {{0, 0, 0}, {2, 2, 2}}},
                                      {{0, 1}, "b", {{0, 0, 0}, {2, 2, 2}}}};
    auto expect = [&](const char *name, bool value)
    {
        const Context ctx{0, "contract", PatchKind::Physical, 0, 0, -1};
        check_equal(summary, ctx, name, true, value);
    };

    expect("illegal rank", rejected([&] { MeshPartition(-1, 2, one, {}, {}, {}); }));
    expect("illegal size", rejected([&] { MeshPartition(0, 0, one, {}, {}, {}); }));
    auto wrong_order = two; wrong_order[1].address.local_block = 2;
    expect("illegal local block order", rejected([&] { MeshPartition(0, 2, wrong_order, {}, {}, {}); }));

    Patch inner = neighbor_patch(PatchKind::Inner, 0, 0, 0, 1);
    MeshPartition(0, 2, two, {inner}, {}, {});
    Patch parallel = neighbor_patch(PatchKind::Parallel, 0, 0, 1, 0);
    MeshPartition(0, 2, one, {parallel, parallel}, {}, {}); // duplicate/periodic adjacency is legal
    Patch physical; physical.kind = PatchKind::Physical; physical.this_block = {0, 0};
    MeshPartition(0, 2, one, {physical}, {}, {});

    Patch bad = inner; bad.this_block.rank = 1;
    expect("patch this rank", rejected([&] { MeshPartition(0, 2, two, {bad}, {}, {}); }));
    bad = inner; bad.this_block.local_block = 2;
    expect("patch this block", rejected([&] { MeshPartition(0, 2, two, {bad}, {}, {}); }));
    bad = parallel; bad.neighbor_block.local_block = -1;
    expect("incomplete neighbor", rejected([&] { MeshPartition(0, 2, one, {bad}, {}, {}); }));
    bad = physical; bad.has_neighbor = true; bad.neighbor_block = {1, 0};
    expect("physical neighbor", rejected([&] { MeshPartition(0, 2, one, {bad}, {}, {}); }));
    bad = inner; bad.has_neighbor = false; bad.neighbor_block = {};
    expect("nonphysical missing neighbor", rejected([&] { MeshPartition(0, 2, two, {bad}, {}, {}); }));
    bad = parallel; bad.transform.reset();
    expect("neighbor missing transform", rejected([&] { MeshPartition(0, 2, one, {bad}, {}, {}); }));
    bad = inner; bad.neighbor_block.rank = 1;
    expect("inner rank relation", rejected([&] { MeshPartition(0, 2, two, {bad}, {}, {}); }));
    bad = parallel; bad.neighbor_block.rank = 0;
    expect("parallel rank relation", rejected([&] { MeshPartition(0, 2, one, {bad}, {}, {}); }));
    return summary;
}
} // namespace

namespace Z0_TEST
{
bool RunMeshPartitionTests(const Grid &grid, const TOPO::Topology &topology, int myid)
{
    int size = 0;
    PARALLEL::mpi_size(&size);
    const MeshPartition partition =
        mpcns::compat::LegacyMeshPartitionAdapter::build(grid, topology, myid, size);
    Summary summary = check_constructor_contracts();

    const Context block_context{myid, "block", PatchKind::Physical, 0, 0, -1};
    check_equal(summary, block_context, "partition.rank", myid, partition.rank());
    check_equal(summary, block_context, "partition.size", size, partition.size());
    check_equal(summary, block_context, "local block count",
                static_cast<std::size_t>(grid.nblock), partition.local_blocks().size());
    Grid &readable = const_cast<Grid &>(grid);
    for (int ib = 0; ib < grid.nblock; ++ib)
    {
        const Block &expected = readable.grids(ib);
        const auto &actual = partition.local_blocks()[static_cast<std::size_t>(ib)];
        Context ctx{myid, "block", PatchKind::Physical, static_cast<std::size_t>(ib), ib, -1};
        check_equal(summary, ctx, "address.rank", myid, actual.address.rank);
        check_equal(summary, ctx, "address.local_block", ib, actual.address.local_block);
        check_equal(summary, ctx, "name", expected.block_name, actual.name);
        const Box3 box{{0, 0, 0}, {expected.mx + 1, expected.my + 1, expected.mz + 1}};
        check_box(summary, ctx, "node_extent", box, actual.node_extent);
    }

    std::size_t out = 0;
    for (std::size_t i = 0; i < topology.inner_patches.size(); ++i)
        check_face(summary, myid, i, topology.inner_patches[i], partition.face_patches()[out++]);
    for (std::size_t i = 0; i < topology.parallel_patches.size(); ++i)
        check_face(summary, myid, i, topology.parallel_patches[i], partition.face_patches()[out++]);
    for (std::size_t i = 0; i < topology.physical_patches.size(); ++i)
        check_physical_face(summary, myid, i, topology.physical_patches[i], partition.face_patches()[out++]);
    check_equal(summary, block_context, "face patch count", out, partition.face_patches().size());

    out = 0;
    for (std::size_t i = 0; i < topology.inner_edge_patches.size(); ++i)
        check_corner(summary, myid, "edge", i, topology.inner_edge_patches[i], partition.edge_patches()[out++], 2);
    for (std::size_t i = 0; i < topology.parallel_edge_patches.size(); ++i)
        check_corner(summary, myid, "edge", i, topology.parallel_edge_patches[i], partition.edge_patches()[out++], 2);
    for (std::size_t i = 0; i < topology.physical_edge_patches.size(); ++i)
        check_physical_corner(summary, myid, "edge", i, topology.physical_edge_patches[i], partition.edge_patches()[out++], 2);
    check_equal(summary, block_context, "edge patch count", out, partition.edge_patches().size());

    out = 0;
    for (std::size_t i = 0; i < topology.inner_vertex_patches.size(); ++i)
        check_corner(summary, myid, "vertex", i, topology.inner_vertex_patches[i], partition.vertex_patches()[out++], 3);
    for (std::size_t i = 0; i < topology.parallel_vertex_patches.size(); ++i)
        check_corner(summary, myid, "vertex", i, topology.parallel_vertex_patches[i], partition.vertex_patches()[out++], 3);
    for (std::size_t i = 0; i < topology.physical_vertex_patches.size(); ++i)
        check_physical_corner(summary, myid, "vertex", i, topology.physical_vertex_patches[i], partition.vertex_patches()[out++], 3);
    check_equal(summary, block_context, "vertex patch count", out, partition.vertex_patches().size());

    const long long global_checked = Z0_TEST::global_sum(summary.checked);
    const long long global_failed = Z0_TEST::global_sum(summary.failed);
    if (myid == 0)
        std::cout << "  [" << (global_failed == 0 ? "PASS" : "FAIL")
                  << "] MeshPartition legacy adapter exact equivalence checked="
                  << global_checked << " failed=" << global_failed << "\n";
    return summary.failed == 0;
}
} // namespace Z0_TEST
