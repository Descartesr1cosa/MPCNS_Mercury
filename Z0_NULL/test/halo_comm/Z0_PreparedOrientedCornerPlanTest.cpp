#include "Z0_Tests.h"
#include "Z0_TestCommon.h"

#include "0_basic/MPI_WRAPPER.h"
#include "1_grid/1_MPCNS_Grid.h"
#include "2_topology/Topology.h"
#include "2_topology/TopologyView.h"
#include "3_field/Field.h"
#include "4_halo/Halo.h"
#include "4_halo/detail/HaloBuildBoxMakers.h"
#include "4_halo/detail/HaloBuildTools.h"
#include "mpcns/compat/legacy_mesh_partition_adapter.hpp"
#include "mpcns/compat/legacy_oriented_corner_plan_adapter.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace
{
using namespace mpcns;

struct Summary
{
    long long checked = 0;
    long long failed = 0;
    long long reported = 0;
};

struct Context
{
    int rank = -1;
    const char *role = "plan";
    int peer = -1;
    int block = -1;
    int axis = -1;
    int tag = -1;
    std::size_t local_order = 0;
    std::size_t peer_order = 0;
    std::uint64_t generation = 0;
};

template <class Expected, class Actual>
void check_equal(Summary &summary, const Context &context,
                 const std::string &field,
                 const Expected &expected, const Actual &actual)
{
    ++summary.checked;
    if (expected == actual)
        return;
    ++summary.failed;
    if (summary.reported++ >= 20)
        return;
    std::cerr << "[PreparedOrientedCornerPlan][rank " << context.rank
              << "][group Eedge][stage Edge]"
              << "[role " << context.role << "]"
              << "[peer " << context.peer << "]"
              << "[block " << context.block << "]"
              << "[axis " << context.axis << "]"
              << "[tag " << context.tag << "]"
              << "[local_order " << context.local_order << "]"
              << "[peer_order " << context.peer_order << "]"
              << "[generation " << context.generation << "]"
              << " field=" << field << std::setprecision(17)
              << " expected=" << expected
              << " actual=" << actual << "\n";
}

StaggerLocation edge_location(int axis)
{
    const StaggerLocation locations[3] = {
        StaggerLocation::EdgeXi,
        StaggerLocation::EdgeEt,
        StaggerLocation::EdgeZe};
    return locations[axis];
}

int volume(const Box3 &box)
{
    return (box.hi.i - box.lo.i) *
           (box.hi.j - box.lo.j) *
           (box.hi.k - box.lo.k);
}

void map_index(const TOPO::IndexTransform &transform,
               int i, int j, int k,
               int &io, int &jo, int &ko)
{
    const int source[3] = {i, j, k};
    int target[3] = {0, 0, 0};
    const int offset[3] = {
        transform.offset.i,
        transform.offset.j,
        transform.offset.k};
    for (int axis = 0; axis < 3; ++axis)
        target[transform.perm[axis]] =
            transform.sign[axis] * source[axis] +
            offset[axis];
    io = target[0];
    jo = target[1];
    ko = target[2];
}

void compare_box(Summary &summary, const Context &context,
                 const std::string &prefix,
                 const Box3 &expected, const Box3 &actual)
{
    const int e[6] = {
        expected.lo.i, expected.lo.j, expected.lo.k,
        expected.hi.i, expected.hi.j, expected.hi.k};
    const int a[6] = {
        actual.lo.i, actual.lo.j, actual.lo.k,
        actual.hi.i, actual.hi.j, actual.hi.k};
    const char *suffix[6] = {
        ".lo.i", ".lo.j", ".lo.k",
        ".hi.i", ".hi.j", ".hi.k"};
    for (int index = 0; index < 6; ++index)
        check_equal(summary, context, prefix + suffix[index],
                    e[index], a[index]);
}

void compare_transform(
    Summary &summary, const Context &context,
    const TOPO::IndexTransform &expected,
    const TOPO::IndexTransform &actual)
{
    for (int axis = 0; axis < 3; ++axis)
    {
        check_equal(
            summary, context,
            "transform.permutation[" + std::to_string(axis) + "]",
            expected.perm[axis], actual.perm[axis]);
        check_equal(
            summary, context,
            "transform.sign[" + std::to_string(axis) + "]",
            expected.sign[axis], actual.sign[axis]);
    }
    check_equal(summary, context, "transform.offset.i",
                expected.offset.i, actual.offset.i);
    check_equal(summary, context, "transform.offset.j",
                expected.offset.j, actual.offset.j);
    check_equal(summary, context, "transform.offset.k",
                expected.offset.k, actual.offset.k);
}

LegacyOrientedCornerTransferSnapshot make_expected_receive(
    const TOPO::EdgePatch &patch,
    int recv_axis, int ncomp, int nghost,
    int tag, std::size_t order)
{
    Direction this_d1 =
        HALO_TOOLS::int_to_direction(patch.dir1);
    Direction this_d2 =
        HALO_TOOLS::int_to_direction(patch.dir2);
    int send_dir1 = HALO_TOOLS::map_dir_to_neighbor(
        patch.dir1, patch.trans, true);
    int send_dir2 = HALO_TOOLS::map_dir_to_neighbor(
        patch.dir2, patch.trans, false);
    if (patch.trans.perm[std::abs(patch.dir1) - 1] + 1 !=
        std::abs(send_dir1))
        std::swap(send_dir1, send_dir2);

    LegacyOrientedCornerTransferSnapshot result;
    result.role = LegacyOrientedCornerRole::Receive;
    result.recv_rank = patch.this_rank;
    result.send_rank = patch.nb_rank;
    result.recv_block = patch.this_block;
    result.send_block = patch.nb_block;
    result.recv_axis = recv_axis;
    result.send_axis = patch.trans.perm[recv_axis];
    result.sign = patch.trans.sign[recv_axis];
    result.recv_box = HALO_BOX::make_2DCorner_ghost_box(
        edge_location(recv_axis), patch.this_box_node,
        this_d1, this_d2, nghost);
    result.send_box = std::nullopt;
    result.node_box_on_send = patch.nb_box_node;
    result.dir1_send = send_dir1;
    result.dir2_send = send_dir2;
    result.trans_recv_to_send = patch.trans;
    result.peer_rank = patch.nb_rank;
    result.tag = tag;
    result.ncomp = ncomp;
    result.nghost = nghost;
    result.message_length = volume(result.recv_box) * ncomp;
    result.local_plan_order = order;
    result.peer_message_order =
        static_cast<std::size_t>(tag);
    return result;
}

std::vector<LegacyOrientedCornerTransferSnapshot>
expected_receive_plan(const TOPO::Topology &topology,
                      int ncomp, int nghost)
{
    std::vector<LegacyOrientedCornerTransferSnapshot> result;
    std::map<int, int> next_tag;
    for (const TOPO::EdgePatch &patch :
         TOPO_VIEW::edge_patches(
             topology, TOPO::PatchKind::Parallel))
    {
        if (patch.is_coupling)
            continue;
        for (int recv_axis = 0; recv_axis < 3; ++recv_axis)
        {
            const int tag = next_tag[patch.nb_rank]++;
            result.push_back(make_expected_receive(
                patch, recv_axis, ncomp, nghost,
                tag, result.size()));
        }
    }
    return result;
}

constexpr int kPackedMetaSize = 36;

std::array<int, kPackedMetaSize> pack_meta(
    const LegacyOrientedCornerTransferSnapshot &m)
{
    std::array<int, kPackedMetaSize> p{};
    int n = 0;
    const auto put_box = [&](const Box3 &box, int *values, int &at)
    {
        values[at++] = box.lo.i;
        values[at++] = box.lo.j;
        values[at++] = box.lo.k;
        values[at++] = box.hi.i;
        values[at++] = box.hi.j;
        values[at++] = box.hi.k;
    };
    p[n++] = m.recv_rank;
    p[n++] = m.send_rank;
    p[n++] = m.recv_block;
    p[n++] = m.send_block;
    p[n++] = m.recv_axis;
    p[n++] = m.send_axis;
    p[n++] = m.sign;
    put_box(m.recv_box, p.data(), n);
    put_box(m.node_box_on_send, p.data(), n);
    p[n++] = m.dir1_send;
    p[n++] = m.dir2_send;
    for (int axis = 0; axis < 3; ++axis)
        p[n++] = m.trans_recv_to_send.perm[axis];
    for (int axis = 0; axis < 3; ++axis)
        p[n++] = m.trans_recv_to_send.sign[axis];
    p[n++] = m.trans_recv_to_send.offset.i;
    p[n++] = m.trans_recv_to_send.offset.j;
    p[n++] = m.trans_recv_to_send.offset.k;
    p[n++] = m.tag;
    p[n++] = m.ncomp;
    p[n++] = m.nghost;
    p[n++] = m.message_length;
    p[n++] = static_cast<int>(m.local_plan_order);
    p[n++] = static_cast<int>(m.peer_message_order);
    if (n != kPackedMetaSize)
        throw std::logic_error("packed metadata size mismatch");
    return p;
}

LegacyOrientedCornerTransferSnapshot unpack_meta(const int *p)
{
    LegacyOrientedCornerTransferSnapshot m;
    int n = 0;
    const auto get_box = [&](Box3 &box, const int *values, int &at)
    {
        box.lo.i = values[at++];
        box.lo.j = values[at++];
        box.lo.k = values[at++];
        box.hi.i = values[at++];
        box.hi.j = values[at++];
        box.hi.k = values[at++];
    };
    m.role = LegacyOrientedCornerRole::Receive;
    m.recv_rank = p[n++];
    m.send_rank = p[n++];
    m.recv_block = p[n++];
    m.send_block = p[n++];
    m.recv_axis = p[n++];
    m.send_axis = p[n++];
    m.sign = p[n++];
    get_box(m.recv_box, p, n);
    get_box(m.node_box_on_send, p, n);
    m.dir1_send = p[n++];
    m.dir2_send = p[n++];
    for (int axis = 0; axis < 3; ++axis)
        m.trans_recv_to_send.perm[axis] = p[n++];
    for (int axis = 0; axis < 3; ++axis)
        m.trans_recv_to_send.sign[axis] = p[n++];
    m.trans_recv_to_send.offset.i = p[n++];
    m.trans_recv_to_send.offset.j = p[n++];
    m.trans_recv_to_send.offset.k = p[n++];
    m.tag = p[n++];
    m.ncomp = p[n++];
    m.nghost = p[n++];
    m.message_length = p[n++];
    m.local_plan_order = static_cast<std::size_t>(p[n++]);
    m.peer_message_order = static_cast<std::size_t>(p[n++]);
    m.peer_rank = m.send_rank;
    if (n != kPackedMetaSize)
        throw std::logic_error("unpacked metadata size mismatch");
    return m;
}

std::vector<LegacyOrientedCornerTransferSnapshot>
expected_send_plan(
    const std::vector<LegacyOrientedCornerTransferSnapshot>
        &local_receive,
    int rank, int size)
{
    std::vector<int> local;
    local.reserve(local_receive.size() * kPackedMetaSize);
    for (const auto &record : local_receive)
    {
        const auto packed = pack_meta(record);
        local.insert(local.end(), packed.begin(), packed.end());
    }

    const int local_count = static_cast<int>(local.size());
    std::vector<int> counts(size, 0);
    MPI_Allgather(&local_count, 1, MPI_INT,
                  counts.data(), 1, MPI_INT,
                  MPI_COMM_WORLD);
    std::vector<int> displacements(size, 0);
    int total = 0;
    for (int peer = 0; peer < size; ++peer)
    {
        displacements[peer] = total;
        total += counts[peer];
    }
    std::vector<int> gathered(total);
    MPI_Allgatherv(
        local.data(), local_count, MPI_INT,
        gathered.data(), counts.data(),
        displacements.data(), MPI_INT, MPI_COMM_WORLD);

    std::vector<LegacyOrientedCornerTransferSnapshot> result;
    for (int origin = 0; origin < size; ++origin)
    {
        if (counts[origin] % kPackedMetaSize != 0)
            throw std::logic_error("gathered metadata is misaligned");
        const int end = displacements[origin] + counts[origin];
        for (int offset = displacements[origin];
             offset < end; offset += kPackedMetaSize)
        {
            auto record = unpack_meta(gathered.data() + offset);
            if (record.send_rank != rank)
                continue;
            record.role = LegacyOrientedCornerRole::Send;
            record.peer_rank = record.recv_rank;
            record.local_plan_order = result.size();
            Direction d1 = HALO_TOOLS::int_to_direction(
                record.dir1_send);
            Direction d2 = HALO_TOOLS::int_to_direction(
                record.dir2_send);
            record.send_box =
                HALO_BOX::make_2DCorner_innerghost_box(
                    edge_location(record.send_axis),
                    record.node_box_on_send,
                    d1, d2, record.nghost);
            result.push_back(record);
        }
    }
    return result;
}

Context context_for(
    int rank,
    const LegacyOrientedCornerTransferSnapshot &record,
    std::size_t order)
{
    const bool receive =
        record.role == LegacyOrientedCornerRole::Receive;
    return {rank,
            receive ? "receive" : "send",
            record.peer_rank,
            receive ? record.recv_block : record.send_block,
            receive ? record.recv_axis : record.send_axis,
            record.tag,
            order,
            record.peer_message_order};
}

void compare_legacy_record(
    Summary &summary, int rank,
    const LegacyOrientedCornerTransferSnapshot &expected,
    const LegacyOrientedCornerTransferSnapshot &actual,
    std::size_t order, std::uint64_t generation = 0)
{
    Context context = context_for(rank, expected, order);
    context.generation = generation;
    check_equal(summary, context, "role",
                static_cast<int>(expected.role),
                static_cast<int>(actual.role));
#define CHECK_MEMBER(name) \
    check_equal(summary, context, #name, expected.name, actual.name)
    CHECK_MEMBER(recv_rank);
    CHECK_MEMBER(send_rank);
    CHECK_MEMBER(recv_block);
    CHECK_MEMBER(send_block);
    CHECK_MEMBER(recv_axis);
    CHECK_MEMBER(send_axis);
    CHECK_MEMBER(sign);
    CHECK_MEMBER(dir1_send);
    CHECK_MEMBER(dir2_send);
    CHECK_MEMBER(peer_rank);
    CHECK_MEMBER(tag);
    CHECK_MEMBER(ncomp);
    CHECK_MEMBER(nghost);
    CHECK_MEMBER(message_length);
    CHECK_MEMBER(local_plan_order);
    CHECK_MEMBER(peer_message_order);
#undef CHECK_MEMBER
    compare_box(summary, context, "recv_box",
                expected.recv_box, actual.recv_box);
    compare_box(summary, context, "node_box_on_send",
                expected.node_box_on_send,
                actual.node_box_on_send);
    check_equal(summary, context, "send_box.has_value",
                expected.send_box.has_value(),
                actual.send_box.has_value());
    if (expected.send_box && actual.send_box)
        compare_box(summary, context, "send_box",
                    *expected.send_box, *actual.send_box);
    compare_transform(summary, context,
                      expected.trans_recv_to_send,
                      actual.trans_recv_to_send);
}

void compare_snapshot(
    Summary &summary, int rank,
    const LegacyOrientedCornerPlanSnapshot &expected,
    const LegacyOrientedCornerPlanSnapshot &actual,
    const char *phase, std::uint64_t generation)
{
    const Context context{
        rank, phase, -1, -1, -1, -1, 0, 0, generation};
    check_equal(summary, context, "group_name",
                expected.group_name, actual.group_name);
    for (int axis = 0; axis < 3; ++axis)
        check_equal(
            summary, context,
            "field_names[" + std::to_string(axis) + "]",
            expected.field_names[axis],
            actual.field_names[axis]);
    check_equal(summary, context, "ncomp",
                expected.ncomp, actual.ncomp);
    check_equal(summary, context, "nghost",
                expected.nghost, actual.nghost);
    check_equal(summary, context, "recv_count",
                expected.recv_plan.size(),
                actual.recv_plan.size());
    check_equal(summary, context, "send_count",
                expected.send_plan.size(),
                actual.send_plan.size());

    const std::size_t receive_count = std::min(
        expected.recv_plan.size(), actual.recv_plan.size());
    for (std::size_t index = 0; index < receive_count; ++index)
        compare_legacy_record(
            summary, rank, expected.recv_plan[index],
            actual.recv_plan[index], index, generation);

    const std::size_t send_count = std::min(
        expected.send_plan.size(), actual.send_plan.size());
    for (std::size_t index = 0; index < send_count; ++index)
        compare_legacy_record(
            summary, rank, expected.send_plan[index],
            actual.send_plan[index], index, generation);
}

Box3 legacy_box(const core::LogicalBox &box)
{
    return {{box.lo.i, box.lo.j, box.lo.k},
            {box.hi.i, box.hi.j, box.hi.k}};
}

TOPO::IndexTransform legacy_transform(
    const mesh::IndexTransform &transform)
{
    TOPO::IndexTransform result;
    for (int axis = 0; axis < 3; ++axis)
    {
        result.perm[axis] = transform.permutation[axis];
        result.sign[axis] = transform.sign[axis];
    }
    result.offset = {
        transform.offset.i,
        transform.offset.j,
        transform.offset.k};
    return result;
}

void compare_public_record(
    Summary &summary, int rank,
    const LegacyOrientedCornerTransferSnapshot &expected,
    const parallel::OrientedCornerTransfer &actual,
    std::size_t order,
    const std::array<std::string, 3> &fields)
{
    const Context context = context_for(rank, expected, order);
    check_equal(summary, context, "stage",
                static_cast<int>(parallel::HaloStage::Edge),
                static_cast<int>(actual.stage));
    check_equal(
        summary, context, "role",
        static_cast<int>(
            expected.role == LegacyOrientedCornerRole::Receive
                ? parallel::OrientedCornerRole::Receive
                : parallel::OrientedCornerRole::Send),
        static_cast<int>(actual.role));
    check_equal(summary, context, "source.rank",
                expected.send_rank, actual.source_block.rank);
    check_equal(summary, context, "source.block",
                expected.send_block,
                actual.source_block.local_block);
    check_equal(summary, context, "destination.rank",
                expected.recv_rank,
                actual.destination_block.rank);
    check_equal(summary, context, "destination.block",
                expected.recv_block,
                actual.destination_block.local_block);
    check_equal(summary, context, "send_axis",
                expected.send_axis, actual.send_axis);
    check_equal(summary, context, "recv_axis",
                expected.recv_axis, actual.recv_axis);
    if (expected.send_axis >= 0 && expected.send_axis < 3)
        check_equal(summary, context, "source_quantity",
                    fields[expected.send_axis],
                    actual.source_quantity);
    if (expected.recv_axis >= 0 && expected.recv_axis < 3)
        check_equal(summary, context, "destination_quantity",
                    fields[expected.recv_axis],
                    actual.destination_quantity);
    check_equal(summary, context, "orientation_sign",
                expected.sign, actual.orientation_sign);
    compare_box(summary, context, "destination_recv_box",
                expected.recv_box,
                legacy_box(actual.destination_recv_box));
    check_equal(summary, context, "source_send_box.has_value",
                expected.send_box.has_value(),
                actual.source_send_box.has_value());
    if (expected.send_box && actual.source_send_box)
        compare_box(summary, context, "source_send_box",
                    *expected.send_box,
                    legacy_box(*actual.source_send_box));
    compare_box(summary, context, "source_node_box",
                expected.node_box_on_send,
                legacy_box(actual.source_node_box));
    check_equal(summary, context, "source_direction[0]",
                expected.dir1_send,
                actual.source_directions[0]);
    check_equal(summary, context, "source_direction[1]",
                expected.dir2_send,
                actual.source_directions[1]);
    compare_transform(
        summary, context,
        expected.trans_recv_to_send,
        legacy_transform(actual.recv_to_send_transform));
    check_equal(summary, context, "peer_rank",
                expected.peer_rank, actual.peer_rank);
    check_equal(summary, context, "message_tag",
                expected.tag, actual.message_tag);
    check_equal(summary, context, "ncomp",
                expected.ncomp, actual.ncomp);
    check_equal(summary, context, "nghost",
                expected.nghost, actual.ghost_layers);
    check_equal(summary, context, "message_length",
                expected.message_length,
                actual.message_length);
    check_equal(summary, context, "local_plan_order",
                expected.local_plan_order,
                actual.local_plan_order);
    check_equal(summary, context, "peer_message_order",
                expected.peer_message_order,
                actual.peer_message_order);
    check_equal(
        summary, context, "packing_order",
        static_cast<int>(
            parallel::OrientedCornerPackingOrder::
                DestinationRecvIJKComponent),
        static_cast<int>(actual.packing_order));
}

using ValueKey =
    std::tuple<int, int, int, int, int, int>;

double marker(int rank, int block, int axis,
              int i, int j, int k, int component)
{
    const long long value =
        100000000000000LL * (rank + 1) +
        1000000000000LL * (block + 1) +
        10000000000LL * (axis + 1) +
        10000000LL * (i + 256) +
        10000LL * (j + 256) +
        10LL * (k + 256) + component;
    return static_cast<double>(value);
}

LegacyHaloPatternSnapshot inner_edge_pattern(
    const Halo &halo, StaggerLocation location, int nghost)
{
    for (const auto &pattern :
         halo.debug_component_copy_patterns(location, nghost))
        if (pattern.kind == LegacyHaloPatternKind::InnerCopy &&
            pattern.stage == HaloLevel::Corner2D)
            return pattern;
    throw std::invalid_argument(
        "missing legacy inner Edge pattern");
}

std::map<ValueKey, double> simulate_inner_edge(
    const Halo &halo, Field &field,
    const std::array<std::string, 3> &fields,
    int rank, int ncomp, int nghost)
{
    std::map<ValueKey, double> state;
    const auto read = [&](const ValueKey &key)
    {
        const auto found = state.find(key);
        if (found != state.end())
            return found->second;
        const auto &[axis, block, i, j, k, component] = key;
        return marker(rank, block, axis, i, j, k, component);
    };

    for (int recv_axis = 0; recv_axis < 3; ++recv_axis)
    {
        const auto pattern =
            inner_edge_pattern(
                halo, edge_location(recv_axis), nghost);
        for (const HaloRegion &region : pattern.regions)
        {
            const int send_axis =
                region.trans.perm[recv_axis];
            const int source_field_id =
                field.field_id(fields[send_axis]);
            const FieldBlock &source =
                field.field(source_field_id,
                            region.neighbor_block);
            const Int3 lo = source.get_lo();
            const Int3 hi = source.get_hi();
            for (int i = region.recv_box.lo.i;
                 i < region.recv_box.hi.i; ++i)
                for (int j = region.recv_box.lo.j;
                     j < region.recv_box.hi.j; ++j)
                    for (int k = region.recv_box.lo.k;
                         k < region.recv_box.hi.k; ++k)
                    {
                        int si = 0, sj = 0, sk = 0;
                        map_index(region.trans, i, j, k,
                                  si, sj, sk);
                        if (si < lo.i || si >= hi.i ||
                            sj < lo.j || sj >= hi.j ||
                            sk < lo.k || sk >= hi.k)
                            continue;
                        for (int component = 0;
                             component < ncomp; ++component)
                        {
                            const ValueKey source_key{
                                send_axis, region.neighbor_block,
                                si, sj, sk, component};
                            const ValueKey target_key{
                                recv_axis, region.this_block,
                                i, j, k, component};
                            state[target_key] =
                                region.trans.sign[recv_axis] *
                                read(source_key);
                        }
                    }
        }
    }
    return state;
}

void seed_fields(
    Field &field,
    const std::array<std::string, 3> &fields,
    int rank, int ncomp)
{
    for (int axis = 0; axis < 3; ++axis)
    {
        const int field_id = field.field_id(fields[axis]);
        for (int block = 0;
             block < field.num_blocks(); ++block)
        {
            FieldBlock &values =
                field.field(field_id, block);
            const Int3 lo = values.get_lo();
            const Int3 hi = values.get_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        for (int component = 0;
                             component < ncomp; ++component)
                            values(i, j, k, component) =
                                marker(rank, block, axis,
                                       i, j, k, component);
        }
    }
}

using PacketKey = std::tuple<int, int, int>;
using PacketMap = std::map<PacketKey, std::vector<double>>;

PacketMap gather_expected_packets(
    const std::vector<LegacyOrientedCornerTransferSnapshot>
        &send_plan,
    const std::map<ValueKey, double> &inner_state,
    const std::array<std::string, 3> &fields,
    Field &field, int rank, int size)
{
    const auto read = [&](const ValueKey &key)
    {
        const auto found = inner_state.find(key);
        if (found != inner_state.end())
            return found->second;
        const auto &[axis, block, i, j, k, component] = key;
        return marker(rank, block, axis, i, j, k, component);
    };

    std::vector<double> local;
    for (const auto &record : send_plan)
    {
        local.push_back(record.send_rank);
        local.push_back(record.recv_rank);
        local.push_back(record.tag);
        local.push_back(record.message_length);

        const int field_id =
            field.field_id(fields[record.send_axis]);
        const FieldBlock &source =
            field.field(field_id, record.send_block);
        const Int3 lo = source.get_lo();
        const Int3 hi = source.get_hi();
        for (int i = record.recv_box.lo.i;
             i < record.recv_box.hi.i; ++i)
            for (int j = record.recv_box.lo.j;
                 j < record.recv_box.hi.j; ++j)
                for (int k = record.recv_box.lo.k;
                     k < record.recv_box.hi.k; ++k)
                {
                    int si = 0, sj = 0, sk = 0;
                    map_index(
                        record.trans_recv_to_send,
                        i, j, k, si, sj, sk);
                    for (int component = 0;
                         component < record.ncomp; ++component)
                    {
                        double value = 0.0;
                        if (si >= lo.i && si < hi.i &&
                            sj >= lo.j && sj < hi.j &&
                            sk >= lo.k && sk < hi.k)
                        {
                            value = record.sign * read(ValueKey{
                                record.send_axis,
                                record.send_block,
                                si, sj, sk, component});
                        }
                        local.push_back(value);
                    }
                }
    }

    const int local_count = static_cast<int>(local.size());
    std::vector<int> counts(size, 0);
    MPI_Allgather(&local_count, 1, MPI_INT,
                  counts.data(), 1, MPI_INT,
                  MPI_COMM_WORLD);
    std::vector<int> displacements(size, 0);
    int total = 0;
    for (int peer = 0; peer < size; ++peer)
    {
        displacements[peer] = total;
        total += counts[peer];
    }
    std::vector<double> gathered(total);
    MPI_Allgatherv(
        local.data(), local_count, MPI_DOUBLE,
        gathered.data(), counts.data(),
        displacements.data(), MPI_DOUBLE,
        MPI_COMM_WORLD);

    PacketMap packets;
    for (int source_rank = 0;
         source_rank < size; ++source_rank)
    {
        int cursor = displacements[source_rank];
        const int end = cursor + counts[source_rank];
        while (cursor < end)
        {
            if (end - cursor < 4)
                throw std::logic_error(
                    "truncated gathered payload header");
            const int sender =
                static_cast<int>(gathered[cursor++]);
            const int receiver =
                static_cast<int>(gathered[cursor++]);
            const int tag =
                static_cast<int>(gathered[cursor++]);
            const int length =
                static_cast<int>(gathered[cursor++]);
            if (length < 0 || cursor + length > end)
                throw std::logic_error(
                    "truncated gathered payload");
            std::vector<double> payload(
                gathered.begin() + cursor,
                gathered.begin() + cursor + length);
            cursor += length;
            const PacketKey key{sender, receiver, tag};
            if (!packets.emplace(key, std::move(payload)).second)
                throw std::logic_error(
                    "duplicate gathered payload key");
        }
    }
    return packets;
}

long long check_final_values(
    Summary &summary, int rank,
    const std::vector<LegacyOrientedCornerTransferSnapshot>
        &receive_plan,
    const PacketMap &packets,
    const std::array<std::string, 3> &fields,
    Field &field)
{
    struct Expected
    {
        double value;
        std::size_t region;
        int peer;
        int tag;
    };
    std::map<ValueKey, Expected> final_values;
    for (std::size_t order = 0;
         order < receive_plan.size(); ++order)
    {
        const auto &record = receive_plan[order];
        const PacketKey packet_key{
            record.send_rank, record.recv_rank, record.tag};
        const auto packet = packets.find(packet_key);
        if (packet == packets.end() ||
            static_cast<int>(packet->second.size()) !=
                record.message_length)
            throw std::logic_error(
                "missing or mismatched expected payload");

        int offset = 0;
        for (int i = record.recv_box.lo.i;
             i < record.recv_box.hi.i; ++i)
            for (int j = record.recv_box.lo.j;
                 j < record.recv_box.hi.j; ++j)
                for (int k = record.recv_box.lo.k;
                     k < record.recv_box.hi.k; ++k)
                    for (int component = 0;
                         component < record.ncomp; ++component)
                    {
                        const ValueKey key{
                            record.recv_axis,
                            record.recv_block,
                            i, j, k, component};
                        final_values[key] = {
                            packet->second[offset++],
                            order, record.send_rank, record.tag};
                    }
    }

    for (const auto &entry : final_values)
    {
        const auto &[axis, block, i, j, k, component] =
            entry.first;
        const Expected &expected = entry.second;
        const int field_id = field.field_id(fields[axis]);
        const double actual =
            field.field(field_id, block)(
                i, j, k, component);
        const Context context{
            rank, "numeric-recv", expected.peer,
            block, axis, expected.tag, expected.region};
        check_equal(summary, context, "ghost_value",
                    expected.value, actual);
    }
    return static_cast<long long>(final_values.size());
}

void compare_plan(
    Summary &summary, int rank,
    const LegacyOrientedCornerPlanSnapshot &captured,
    const parallel::OrientedCornerTransferPlan &plan)
{
    const Context context{rank};
    check_equal(summary, context, "plan.rank",
                rank, plan.rank());
    check_equal(summary, context, "plan.group",
                captured.group_name, plan.sync_group());
    check_equal(summary, context, "plan.form_kind",
                static_cast<int>(
                    parallel::OrientedFormKind::Edge1Form),
                static_cast<int>(plan.form_kind()));
    check_equal(summary, context, "plan.stage",
                static_cast<int>(parallel::HaloStage::Edge),
                static_cast<int>(plan.stage()));
    check_equal(summary, context, "plan.ncomp",
                captured.ncomp, plan.ncomp());
    check_equal(summary, context, "plan.nghost",
                captured.nghost, plan.ghost_layers());
    for (int axis = 0; axis < 3; ++axis)
    {
        check_equal(summary, context,
                    "quantity.name[" + std::to_string(axis) + "]",
                    captured.field_names[axis],
                    plan.quantities()[axis].name);
        check_equal(summary, context,
                    "quantity.axis[" + std::to_string(axis) + "]",
                    axis, plan.quantities()[axis].axis);
    }

    check_equal(summary, context, "receive_count",
                captured.recv_plan.size(),
                plan.receive_transfers().size());
    check_equal(summary, context, "send_count",
                captured.send_plan.size(),
                plan.send_transfers().size());
    const std::size_t nrecv = std::min(
        captured.recv_plan.size(),
        plan.receive_transfers().size());
    for (std::size_t index = 0; index < nrecv; ++index)
        compare_public_record(
            summary, rank, captured.recv_plan[index],
            plan.receive_transfers()[index], index,
            captured.field_names);
    const std::size_t nsend = std::min(
        captured.send_plan.size(),
        plan.send_transfers().size());
    for (std::size_t index = 0; index < nsend; ++index)
        compare_public_record(
            summary, rank, captured.send_plan[index],
            plan.send_transfers()[index], index,
            captured.field_names);
}

template <class Operation>
bool rejected(Operation operation)
{
    try
    {
        operation();
    }
    catch (const std::invalid_argument &)
    {
        return true;
    }
    return false;
}

void run_contract_tests(Summary &summary, int rank, int size)
{
    using parallel::FieldLocation;
    using parallel::FormQuantityIdentity;
    using parallel::OrientedCornerRole;
    using parallel::OrientedCornerTransfer;
    using parallel::OrientedCornerTransferPlan;
    using parallel::OrientedFormKind;

    std::array<FormQuantityIdentity, 3> quantities = {{
        {"E_xi", OrientedFormKind::Edge1Form,
         FieldLocation::EdgeXi, 0},
        {"E_eta", OrientedFormKind::Edge1Form,
         FieldLocation::EdgeEta, 1},
        {"E_zeta", OrientedFormKind::Edge1Form,
         FieldLocation::EdgeZeta, 2}}};

    const int peer = (rank + 1) % size;
    OrientedCornerTransfer receive;
    receive.role = OrientedCornerRole::Receive;
    receive.source_block = {peer, 0};
    receive.destination_block = {rank, 0};
    receive.send_axis = 1;
    receive.recv_axis = 0;
    receive.source_quantity = "E_eta";
    receive.destination_quantity = "E_xi";
    receive.orientation_sign = -1;
    receive.destination_recv_box = {{0, 0, 0}, {1, 1, 1}};
    receive.source_node_box = {{0, 0, 0}, {1, 1, 1}};
    receive.source_directions = {1, 2};
    receive.recv_to_send_transform.permutation = {1, 2, 0};
    receive.recv_to_send_transform.sign = {-1, 1, 1};
    receive.peer_rank = peer;
    receive.message_tag = 0;
    receive.ncomp = 1;
    receive.ghost_layers = 1;
    receive.message_length = 1;

    OrientedCornerTransferPlan valid(
        rank, size, "synthetic-edge", quantities,
        1, 1, {receive}, {});
    const Context context{rank, "synthetic", peer, 0, 0, 0, 0};
    check_equal(summary, context, "direct_permutation_send_axis",
                1, valid.receive_transfers()[0].send_axis);
    check_equal(summary, context, "differs_from_face_inverse_axis",
                true,
                valid.receive_transfers()[0].send_axis != 2);
    check_equal(summary, context, "negative_sign",
                -1,
                valid.receive_transfers()[0].orientation_sign);

    auto invalid_axis = receive;
    invalid_axis.send_axis = 2;
    check_equal(summary, context, "invalid_axis_rejected",
                true, rejected([&]
                {
                    OrientedCornerTransferPlan(
                        rank, size, "bad-axis", quantities,
                        1, 1, {invalid_axis}, {});
                }));

    auto invalid_transform = receive;
    invalid_transform.recv_to_send_transform.permutation =
        {1, 1, 0};
    check_equal(summary, context, "invalid_transform_rejected",
                true, rejected([&]
                {
                    OrientedCornerTransferPlan(
                        rank, size, "bad-transform",
                        quantities, 1, 1,
                        {invalid_transform}, {});
                }));

    auto invalid_rank = receive;
    invalid_rank.source_block.rank = size;
    invalid_rank.peer_rank = size;
    check_equal(summary, context, "invalid_rank_rejected",
                true, rejected([&]
                {
                    OrientedCornerTransferPlan(
                        rank, size, "bad-rank", quantities,
                        1, 1, {invalid_rank}, {});
                }));

    auto missing_send_box = receive;
    missing_send_box.role = OrientedCornerRole::Send;
    missing_send_box.source_block = {rank, 0};
    missing_send_box.destination_block = {peer, 0};
    missing_send_box.peer_rank = peer;
    check_equal(summary, context, "missing_send_box_rejected",
                true, rejected([&]
                {
                    OrientedCornerTransferPlan(
                        rank, size, "missing-send-box",
                        quantities, 1, 1, {},
                        {missing_send_box});
                }));

    auto invalid_length = receive;
    invalid_length.message_length = 2;
    check_equal(summary, context, "invalid_length_rejected",
                true, rejected([&]
                {
                    OrientedCornerTransferPlan(
                        rank, size, "bad-length", quantities,
                        1, 1, {invalid_length}, {});
                }));

    auto invalid_tag = receive;
    invalid_tag.message_tag = 1;
    check_equal(summary, context, "invalid_tag_rejected",
                true, rejected([&]
                {
                    OrientedCornerTransferPlan(
                        rank, size, "bad-tag", quantities,
                        1, 1, {invalid_tag}, {});
                }));

    auto invalid_role = receive;
    invalid_role.role = OrientedCornerRole::Send;
    check_equal(summary, context, "invalid_receive_role_rejected",
                true, rejected([&]
                {
                    OrientedCornerTransferPlan(
                        rank, size, "bad-role", quantities,
                        1, 1, {invalid_role}, {});
                }));
}
} // namespace

namespace Z0_TEST
{
bool RunPreparedOrientedCornerPlanTests(
    const Grid &grid,
    const TOPO::Topology &topology,
    Field &field,
    Halo &halo,
    int rank)
{
    int size = 0;
    PARALLEL::mpi_size(&size);
    Summary summary;
    long long prepared_plans = 0;
    long long receive_records = 0;
    long long send_records = 0;
    long long matched_pairs = 0;
    long long value_checks = 0;
    long long runtime_syncs = 0;
    long long rebuilds = 0;
    long long registry_invalidations = 0;
    long long cache_reuses = 0;
    long long level_upgrades = 0;
    long long zero_local_parallel_patch_ranks = 0;
    long long multiple_peer_ranks = 0;
    long long multiple_same_peer_message_ranks = 0;

    run_contract_tests(summary, rank, size);

    const auto requests = field.halo_requests();
    const auto request = std::find_if(
        requests.begin(), requests.end(),
        [](const FieldHaloRequest &candidate)
        {
            const std::string group =
                candidate.sync_group.empty()
                    ? candidate.field_name
                    : candidate.sync_group;
            return candidate.do_halo &&
                   candidate.orientation_aware &&
                   candidate.value_kind ==
                       FieldValueKind::EdgeCovariant1Form &&
                   static_cast<int>(candidate.level) >=
                       static_cast<int>(HaloLevel::Corner2D) &&
                   group == "Eedge";
        });

    const Context initial_context{rank};
    check_equal(
        summary, initial_context, "prepared_absent_before_prepare",
        false,
        halo.debug_prepared_edge_1form_corner_plan(
            "Eedge").has_value());
    check_equal(
        summary, initial_context, "prepared_adapter_rejects_absent",
        true,
        rejected([&]
        {
            const auto partition =
                compat::LegacyMeshPartitionAdapter::build(
                    grid, topology, rank, size);
            FieldHaloRequest absent_request;
            absent_request.field_name = "E_xi";
            absent_request.sync_group = "Eedge";
            absent_request.location = StaggerLocation::EdgeXi;
            absent_request.value_kind =
                FieldValueKind::EdgeCovariant1Form;
            absent_request.ncomp = 1;
            absent_request.nghost = 1;
            absent_request.do_halo = true;
            absent_request.level = HaloLevel::Corner2D;
            absent_request.orientation_aware = true;
            compat::LegacyOrientedCornerPlanAdapter::build_prepared(
                partition, halo, absent_request);
        }));

    const std::uint64_t exchanges_before =
        halo.debug_prepared_edge_1form_corner_metadata_exchange_count();
    check_equal(
        summary, initial_context, "initial_prepare_succeeds", true,
        halo.prepare_edge_1form_corner_metadata());
    const std::uint64_t exchanges_after_first_prepare =
        halo.debug_prepared_edge_1form_corner_metadata_exchange_count();

    if (request == requests.end())
    {
        check_equal(
            summary, initial_context,
            "inactive_group_has_no_prepared_plan", false,
            halo.debug_prepared_edge_1form_corner_plan(
                "Eedge").has_value());
        check_equal(
            summary, initial_context,
            "inactive_group_has_no_metadata_exchange",
            exchanges_before, exchanges_after_first_prepare);
    }
    else
    {
        const std::array<std::string, 3> fields = {
            "E_xi", "E_eta", "E_zeta"};
        const auto expected_receive = expected_receive_plan(
            topology, request->ncomp, request->nghost);
        const auto expected_send = expected_send_plan(
            expected_receive, rank, size);

        check_equal(
            summary, initial_context,
            "first_prepare_metadata_exchange_count",
            exchanges_before + 1, exchanges_after_first_prepare);

        const auto prepared_first =
            halo.debug_prepared_edge_1form_corner_plan("Eedge");
        const auto generation_first =
            halo.debug_prepared_edge_1form_corner_plan_generation();
        check_equal(summary, initial_context, "prepared_present",
                    true, prepared_first.has_value());
        check_equal(summary, initial_context, "generation_present",
                    true, generation_first.has_value());

        const auto partition =
            compat::LegacyMeshPartitionAdapter::build(
                grid, topology, rank, size);
        if (prepared_first && generation_first)
        {
            ++prepared_plans;
            compare_snapshot(
                summary, rank,
                LegacyOrientedCornerPlanSnapshot{
                    "Eedge", fields, request->ncomp,
                    request->nghost, expected_receive,
                    expected_send},
                *prepared_first, "prepared", *generation_first);
            receive_records +=
                static_cast<long long>(
                    prepared_first->recv_plan.size());
            send_records +=
                static_cast<long long>(
                    prepared_first->send_plan.size());
            matched_pairs +=
                static_cast<long long>(
                    prepared_first->recv_plan.size());
            if (prepared_first->recv_plan.empty())
                ++zero_local_parallel_patch_ranks;
            std::map<int, int> messages_by_peer;
            for (const auto &transfer : prepared_first->recv_plan)
                ++messages_by_peer[transfer.peer_rank];
            if (messages_by_peer.size() >= 2)
                ++multiple_peer_ranks;
            if (std::any_of(
                    messages_by_peer.begin(), messages_by_peer.end(),
                    [](const auto &entry)
                    {
                        return entry.second >= 2;
                    }))
                ++multiple_same_peer_message_ranks;

            const auto public_prepared =
                compat::LegacyOrientedCornerPlanAdapter::build_prepared(
                    partition, halo, *request);
            compare_plan(
                summary, rank, *prepared_first, public_prepared);

            FieldHaloRequest wrong_ghost = *request;
            ++wrong_ghost.nghost;
            check_equal(
                summary, initial_context,
                "prepared_adapter_rejects_ghost_mismatch", true,
                rejected([&]
                {
                    compat::LegacyOrientedCornerPlanAdapter::
                        build_prepared(
                            partition, halo, wrong_ghost);
                }));
        }

        const auto inner_state =
            simulate_inner_edge(
                halo, field, fields, rank,
                request->ncomp, request->nghost);
        const PacketMap packets =
            gather_expected_packets(
                expected_send, inner_state, fields,
                field, rank, size);

        const auto capture_runtime =
            [&]()
            {
                halo.enable_edge_1form_corner_plan_capture(true);
                halo.sync_group("Eedge", HaloLevel::Corner2D);
                halo.enable_edge_1form_corner_plan_capture(false);
                ++runtime_syncs;
                return halo.debug_edge_1form_corner_plan(
                    "Eedge");
            };

        seed_fields(field, fields, rank, request->ncomp);
        const auto runtime_first = capture_runtime();
        check_equal(summary, initial_context,
                    "runtime_first_present", true,
                    runtime_first.has_value());
        if (prepared_first && generation_first && runtime_first)
        {
            compare_snapshot(
                summary, rank, *prepared_first, *runtime_first,
                "runtime-first", *generation_first);
            const auto runtime_public =
                compat::LegacyOrientedCornerPlanAdapter::build(
                    partition, halo, *request);
            compare_plan(
                summary, rank, *runtime_first, runtime_public);
            value_checks += check_final_values(
                summary, rank, expected_receive, packets,
                fields, field);
        }
        check_equal(
            summary, initial_context,
            "runtime_sync_does_not_prepare_metadata",
            exchanges_after_first_prepare,
            halo.debug_prepared_edge_1form_corner_metadata_exchange_count());

        check_equal(
            summary, initial_context,
            "repeated_prepare_uses_cache", true,
            halo.prepare_edge_1form_corner_metadata());
        const auto prepared_second =
            halo.debug_prepared_edge_1form_corner_plan("Eedge");
        const auto generation_second =
            halo.debug_prepared_edge_1form_corner_plan_generation();
        check_equal(
            summary, initial_context,
            "repeated_prepare_no_new_metadata_exchange",
            exchanges_after_first_prepare,
            halo.debug_prepared_edge_1form_corner_metadata_exchange_count());
        check_equal(summary, initial_context,
                    "repeated_prepare_generation_stable", true,
                    generation_first.has_value() &&
                        generation_second.has_value() &&
                        *generation_first == *generation_second);
        if (prepared_first && prepared_second && generation_first)
        {
            ++cache_reuses;
            compare_snapshot(
                summary, rank, *prepared_first, *prepared_second,
                "prepared-reuse", *generation_first);
        }

        seed_fields(field, fields, rank, request->ncomp);
        const auto runtime_second = capture_runtime();
        check_equal(summary, initial_context,
                    "runtime_second_present", true,
                    runtime_second.has_value());
        if (prepared_first && generation_first && runtime_second)
        {
            compare_snapshot(
                summary, rank, *prepared_first, *runtime_second,
                "runtime-second", *generation_first);
            value_checks += check_final_values(
                summary, rank, expected_receive, packets,
                fields, field);
        }

        halo.build_registered_patterns();
        ++rebuilds;
        check_equal(
            summary, initial_context,
            "rebuild_invalidates_prepared", false,
            halo.debug_prepared_edge_1form_corner_plan(
                "Eedge").has_value());
        check_equal(
            summary, initial_context,
            "rebuild_invalidates_generation", false,
            halo.debug_prepared_edge_1form_corner_plan_generation()
                .has_value());
        check_equal(
            summary, initial_context,
            "stale_prepared_adapter_rejected", true,
            rejected([&]
            {
                compat::LegacyOrientedCornerPlanAdapter::
                    build_prepared(partition, halo, *request);
            }));
        check_equal(
            summary, initial_context,
            "prepare_after_rebuild_succeeds", true,
            halo.prepare_edge_1form_corner_metadata());
        const std::uint64_t exchanges_after_rebuild =
            halo.debug_prepared_edge_1form_corner_metadata_exchange_count();
        check_equal(
            summary, initial_context,
            "rebuild_prepare_metadata_exchange_count",
            exchanges_after_first_prepare + 1,
            exchanges_after_rebuild);
        const auto prepared_rebuilt =
            halo.debug_prepared_edge_1form_corner_plan("Eedge");
        const auto generation_rebuilt =
            halo.debug_prepared_edge_1form_corner_plan_generation();
        check_equal(summary, initial_context,
                    "rebuilt_prepared_present", true,
                    prepared_rebuilt.has_value());
        check_equal(summary, initial_context,
                    "rebuilt_generation_changed", true,
                    generation_rebuilt.has_value() &&
                        generation_first.has_value() &&
                        *generation_rebuilt != *generation_first);

        seed_fields(field, fields, rank, request->ncomp);
        const auto runtime_rebuilt = capture_runtime();
        if (prepared_rebuilt && generation_rebuilt &&
            runtime_rebuilt)
        {
            compare_snapshot(
                summary, rank, *prepared_rebuilt, *runtime_rebuilt,
                "runtime-rebuild", *generation_rebuilt);
            value_checks += check_final_values(
                summary, rank, expected_receive, packets,
                fields, field);
        }

        std::vector<FieldHaloRequest> mutated_triplet;
        for (const FieldHaloRequest &candidate : requests)
        {
            const std::string group =
                candidate.sync_group.empty()
                    ? candidate.field_name
                    : candidate.sync_group;
            if (!candidate.do_halo ||
                !candidate.orientation_aware ||
                candidate.value_kind !=
                    FieldValueKind::EdgeCovariant1Form ||
                group != "Eedge")
                continue;

            FieldHaloRequest mutated = candidate;
            if (static_cast<int>(mutated.level) <
                static_cast<int>(HaloLevel::Corner3D))
            {
                mutated.level = HaloLevel::Corner3D;
                ++level_upgrades;
            }
            mutated_triplet.push_back(std::move(mutated));
        }
        check_equal(summary, initial_context,
                    "registry_mutation_has_complete_triplet",
                    std::size_t(3), mutated_triplet.size());
        halo.register_halo_fields(mutated_triplet);
        ++registry_invalidations;
        check_equal(
            summary, initial_context,
            "registry_mutation_invalidates_prepared", false,
            halo.debug_prepared_edge_1form_corner_plan(
                "Eedge").has_value());
        check_equal(
            summary, initial_context,
            "prepare_before_pattern_rebuild_rejected", false,
            halo.prepare_edge_1form_corner_metadata());
        check_equal(
            summary, initial_context,
            "unbuilt_prepare_has_no_metadata_exchange",
            exchanges_after_rebuild,
            halo.debug_prepared_edge_1form_corner_metadata_exchange_count());

        halo.build_registered_patterns();
        ++rebuilds;
        check_equal(
            summary, initial_context,
            "prepare_after_registry_rebuild_succeeds", true,
            halo.prepare_edge_1form_corner_metadata());
        const std::uint64_t exchanges_after_mutation =
            halo.debug_prepared_edge_1form_corner_metadata_exchange_count();
        check_equal(
            summary, initial_context,
            "registry_rebuild_metadata_exchange_count",
            exchanges_after_rebuild + 1,
            exchanges_after_mutation);
        const auto prepared_mutated =
            halo.debug_prepared_edge_1form_corner_plan("Eedge");
        const auto generation_mutated =
            halo.debug_prepared_edge_1form_corner_plan_generation();
        check_equal(summary, initial_context,
                    "mutated_prepared_present", true,
                    prepared_mutated.has_value());

        seed_fields(field, fields, rank, request->ncomp);
        const auto runtime_mutated = capture_runtime();
        if (prepared_mutated && generation_mutated &&
            runtime_mutated)
        {
            compare_snapshot(
                summary, rank, *prepared_mutated, *runtime_mutated,
                "runtime-mutation", *generation_mutated);
            value_checks += check_final_values(
                summary, rank, expected_receive, packets,
                fields, field);
        }
    }

    const long long global_prepared_plans =
        global_sum(prepared_plans);
    const long long global_receive_records =
        global_sum(receive_records);
    const long long global_send_records =
        global_sum(send_records);
    const long long global_pairs =
        global_sum(matched_pairs);
    const long long global_values =
        global_sum(value_checks);
    const long long global_runtime_syncs =
        global_sum(runtime_syncs);
    const long long global_rebuilds = global_sum(rebuilds);
    const long long global_registry_invalidations =
        global_sum(registry_invalidations);
    const long long global_cache_reuses =
        global_sum(cache_reuses);
    const long long global_level_upgrades =
        global_sum(level_upgrades);
    const long long global_zero_local_parallel_patch_ranks =
        global_sum(zero_local_parallel_patch_ranks);
    const long long global_multiple_peer_ranks =
        global_sum(multiple_peer_ranks);
    const long long global_multiple_same_peer_message_ranks =
        global_sum(multiple_same_peer_message_ranks);
    const long long global_checked =
        global_sum(summary.checked);
    const long long global_failed =
        global_sum(summary.failed);
    const long long global_metadata_exchanges =
        global_sum(static_cast<long long>(
            halo.debug_prepared_edge_1form_corner_metadata_exchange_count()));

    std::ostringstream detail;
    detail << "prepared_plans=" << global_prepared_plans
           << " recv_records=" << global_receive_records
           << " send_records=" << global_send_records
           << " matched_pairs=" << global_pairs
           << " checks=" << global_checked
           << " failed=" << global_failed
           << " value_checks=" << global_values
           << " runtime_syncs=" << global_runtime_syncs
           << " cache_reuses=" << global_cache_reuses
           << " rebuilds=" << global_rebuilds
           << " registry_invalidations="
           << global_registry_invalidations
           << " level_upgrades=" << global_level_upgrades
           << " zero_local_parallel_patch_ranks="
           << global_zero_local_parallel_patch_ranks
           << " multiple_peer_ranks="
           << global_multiple_peer_ranks
           << " multiple_same_peer_message_ranks="
           << global_multiple_same_peer_message_ranks
           << " prepare_metadata_exchange_groups="
           << global_metadata_exchanges
           << " synthetic_contract_cases=10";
    return print_result(
        "Prepared Edge-stage Edge1Form metadata parity",
        summary.failed == 0, detail.str());
}

bool RunPreparedEdge1FormEdgeExecutionTests(
    const Grid &grid,
    const TOPO::Topology &topology,
    Field &field,
    Halo &halo,
    int rank)
{
    int size = 0;
    PARALLEL::mpi_size(&size);
    Summary summary;
    const Context context{rank};
    const std::array<std::string, 3> fields = {
        "E_xi", "E_eta", "E_zeta"};

    const auto requests = field.halo_requests();
    const auto request = std::find_if(
        requests.begin(), requests.end(),
        [](const FieldHaloRequest &candidate)
        {
            const std::string group =
                candidate.sync_group.empty()
                    ? candidate.field_name
                    : candidate.sync_group;
            return candidate.do_halo &&
                   candidate.orientation_aware &&
                   candidate.value_kind ==
                       FieldValueKind::EdgeCovariant1Form &&
                   static_cast<int>(candidate.level) >=
                       static_cast<int>(HaloLevel::Corner2D) &&
                   group == "Eedge";
        });

    check_equal(
        summary, context, "explicit_legacy_deactivation", true,
        halo.deactivate_prepared_edge_1form_edge_execution());
    check_equal(
        summary, context, "legacy_mode_after_deactivation",
        static_cast<int>(Edge1FormEdgeExecutionMode::Legacy),
        static_cast<int>(halo.debug_edge_1form_edge_execution_mode()));
    check_equal(
        summary, context, "legacy_mode_is_not_active", false,
        halo.debug_prepared_edge_1form_edge_execution_active());

    const std::uint64_t prepare_before =
        halo.debug_prepared_edge_1form_corner_metadata_exchange_count();
    check_equal(
        summary, context, "prepare_before_activation", true,
        halo.prepare_edge_1form_corner_metadata());

    if (request == requests.end())
    {
        check_equal(
            summary, context, "inactive_orientation_has_no_snapshot",
            false,
            halo.debug_prepared_edge_1form_corner_plan(
                "Eedge").has_value());
        check_equal(
            summary, context, "inactive_orientation_activation_rejected",
            false,
            halo.activate_prepared_edge_1form_edge_execution());
        check_equal(
            summary, context, "inactive_orientation_stays_legacy",
            static_cast<int>(Edge1FormEdgeExecutionMode::Legacy),
            static_cast<int>(
                halo.debug_edge_1form_edge_execution_mode()));

        const long long global_checked = global_sum(summary.checked);
        const long long global_failed = global_sum(summary.failed);
        std::ostringstream detail;
        detail << "orientation_aware=0 checks=" << global_checked
               << " failed=" << global_failed
               << " prepare_metadata_exchange_groups="
               << (halo.debug_prepared_edge_1form_corner_metadata_exchange_count() -
                   prepare_before);
        return print_result(
            "Prepared Edge-stage Edge1Form execution",
            summary.failed == 0, detail.str());
    }

    const auto prepared =
        halo.debug_prepared_edge_1form_corner_plan("Eedge");
    check_equal(summary, context, "prepared_snapshot_present", true,
                prepared.has_value());
    if (!prepared)
    {
        const long long global_checked = global_sum(summary.checked);
        const long long global_failed = global_sum(summary.failed);
        std::ostringstream detail;
        detail << "prepared_snapshot_missing checks=" << global_checked
               << " failed=" << global_failed;
        return print_result(
            "Prepared Edge-stage Edge1Form execution", false,
            detail.str());
    }

    const auto expected_receive = expected_receive_plan(
        topology, request->ncomp, request->nghost);
    const auto expected_send = expected_send_plan(
        expected_receive, rank, size);
    const auto inner_state = simulate_inner_edge(
        halo, field, fields, rank,
        request->ncomp, request->nghost);
    const PacketMap packets = gather_expected_packets(
        expected_send, inner_state, fields, field, rank, size);

    const auto capture_receive_values =
        [&]()
        {
            std::map<ValueKey, double> values;
            for (const auto &record : expected_receive)
                for (int i = record.recv_box.lo.i;
                     i < record.recv_box.hi.i; ++i)
                    for (int j = record.recv_box.lo.j;
                         j < record.recv_box.hi.j; ++j)
                        for (int k = record.recv_box.lo.k;
                             k < record.recv_box.hi.k; ++k)
                            for (int component = 0;
                                 component < record.ncomp; ++component)
                            {
                                const ValueKey key{
                                    record.recv_axis,
                                    record.recv_block,
                                    i, j, k, component};
                                values[key] =
                                    field.field(
                                        field.field_id(
                                            fields[record.recv_axis]),
                                        record.recv_block)(
                                        i, j, k, component);
                            }
            return values;
        };

    const auto compare_values =
        [&](const std::string &role,
            const std::map<ValueKey, double> &expected,
            const std::map<ValueKey, double> &actual)
        {
            Context value_context = context;
            value_context.role = role.c_str();
            check_equal(
                summary, value_context, "value_count",
                expected.size(), actual.size());
            for (const auto &entry : expected)
            {
                const auto found = actual.find(entry.first);
                check_equal(
                    summary, value_context, "value_present",
                    true, found != actual.end());
                if (found != actual.end())
                    check_equal(
                        summary, value_context, "ghost_value",
                        entry.second, found->second);
            }
        };

    const auto run_edge_sync = [&]()
    {
        halo.sync_group("Eedge", HaloLevel::Corner2D);
    };

    seed_fields(field, fields, rank, request->ncomp);
    const Edge1FormEdgeExecutionCounters legacy_before =
        halo.debug_edge_1form_edge_execution_counters();
    run_edge_sync();
    const Edge1FormEdgeExecutionCounters legacy_after =
        halo.debug_edge_1form_edge_execution_counters();
    const auto legacy_values = capture_receive_values();
    const long long legacy_value_checks = check_final_values(
        summary, rank, expected_receive, packets, fields, field);
    check_equal(
        summary, context, "legacy_metadata_alltoall_once",
        legacy_before.legacy_metadata_alltoall + 1,
        legacy_after.legacy_metadata_alltoall);
    check_equal(
        summary, context, "legacy_metadata_alltoallv_once",
        legacy_before.legacy_metadata_alltoallv + 1,
        legacy_after.legacy_metadata_alltoallv);
    check_equal(
        summary, context, "legacy_execution_once",
        legacy_before.legacy_executions + 1,
        legacy_after.legacy_executions);

    if (size > 1)
    {
        const Edge1FormEdgeExecutionMode mismatched_mode =
            rank == 0
                ? Edge1FormEdgeExecutionMode::Prepared
                : Edge1FormEdgeExecutionMode::Legacy;
        check_equal(
            summary, context, "mixed_rank_activation_rejected", false,
            halo.set_edge_1form_edge_execution_mode(mismatched_mode));
        check_equal(
            summary, context, "mixed_rank_activation_returns_legacy",
            static_cast<int>(Edge1FormEdgeExecutionMode::Legacy),
            static_cast<int>(
                halo.debug_edge_1form_edge_execution_mode()));
    }

    check_equal(
        summary, context, "collective_prepared_activation", true,
        halo.activate_prepared_edge_1form_edge_execution());
    check_equal(
        summary, context, "prepared_mode_after_activation",
        static_cast<int>(Edge1FormEdgeExecutionMode::Prepared),
        static_cast<int>(halo.debug_edge_1form_edge_execution_mode()));
    check_equal(
        summary, context, "prepared_execution_is_active", true,
        halo.debug_prepared_edge_1form_edge_execution_active());

    const std::uint64_t prepared_metadata_before =
        halo.debug_prepared_edge_1form_corner_metadata_exchange_count();
    const Edge1FormEdgeExecutionCounters prepared_before =
        halo.debug_edge_1form_edge_execution_counters();
    seed_fields(field, fields, rank, request->ncomp);
    run_edge_sync();
    const auto prepared_first_values = capture_receive_values();
    const long long prepared_first_value_checks = check_final_values(
        summary, rank, expected_receive, packets, fields, field);
    seed_fields(field, fields, rank, request->ncomp);
    run_edge_sync();
    const auto prepared_second_values = capture_receive_values();
    const long long prepared_second_value_checks = check_final_values(
        summary, rank, expected_receive, packets, fields, field);
    const Edge1FormEdgeExecutionCounters prepared_after =
        halo.debug_edge_1form_edge_execution_counters();

    compare_values(
        "prepared-first", legacy_values, prepared_first_values);
    compare_values(
        "prepared-second", legacy_values, prepared_second_values);
    check_equal(
        summary, context, "prepared_execution_count_twice",
        prepared_before.prepared_executions + 2,
        prepared_after.prepared_executions);
    check_equal(
        summary, context, "prepared_runtime_metadata_alltoall_zero",
        prepared_before.legacy_metadata_alltoall,
        prepared_after.legacy_metadata_alltoall);
    check_equal(
        summary, context, "prepared_runtime_metadata_alltoallv_zero",
        prepared_before.legacy_metadata_alltoallv,
        prepared_after.legacy_metadata_alltoallv);
    check_equal(
        summary, context, "prepared_reuses_prepare_snapshot",
        prepared_metadata_before,
        halo.debug_prepared_edge_1form_corner_metadata_exchange_count());
    check_equal(
        summary, context, "prepared_reuses_active_view", true,
        halo.debug_prepared_edge_1form_edge_execution_active());
    check_equal(
        summary, context, "prepared_send_message_order_count",
        2 * (legacy_after.send_messages -
             legacy_before.send_messages),
        prepared_after.send_messages -
            prepared_before.send_messages);
    check_equal(
        summary, context, "prepared_recv_message_order_count",
        2 * (legacy_after.recv_messages -
             legacy_before.recv_messages),
        prepared_after.recv_messages -
            prepared_before.recv_messages);
    check_equal(
        summary, context, "prepared_sent_payload_count",
        2 * (legacy_after.sent_doubles -
             legacy_before.sent_doubles),
        prepared_after.sent_doubles -
            prepared_before.sent_doubles);
    check_equal(
        summary, context, "prepared_recv_payload_count",
        2 * (legacy_after.received_doubles -
             legacy_before.received_doubles),
        prepared_after.received_doubles -
            prepared_before.received_doubles);

    const auto fill_nan_fields = [&]()
    {
        const double nan = std::nan("");
        for (int axis = 0; axis < 3; ++axis)
        {
            const int field_id = field.field_id(fields[axis]);
            for (int block = 0; block < field.num_blocks(); ++block)
            {
                FieldBlock &values = field.field(field_id, block);
                const Int3 lo = values.get_lo();
                const Int3 hi = values.get_hi();
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                            for (int component = 0;
                                 component < request->ncomp;
                                 ++component)
                                values(i, j, k, component) = nan;
            }
        }
    };
    const auto check_zero_values =
        [&](const std::string &role,
            const std::map<ValueKey, double> &values)
        {
            Context value_context = context;
            value_context.role = role.c_str();
            for (const auto &entry : values)
                check_equal(
                    summary, value_context, "nonfinite_packed_as_zero",
                    0.0, entry.second);
        };

    check_equal(
        summary, context, "deactivate_before_nonfinite_legacy", true,
        halo.deactivate_prepared_edge_1form_edge_execution());
    fill_nan_fields();
    run_edge_sync();
    const auto legacy_nan_values = capture_receive_values();
    check_zero_values("legacy-nonfinite", legacy_nan_values);
    check_equal(
        summary, context, "prepare_before_nonfinite_prepared", true,
        halo.prepare_edge_1form_corner_metadata());
    check_equal(
        summary, context, "activate_before_nonfinite_prepared", true,
        halo.activate_prepared_edge_1form_edge_execution());
    fill_nan_fields();
    run_edge_sync();
    const auto prepared_nan_values = capture_receive_values();
    check_zero_values("prepared-nonfinite", prepared_nan_values);
    compare_values(
        "nonfinite-legacy-prepared",
        legacy_nan_values, prepared_nan_values);

    halo.build_registered_patterns();
    check_equal(
        summary, context, "rebuild_keeps_prepared_mode_stale",
        static_cast<int>(Edge1FormEdgeExecutionMode::Prepared),
        static_cast<int>(halo.debug_edge_1form_edge_execution_mode()));
    check_equal(
        summary, context, "rebuild_invalidates_active_view", false,
        halo.debug_prepared_edge_1form_edge_execution_active());
    check_equal(
        summary, context, "explicit_deactivate_after_rebuild", true,
        halo.deactivate_prepared_edge_1form_edge_execution());
    check_equal(
        summary, context, "prepare_after_rebuild", true,
        halo.prepare_edge_1form_corner_metadata());
    check_equal(
        summary, context, "activate_after_rebuild", true,
        halo.activate_prepared_edge_1form_edge_execution());

    std::vector<FieldHaloRequest> triplet;
    for (const FieldHaloRequest &candidate : requests)
    {
        const std::string group =
            candidate.sync_group.empty()
                ? candidate.field_name
                : candidate.sync_group;
        if (candidate.do_halo &&
            candidate.orientation_aware &&
            candidate.value_kind ==
                FieldValueKind::EdgeCovariant1Form &&
            group == "Eedge")
            triplet.push_back(candidate);
    }
    check_equal(summary, context, "reregister_triplet_size",
                std::size_t(3), triplet.size());
    halo.register_halo_fields(triplet);
    check_equal(
        summary, context, "reregister_keeps_prepared_mode_stale",
        static_cast<int>(Edge1FormEdgeExecutionMode::Prepared),
        static_cast<int>(halo.debug_edge_1form_edge_execution_mode()));
    check_equal(
        summary, context, "reregister_invalidates_active_view", false,
        halo.debug_prepared_edge_1form_edge_execution_active());
    check_equal(
        summary, context, "explicit_deactivate_after_reregister", true,
        halo.deactivate_prepared_edge_1form_edge_execution());
    halo.build_registered_patterns();
    check_equal(
        summary, context, "prepare_after_reregister", true,
        halo.prepare_edge_1form_corner_metadata());
    check_equal(
        summary, context, "activate_after_reregister", true,
        halo.activate_prepared_edge_1form_edge_execution());
    check_equal(
        summary, context, "final_deactivate", true,
        halo.deactivate_prepared_edge_1form_edge_execution());

    const long long global_checked = global_sum(summary.checked);
    const long long global_failed = global_sum(summary.failed);
    const long long global_legacy_value_checks =
        global_sum(legacy_value_checks);
    const long long global_prepared_value_checks =
        global_sum(
            prepared_first_value_checks +
            prepared_second_value_checks);
    const long long global_recv_records =
        global_sum(static_cast<long long>(expected_receive.size()));
    const long long global_send_records =
        global_sum(static_cast<long long>(expected_send.size()));
    std::ostringstream detail;
    detail << "checks=" << global_checked
           << " failed=" << global_failed
           << " recv_records=" << global_recv_records
           << " send_records=" << global_send_records
           << " legacy_value_checks=" << global_legacy_value_checks
           << " prepared_value_checks="
           << global_prepared_value_checks
           << " prepared_runtime_metadata_collectives=0"
           << " repeated_prepared_syncs=2";
    return print_result(
        "Prepared Edge-stage Edge1Form execution",
        summary.failed == 0, detail.str());
}

} // namespace Z0_TEST
