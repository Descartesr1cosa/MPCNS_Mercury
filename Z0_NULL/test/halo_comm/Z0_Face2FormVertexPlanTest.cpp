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
#include "mpcns/compat/legacy_face_2form_vertex_plan_adapter.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <limits>
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
    std::size_t region = 0;
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
    std::cerr << "[Face2FormVertexPlan][rank " << context.rank
              << "][group Bface][stage Vertex]"
              << "[role " << context.role << "]"
              << "[peer " << context.peer << "]"
              << "[block " << context.block << "]"
              << "[axis " << context.axis << "]"
              << "[tag " << context.tag << "]"
              << "[region " << context.region << "]"
              << " field=" << field << std::setprecision(17)
              << " expected=" << expected
              << " actual=" << actual << "\n";
}

StaggerLocation face_location(int axis)
{
    const StaggerLocation locations[3] = {
        StaggerLocation::FaceXi,
        StaggerLocation::FaceEt,
        StaggerLocation::FaceZe};
    return locations[axis];
}

int volume(const Box3 &box)
{
    return (box.hi.i - box.lo.i) *
           (box.hi.j - box.lo.j) *
           (box.hi.k - box.lo.k);
}

int inverse_axis(const TOPO::IndexTransform &transform,
                 int destination_axis)
{
    for (int source_axis = 0; source_axis < 3; ++source_axis)
        if (transform.perm[source_axis] == destination_axis)
            return source_axis;
    throw std::logic_error(
        "test oracle found invalid permutation");
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

LegacyFace2FormVertexTransferSnapshot make_expected_receive(
    const TOPO::VertexPatch &patch,
    int recv_axis, int ncomp, int nghost,
    int tag, std::size_t order)
{
    Direction this_d1 =
        HALO_TOOLS::int_to_direction(patch.dir1);
    Direction this_d2 =
        HALO_TOOLS::int_to_direction(patch.dir2);
    Direction this_d3 =
        HALO_TOOLS::int_to_direction(patch.dir3);
    int send_dir1 = HALO_TOOLS::map_dir_to_neighbor(
        patch.dir1, patch.trans, true);
    int send_dir2 = HALO_TOOLS::map_dir_to_neighbor(
        patch.dir2, patch.trans, false);
    int send_dir3 = HALO_TOOLS::map_dir_to_neighbor(
        patch.dir3, patch.trans, false);
    if (patch.trans.perm[std::abs(patch.dir1) - 1] + 1 !=
        std::abs(send_dir1))
    {
        if (patch.trans.perm[std::abs(patch.dir1) - 1] + 1 ==
            std::abs(send_dir2))
            std::swap(send_dir1, send_dir2);
        else if (patch.trans.perm[std::abs(patch.dir1) - 1] + 1 ==
                 std::abs(send_dir3))
            std::swap(send_dir1, send_dir3);
        else
            throw std::logic_error(
                "test oracle Face2Form Vertex direction mapping broken");
    }

    LegacyFace2FormVertexTransferSnapshot result;
    result.role = LegacyOrientedCornerRole::Receive;
    result.recv_rank = patch.this_rank;
    result.send_rank = patch.nb_rank;
    result.recv_block = patch.this_block;
    result.send_block = patch.nb_block;
    result.recv_axis = recv_axis;
    result.send_axis = inverse_axis(
        patch.trans, recv_axis);
    result.sign =
        HALO_TOOLS::face_2form_orientation_sign(
            patch.trans, result.send_axis);
    result.recv_box = HALO_BOX::make_3DCorner_ghost_box(
        face_location(recv_axis), patch.this_box_node,
        this_d1, this_d2, this_d3, nghost);
    result.send_box = std::nullopt;
    result.node_box_on_send = patch.nb_box_node;
    result.dir1_send = send_dir1;
    result.dir2_send = send_dir2;
    result.dir3_send = send_dir3;
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

std::vector<LegacyFace2FormVertexTransferSnapshot>
expected_receive_plan(const TOPO::Topology &topology,
                      int ncomp, int nghost)
{
    std::vector<LegacyFace2FormVertexTransferSnapshot> result;
    std::map<int, int> next_tag;
    for (const TOPO::VertexPatch &patch :
         TOPO_VIEW::vertex_patches(
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

constexpr int kPackedMetaSize = 37;

std::array<int, kPackedMetaSize> pack_meta(
    const LegacyFace2FormVertexTransferSnapshot &m)
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
    p[n++] = m.dir3_send;
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

LegacyFace2FormVertexTransferSnapshot unpack_meta(const int *p)
{
    LegacyFace2FormVertexTransferSnapshot m;
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
    m.dir3_send = p[n++];
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

std::vector<LegacyFace2FormVertexTransferSnapshot>
expected_send_plan(
    const std::vector<LegacyFace2FormVertexTransferSnapshot>
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

    std::vector<LegacyFace2FormVertexTransferSnapshot> result;
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
            Direction d3 = HALO_TOOLS::int_to_direction(
                record.dir3_send);
            record.send_box =
                HALO_BOX::make_3DCorner_innerghost_box(
                    face_location(record.send_axis),
                    record.node_box_on_send,
                    d1, d2, d3, record.nghost);
            result.push_back(record);
        }
    }
    return result;
}

Context context_for(
    int rank,
    const LegacyFace2FormVertexTransferSnapshot &record,
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
            order};
}

void compare_legacy_record(
    Summary &summary, int rank,
    const LegacyFace2FormVertexTransferSnapshot &expected,
    const LegacyFace2FormVertexTransferSnapshot &actual,
    std::size_t order)
{
    const Context context = context_for(rank, expected, order);
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
    CHECK_MEMBER(dir3_send);
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
    const LegacyFace2FormVertexTransferSnapshot &expected,
    const parallel::Face2FormVertexTransfer &actual,
    std::size_t order,
    const std::array<std::string, 3> &fields)
{
    const Context context = context_for(rank, expected, order);
    check_equal(summary, context, "stage",
                static_cast<int>(parallel::HaloStage::Vertex),
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
    check_equal(summary, context, "source_direction[2]",
                expected.dir3_send,
                actual.source_directions[2]);
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

LegacyHaloPatternSnapshot inner_face2_vertex_pattern(
    const Halo &halo, StaggerLocation location, int nghost)
{
    for (const auto &pattern :
         halo.debug_component_copy_patterns(location, nghost))
        if (pattern.kind == LegacyHaloPatternKind::InnerCopy &&
            pattern.stage == HaloLevel::Corner3D)
            return pattern;
    throw std::invalid_argument(
        "missing legacy inner Face2Form Vertex pattern");
}

std::map<ValueKey, double> simulate_inner_face2_vertex(
    const Halo &halo, Field &field,
    const std::array<std::string, 3> &fields,
    int rank, int ncomp, int nghost,
    const std::map<ValueKey, double> &initial_state)
{
    std::map<ValueKey, double> state = initial_state;
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
            inner_face2_vertex_pattern(
                halo, face_location(recv_axis), nghost);
        for (const HaloRegion &region : pattern.regions)
        {
            const int send_axis =
                inverse_axis(region.trans, recv_axis);
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
                                HALO_TOOLS::
                                    face_2form_orientation_sign(
                                        region.trans, send_axis) *
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

std::map<ValueKey, double> snapshot_fields(
    Field &field,
    const std::array<std::string, 3> &fields,
    int ncomp)
{
    std::map<ValueKey, double> result;
    for (int axis = 0; axis < 3; ++axis)
    {
        const int field_id = field.field_id(fields[axis]);
        for (int block = 0;
             block < field.num_blocks(); ++block)
        {
            const FieldBlock &values =
                field.field(field_id, block);
            const Int3 lo = values.get_lo();
            const Int3 hi = values.get_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        for (int component = 0;
                             component < ncomp; ++component)
                            result[{axis, block, i, j, k,
                                    component}] =
                                values(i, j, k, component);
        }
    }
    return result;
}

void fill_fields(
    Field &field,
    const std::array<std::string, 3> &fields,
    int ncomp, double value)
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
                            values(i, j, k, component) = value;
        }
    }
}

using PacketKey = std::tuple<int, int, int>;
using PacketMap = std::map<PacketKey, std::vector<double>>;

PacketMap gather_expected_packets(
    const std::vector<LegacyFace2FormVertexTransferSnapshot>
        &send_plan,
    const std::map<ValueKey, double> &inner_state,
    const std::array<std::string, 3> &fields,
    Field &field, int rank, int size,
    long long &out_of_bounds_values,
    long long &nonfinite_values)
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
                        if (si < lo.i || si >= hi.i ||
                            sj < lo.j || sj >= hi.j ||
                            sk < lo.k || sk >= hi.k)
                        {
                            ++out_of_bounds_values;
                        }
                        else
                        {
                            const double source = read(ValueKey{
                                record.send_axis,
                                record.send_block,
                                si, sj, sk, component});
                            if (std::isfinite(source))
                                value = record.sign * source;
                            else
                                ++nonfinite_values;
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
    const std::vector<LegacyFace2FormVertexTransferSnapshot>
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


long long check_zero_final_values(
    Summary &summary, int rank,
    const std::vector<LegacyFace2FormVertexTransferSnapshot>
        &receive_plan,
    const std::array<std::string, 3> &fields,
    Field &field)
{
    struct Writer
    {
        std::size_t region;
        int peer;
        int tag;
    };
    std::map<ValueKey, Writer> final_writers;
    for (std::size_t order = 0;
         order < receive_plan.size(); ++order)
    {
        const auto &record = receive_plan[order];
        for (int i = record.recv_box.lo.i;
             i < record.recv_box.hi.i; ++i)
            for (int j = record.recv_box.lo.j;
                 j < record.recv_box.hi.j; ++j)
                for (int k = record.recv_box.lo.k;
                     k < record.recv_box.hi.k; ++k)
                    for (int component = 0;
                         component < record.ncomp; ++component)
                        final_writers[{
                            record.recv_axis,
                            record.recv_block,
                            i, j, k, component}] = {
                                order, record.send_rank,
                                record.tag};
    }

    for (const auto &entry : final_writers)
    {
        const auto &[axis, block, i, j, k, component] =
            entry.first;
        const Writer &writer = entry.second;
        const double actual =
            field.field(field.field_id(fields[axis]), block)(
                i, j, k, component);
        const Context context{
            rank, "nonfinite-recv", writer.peer,
            block, axis, writer.tag, writer.region};
        check_equal(summary, context,
                    "nonfinite_source_becomes_zero",
                    0.0, actual);
    }
    return static_cast<long long>(final_writers.size());
}

void compare_plan(
    Summary &summary, int rank,
    const LegacyFace2FormVertexPlanSnapshot &captured,
    const parallel::Face2FormVertexTransferPlan &plan)
{
    const Context context{rank};
    check_equal(summary, context, "plan.rank",
                rank, plan.rank());
    check_equal(summary, context, "plan.group",
                captured.group_name, plan.sync_group());
    check_equal(summary, context, "plan.form_kind",
                static_cast<int>(
                    parallel::OrientedFormKind::Face2Form),
                static_cast<int>(plan.form_kind()));
    check_equal(summary, context, "plan.stage",
                static_cast<int>(parallel::HaloStage::Vertex),
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
    using parallel::Face2FormVertexTransfer;
    using parallel::Face2FormVertexTransferPlan;
    using parallel::OrientedFormKind;

    std::array<FormQuantityIdentity, 3> quantities = {{
        {"B_xi", OrientedFormKind::Face2Form,
         FieldLocation::FaceXi, 0},
        {"B_eta", OrientedFormKind::Face2Form,
         FieldLocation::FaceEta, 1},
        {"B_zeta", OrientedFormKind::Face2Form,
         FieldLocation::FaceZeta, 2}}};

    const int peer = (rank + 1) % size;
    Face2FormVertexTransfer receive;
    receive.role = OrientedCornerRole::Receive;
    receive.source_block = {peer, 0};
    receive.destination_block = {rank, 0};
    receive.send_axis = 2;
    receive.recv_axis = 0;
    receive.source_quantity = "B_zeta";
    receive.destination_quantity = "B_xi";
    receive.orientation_sign = -1;
    receive.destination_recv_box = {{0, 0, 0}, {1, 1, 1}};
    receive.source_node_box = {{0, 0, 0}, {1, 1, 1}};
    receive.source_directions = {1, 2, 3};
    receive.recv_to_send_transform.permutation = {1, 2, 0};
    receive.recv_to_send_transform.sign = {-1, 1, 1};
    receive.peer_rank = peer;
    receive.message_tag = 0;
    receive.ncomp = 1;
    receive.ghost_layers = 1;
    receive.message_length = 1;

    Face2FormVertexTransferPlan valid(
        rank, size, "synthetic-face2-vertex", quantities,
        1, 1, {receive}, {});
    const Context context{rank, "synthetic", peer, 0, 0, 0, 0};
    check_equal(summary, context, "inverse_permutation_send_axis",
                2, valid.receive_transfers()[0].send_axis);
    check_equal(summary, context, "differs_from_direct_permutation",
                true,
                valid.receive_transfers()[0].send_axis != 1);
    check_equal(summary, context, "negative_sign",
                -1,
                valid.receive_transfers()[0].orientation_sign);

    const TOPO::IndexTransform direction_transform{
        {1, 2, 0}, {-1, 1, 1}, {0, 0, 0}};
    int mapped_dir1 = HALO_TOOLS::map_dir_to_neighbor(
        1, direction_transform, true);
    int mapped_dir2 = HALO_TOOLS::map_dir_to_neighbor(
        2, direction_transform, false);
    int mapped_dir3 = HALO_TOOLS::map_dir_to_neighbor(
        3, direction_transform, false);
    if (std::abs(mapped_dir1) != 2)
    {
        if (std::abs(mapped_dir2) == 2)
            std::swap(mapped_dir1, mapped_dir2);
        else
            std::swap(mapped_dir1, mapped_dir3);
    }
    check_equal(summary, context, "mapped_normal_axis",
                2, std::abs(mapped_dir1));
    check_equal(summary, context, "mapped_directions_distinct",
                true,
                std::abs(mapped_dir1) != std::abs(mapped_dir2) &&
                    std::abs(mapped_dir1) != std::abs(mapped_dir3) &&
                    std::abs(mapped_dir2) != std::abs(mapped_dir3));
    check_equal(summary, context, "three_direction_reordered",
                true,
                std::array<int, 3>{
                    std::abs(mapped_dir1),
                    std::abs(mapped_dir2),
                    std::abs(mapped_dir3)} !=
                    std::array<int, 3>{1, 2, 3});

    auto invalid_axis = receive;
    invalid_axis.send_axis = 1;
    check_equal(summary, context, "invalid_axis_rejected",
                true, rejected([&]
                {
                    Face2FormVertexTransferPlan(
                        rank, size, "bad-axis", quantities,
                        1, 1, {invalid_axis}, {});
                }));

    auto invalid_transform = receive;
    invalid_transform.recv_to_send_transform.permutation =
        {1, 1, 0};
    check_equal(summary, context, "invalid_transform_rejected",
                true, rejected([&]
                {
                    Face2FormVertexTransferPlan(
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
                    Face2FormVertexTransferPlan(
                        rank, size, "bad-rank", quantities,
                        1, 1, {invalid_rank}, {});
                }));

    auto zero_length = receive;
    zero_length.destination_recv_box =
        {{0, 0, 0}, {0, 1, 1}};
    zero_length.message_length = 0;
    Face2FormVertexTransferPlan zero_length_plan(
        rank, size, "zero-length", quantities,
        1, 1, {zero_length}, {});
    check_equal(summary, context, "zero_length_accepted",
                0, zero_length_plan
                       .receive_transfers()[0]
                       .message_length);

    auto wrong_formula_sign = receive;
    wrong_formula_sign.orientation_sign = 1;
    check_equal(summary, context,
                "wrong_2form_sign_rejected",
                true, rejected([&]
                {
                    Face2FormVertexTransferPlan(
                        rank, size, "wrong-sign",
                        quantities, 1, 1,
                        {wrong_formula_sign}, {});
                }));

    auto invalid_sign = receive;
    invalid_sign.orientation_sign = 0;
    check_equal(summary, context, "invalid_sign_rejected",
                true, rejected([&]
                {
                    Face2FormVertexTransferPlan(
                        rank, size, "zero-sign",
                        quantities, 1, 1,
                        {invalid_sign}, {});
                }));

    auto duplicate_directions = receive;
    duplicate_directions.source_directions = {1, -1, 3};
    check_equal(summary, context,
                "duplicate_directions_rejected",
                true, rejected([&]
                {
                    Face2FormVertexTransferPlan(
                        rank, size, "duplicate-directions",
                        quantities, 1, 1,
                        {duplicate_directions}, {});
                }));

    auto missing_send_box = receive;
    missing_send_box.role = OrientedCornerRole::Send;
    missing_send_box.source_block = {rank, 0};
    missing_send_box.destination_block = {peer, 0};
    missing_send_box.peer_rank = peer;
    check_equal(summary, context, "missing_send_box_rejected",
                true, rejected([&]
                {
                    Face2FormVertexTransferPlan(
                        rank, size, "missing-send-box",
                        quantities, 1, 1, {},
                        {missing_send_box});
                }));
}
} // namespace

namespace Z0_TEST
{
bool RunFace2FormVertexPlanTests(
    const Grid &grid,
    const TOPO::Topology &topology,
    Field &field,
    Halo &halo,
    int rank)
{
    int size = 0;
    PARALLEL::mpi_size(&size);
    Summary summary;
    long long plans = 0;
    long long recv_records = 0;
    long long send_records = 0;
    long long pair_checks = 0;
    long long value_checks = 0;
    long long nonidentity = 0;
    long long negative = 0;
    long long zero_length_messages = 0;
    long long multi_block_ranks = 0;
    long long reordered_direction_records = 0;
    long long out_of_bounds_values = 0;
    long long oracle_nonfinite_values = 0;
    long long nonfinite_value_checks = 0;
    long long unoriented_path_checks = 0;

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
                       FieldValueKind::FaceContravariant2Form &&
                   static_cast<int>(candidate.level) >=
                       static_cast<int>(HaloLevel::Corner3D) &&
                   group == "Bface";
        });

    if (request != requests.end())
    {
        const Context initial_context{rank};
        check_equal(
            summary, initial_context,
            "snapshot.disabled_by_default", false,
            halo.debug_vertex_face_2form_plan(
                "Bface").has_value());

        const std::array<std::string, 3> fields = {
            "B_xi", "B_eta", "B_zeta"};
        const auto expected_recv =
            expected_receive_plan(
                topology, request->ncomp, request->nghost);
        const auto expected_send =
            expected_send_plan(
                expected_recv, rank, size);

        seed_fields(
            field, fields, rank, request->ncomp);
        halo.sync_group("Bface", HaloLevel::Corner1D);
        halo.sync_group("Bface", HaloLevel::Corner2D);
        const auto face_edge_state =
            snapshot_fields(
                field, fields, request->ncomp);
        const auto inner_state =
            simulate_inner_face2_vertex(
                halo, field, fields, rank,
                request->ncomp, request->nghost,
                face_edge_state);
        const PacketMap packets =
            gather_expected_packets(
                expected_send, inner_state, fields,
                field, rank, size,
                out_of_bounds_values,
                oracle_nonfinite_values);

        halo.enable_vertex_face_2form_plan_capture(true);
        halo.sync_group("Bface", HaloLevel::Corner3D);
        halo.enable_vertex_face_2form_plan_capture(false);
        check_equal(
            summary, initial_context,
            "vertex_face2_capture_does_not_create_edge_face2_snapshot",
            false,
            halo.debug_edge_face_2form_plan(
                "Bface").has_value());
        check_equal(
            summary, initial_context,
            "vertex_face2_capture_does_not_create_edge1_snapshot",
            false,
            halo.debug_edge_1form_corner_plan(
                "Bface").has_value());
        check_equal(
            summary, initial_context,
            "vertex_face2_capture_does_not_create_vertex1_snapshot",
            false,
            halo.debug_vertex_edge_1form_plan(
                "Bface").has_value());

        const auto captured =
            halo.debug_vertex_face_2form_plan("Bface");
        const Context plan_context{rank};
        check_equal(summary, plan_context,
                    "snapshot.present", true,
                    captured.has_value());
        if (captured)
        {
            ++plans;
            check_equal(summary, plan_context,
                        "snapshot.group", std::string("Bface"),
                        captured->group_name);
            check_equal(summary, plan_context,
                        "snapshot.fields.xi", fields[0],
                        captured->field_names[0]);
            check_equal(summary, plan_context,
                        "snapshot.fields.eta", fields[1],
                        captured->field_names[1]);
            check_equal(summary, plan_context,
                        "snapshot.fields.zeta", fields[2],
                        captured->field_names[2]);
            check_equal(summary, plan_context,
                        "snapshot.ncomp", request->ncomp,
                        captured->ncomp);
            check_equal(summary, plan_context,
                        "snapshot.nghost", request->nghost,
                        captured->nghost);
            check_equal(summary, plan_context,
                        "snapshot.recv_count",
                        expected_recv.size(),
                        captured->recv_plan.size());
            check_equal(summary, plan_context,
                        "snapshot.send_count",
                        expected_send.size(),
                        captured->send_plan.size());

            const std::size_t nrecv = std::min(
                expected_recv.size(),
                captured->recv_plan.size());
            for (std::size_t index = 0;
                 index < nrecv; ++index)
            {
                compare_legacy_record(
                    summary, rank, expected_recv[index],
                    captured->recv_plan[index], index);
                ++recv_records;
                ++pair_checks;
                for (int axis = 0; axis < 3; ++axis)
                    if (expected_recv[index]
                            .trans_recv_to_send.perm[axis] != axis)
                        ++nonidentity;
                if (expected_recv[index].sign < 0)
                    ++negative;
                if (std::array<int, 3>{
                        std::abs(expected_recv[index].dir1_send),
                        std::abs(expected_recv[index].dir2_send),
                        std::abs(expected_recv[index].dir3_send)} !=
                    std::array<int, 3>{1, 2, 3})
                    ++reordered_direction_records;
                if (expected_recv[index].message_length == 0)
                    ++zero_length_messages;
            }
            const std::size_t nsend = std::min(
                expected_send.size(),
                captured->send_plan.size());
            for (std::size_t index = 0;
                 index < nsend; ++index)
            {
                compare_legacy_record(
                    summary, rank, expected_send[index],
                    captured->send_plan[index], index);
                ++send_records;
            }

            const auto partition =
                compat::LegacyMeshPartitionAdapter::build(
                    grid, topology, rank, size);
            if (partition.local_blocks().size() >= 2)
                ++multi_block_ranks;
            const auto plan =
                compat::LegacyFace2FormVertexPlanAdapter::build(
                    partition, halo, *request);
            compare_plan(summary, rank, *captured, plan);

            value_checks =
                check_final_values(
                    summary, rank, expected_recv,
                    packets, fields, field);

            FieldHaloRequest lower_level = *request;
            lower_level.level = HaloLevel::Corner2D;
            check_equal(
                summary, plan_context,
                "non_vertex_request_rejected", true,
                rejected([&]
                {
                    compat::LegacyFace2FormVertexPlanAdapter::build(
                        partition, halo, lower_level);
                }));

            FieldHaloRequest component_request = *request;
            component_request.orientation_aware = false;
            check_equal(
                summary, plan_context,
                "unoriented_request_rejected", true,
                rejected([&]
                {
                    compat::LegacyFace2FormVertexPlanAdapter::build(
                        partition, halo, component_request);
                }));

            FieldHaloRequest missing = *request;
            missing.field_name = "missing_edge_form";
            missing.sync_group = "missing_edge_group";
            check_equal(
                summary, plan_context,
                "missing_metadata_rejected", true,
                rejected([&]
                {
                    compat::LegacyFace2FormVertexPlanAdapter::build(
                        partition, halo, missing);
                }));

        }

        seed_fields(
            field, fields, rank, request->ncomp);
        halo.sync_group("Bface", HaloLevel::Corner1D);
        halo.sync_group("Bface", HaloLevel::Corner2D);
        fill_fields(
            field, fields, request->ncomp,
            std::numeric_limits<double>::quiet_NaN());
        halo.sync_group("Bface", HaloLevel::Corner3D);
        nonfinite_value_checks =
            check_zero_final_values(
                summary, rank, expected_recv,
                fields, field);
        seed_fields(
            field, fields, rank, request->ncomp);

        halo.build_registered_patterns();
        check_equal(
            summary, plan_context,
            "snapshot_invalidated_on_rebuild", false,
            halo.debug_vertex_face_2form_plan(
                "Bface").has_value());

        halo.enable_vertex_face_2form_plan_capture(true);
        const std::vector<std::string> legacy_fields{
            fields[0], fields[1], fields[2]};
        halo.data_trans_face_axis_triplet(
            legacy_fields, HaloLevel::Corner3D);
        halo.enable_vertex_face_2form_plan_capture(false);
        check_equal(
            summary, plan_context,
            "unoriented_legacy_path_has_no_snapshot",
            false,
            halo.debug_vertex_face_2form_plan(
                "Bface").has_value());
        ++unoriented_path_checks;
        seed_fields(
            field, fields, rank, request->ncomp);
    }

    const long long global_plans = global_sum(plans);
    const long long global_recv = global_sum(recv_records);
    const long long global_send = global_sum(send_records);
    const long long global_pairs = global_sum(pair_checks);
    const long long global_values = global_sum(value_checks);
    const long long global_checked = global_sum(summary.checked);
    const long long global_failed = global_sum(summary.failed);
    const long long global_nonidentity = global_sum(nonidentity);
    const long long global_negative = global_sum(negative);
    const long long global_zero_length =
        global_sum(zero_length_messages);
    const long long global_multi_block =
        global_sum(multi_block_ranks);
    const long long global_reordered_directions =
        global_sum(reordered_direction_records);
    const long long global_out_of_bounds =
        global_sum(out_of_bounds_values);
    const long long global_oracle_nonfinite =
        global_sum(oracle_nonfinite_values);
    const long long global_nonfinite_checks =
        global_sum(nonfinite_value_checks);
    const long long global_unoriented_checks =
        global_sum(unoriented_path_checks);

    std::ostringstream detail;
    detail << "plans=" << global_plans
           << " recv_records=" << global_recv
           << " send_records=" << global_send
           << " matched_pairs=" << global_pairs
           << " checks=" << global_checked
           << " failed=" << global_failed
           << " value_checks=" << global_values
           << " nonidentity_axis_entries="
           << global_nonidentity
           << " negative_signs=" << global_negative
           << " zero_length_messages="
           << global_zero_length
           << " multi_local_block_ranks="
           << global_multi_block
           << " reordered_direction_records="
           << global_reordered_directions
           << " out_of_bounds_values="
           << global_out_of_bounds
           << " oracle_nonfinite_values="
           << global_oracle_nonfinite
           << " nonfinite_value_checks="
           << global_nonfinite_checks
           << " unoriented_path_checks="
           << global_unoriented_checks
           << " synthetic_contract_cases=14";
    return print_result(
        "Vertex-stage Face2Form metadata snapshot parity",
        summary.failed == 0, detail.str());
}
} // namespace Z0_TEST
