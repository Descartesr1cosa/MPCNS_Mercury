#include "Z0_Tests.h"
#include "Z0_TestCommon.h"

#include "0_basic/MPI_WRAPPER.h"
#include "1_grid/1_MPCNS_Grid.h"
#include "2_topology/Topology.h"
#include "3_field/Field.h"
#include "4_halo/Halo.h"
#include "mpcns/compat/legacy_mesh_partition_adapter.hpp"
#include "mpcns/compat/legacy_oriented_form_plan_adapter.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>

namespace
{
using namespace mpcns;

struct Summary
{
    long long checked = 0;
    long long failed = 0;
    long long reported = 0;
};

const char *kind_name(parallel::OrientedFormKind kind)
{
    return kind == parallel::OrientedFormKind::Edge1Form
               ? "Edge1Form"
               : "Face2Form";
}

template <class Expected, class Actual>
void check_equal(Summary &summary, int rank, const std::string &group,
                 const char *kind, const char *role, std::size_t region,
                 const std::string &field, const Expected &expected,
                 const Actual &actual)
{
    ++summary.checked;
    if (expected == actual)
        return;
    ++summary.failed;
    if (summary.reported++ >= 16)
        return;
    std::cerr << "[OrientedFormPlan][rank " << rank << "][stage Face]"
              << "[group " << group << "][kind " << kind << "]"
              << "[role " << role << "][region " << region << "]"
              << " field=" << field << std::setprecision(17)
              << " expected=" << expected
              << " actual=" << actual << "\n";
}

parallel::OrientedFormKind form_kind(LegacyOrientedFormKind kind)
{
    return kind == LegacyOrientedFormKind::Edge1Form
               ? parallel::OrientedFormKind::Edge1Form
               : parallel::OrientedFormKind::Face2Form;
}

parallel::FieldLocation location(StaggerLocation value)
{
    switch (value)
    {
    case StaggerLocation::EdgeXi: return parallel::FieldLocation::EdgeXi;
    case StaggerLocation::EdgeEt: return parallel::FieldLocation::EdgeEta;
    case StaggerLocation::EdgeZe: return parallel::FieldLocation::EdgeZeta;
    case StaggerLocation::FaceXi: return parallel::FieldLocation::FaceXi;
    case StaggerLocation::FaceEt: return parallel::FieldLocation::FaceEta;
    case StaggerLocation::FaceZe: return parallel::FieldLocation::FaceZeta;
    default: throw std::invalid_argument("not an oriented-form location");
    }
}

LegacyHaloPatternSnapshot face_pattern(
    const Halo &halo, StaggerLocation location_value, int ghost_layers,
    LegacyHaloPatternKind wanted)
{
    for (const auto &pattern :
         halo.debug_component_copy_patterns(location_value, ghost_layers))
        if (pattern.stage == HaloLevel::Corner1D &&
            pattern.kind == wanted)
            return pattern;
    throw std::invalid_argument("missing Face pattern in test oracle");
}

int legacy_source_axis(const TOPO::IndexTransform &transform,
                       int destination_axis)
{
    for (int axis = 0; axis < 3; ++axis)
        if (transform.perm[axis] == destination_axis)
            return axis;
    throw std::invalid_argument("invalid legacy permutation");
}

int legacy_sign(LegacyOrientedFormKind kind,
                const TOPO::IndexTransform &transform,
                int source_axis)
{
    if (kind == LegacyOrientedFormKind::Edge1Form)
        return transform.sign[source_axis];
    const int determinant =
        transform.sign[0] * transform.sign[1] * transform.sign[2];
    return determinant * transform.sign[source_axis] >= 0 ? 1 : -1;
}

void map_index(const TOPO::IndexTransform &transform,
               int i, int j, int k, int &io, int &jo, int &ko)
{
    const int source[3] = {i, j, k};
    int destination[3] = {0, 0, 0};
    const int offset[3] = {
        transform.offset.i, transform.offset.j, transform.offset.k};
    for (int axis = 0; axis < 3; ++axis)
        destination[transform.perm[axis]] =
            transform.sign[axis] * source[axis] + offset[axis];
    io = destination[0];
    jo = destination[1];
    ko = destination[2];
}

void inverse_map_index(const TOPO::IndexTransform &transform,
                       int i, int j, int k, int &io, int &jo, int &ko)
{
    const int destination[3] = {i, j, k};
    int source[3] = {0, 0, 0};
    const int offset[3] = {
        transform.offset.i, transform.offset.j, transform.offset.k};
    for (int axis = 0; axis < 3; ++axis)
        source[axis] = transform.sign[axis] *
            (destination[transform.perm[axis]] - offset[axis]);
    io = source[0];
    jo = source[1];
    ko = source[2];
}

double marker(int rank, int block, int axis,
              int i, int j, int k, int component)
{
    (void)i;
    (void)j;
    (void)k;
    return 100000000.0 * (rank + 1) +
           1000000.0 * (block + 1) +
           10000.0 * (axis + 1) + component;
}

void seed_triplet(Field &field,
                  const LegacyOrientedTripletSnapshot &group,
                  int rank)
{
    for (int axis = 0; axis < 3; ++axis)
    {
        const int field_id = field.field_id(group.field_names[axis]);
        for (int block = 0; block < field.num_blocks(); ++block)
        {
            FieldBlock &values = field.field(field_id, block);
            if (!values.is_allocated())
                continue;
            const Int3 lo = values.get_lo();
            const Int3 hi = values.get_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        for (int component = 0;
                             component < group.ncomp[axis]; ++component)
                            values(i, j, k, component) =
                                marker(rank, block, axis,
                                       i, j, k, component);
        }
    }
}

struct ValueCoverage
{
    long long total = 0;
    long long inner = 0;
    long long parallel = 0;
};

ValueCoverage check_transferred_values(
    Summary &summary, int rank, Field &field, const Halo &halo,
    const LegacyOrientedTripletSnapshot &group)
{
    using Target = std::tuple<int, int, int, int, int, int>;
    struct ExpectedValue
    {
        double value;
        std::size_t region;
        bool inner;
    };
    std::map<Target, double> expected_state;
    std::map<Target, ExpectedValue> expected_final;
    const char *form = kind_name(form_kind(group.form_kind));

    const auto expected_read = [&](const Target &key)
    {
        const auto it = expected_state.find(key);
        if (it != expected_state.end())
            return it->second;
        const auto &[axis, block, i, j, k, component] = key;
        return marker(rank, block, axis, i, j, k, component);
    };

    // Simulate the directed in-place Inner copies in the exact
    // legacy axis/region/ijk/component order.
    for (int destination_axis = 0; destination_axis < 3;
         ++destination_axis)
    {
        const auto destination_pattern =
            face_pattern(halo, group.locations[destination_axis],
                         group.nghost[0],
                        LegacyHaloPatternKind::InnerCopy);
        for (std::size_t region_index = 0;
             region_index < destination_pattern.regions.size();
             ++region_index)
        {
            const HaloRegion &destination_region =
                destination_pattern.regions[region_index];
            const int source_axis =
                legacy_source_axis(destination_region.trans,
                                 destination_axis);
            const auto source_pattern =
                face_pattern(halo, group.locations[source_axis],
                             group.nghost[0],
                             LegacyHaloPatternKind::InnerCopy);
            if (region_index >= source_pattern.regions.size())
                continue;

            const HaloRegion &source_region =
                source_pattern.regions[region_index];
            const int source_block = source_region.this_block;
            const int target_block =
                destination_region.neighbor_block;
            if (source_block < 0 || source_block >= field.num_blocks() ||
                target_block < 0 || target_block >= field.num_blocks())
                continue;

            const int source_field_id =
                field.field_id(group.field_names[source_axis]);
            const FieldBlock &source_values =
                field.field(source_field_id, source_block);
            const Int3 source_lo = source_values.get_lo();
            const Int3 source_hi = source_values.get_hi();
            const int sign =
                legacy_sign(group.form_kind, destination_region.trans,
                                source_axis);
            const Box3 &box = destination_region.recv_box;
            for (int i = box.lo.i; i < box.hi.i; ++i)
                for (int j = box.lo.j; j < box.hi.j; ++j)
                    for (int k = box.lo.k; k < box.hi.k; ++k)
                    {
                        int si = 0, sj = 0, sk = 0;
                        inverse_map_index(destination_region.trans,
                                        i, j, k, si, sj, sk);
                         if (si < source_lo.i || si >= source_hi.i ||
                           sj < source_lo.j || sj >= source_hi.j ||
                            sk < source_lo.k || sk >= source_hi.k)
                          continue;
                        for (int component = 0;
                               component < group.ncomp[0]; ++component)
                        {
                            const Target source_key{
                                source_axis, source_block, si, sj, sk, component};
                          const Target target_key{
                              destination_axis, target_block,
                                i, j, k, component};
                          const double value =
                              sign * expected_read(source_key);
                          expected_state[target_key] = value;
                          expected_final[target_key] =
                                {value, region_index, true};
                        }
                    }
        }
    }

    // Parallel Face unpacks after all Inner Face copies. The
    // peer value identity is dergived from the reverse directed
     // legacy pattern; later regions are the final writers.
    for (int destination_axis = 0; destination_axis < 3;
         ++destination_axis)
    {
        const auto pattern =
            face_pattern(halo, group.locations[destination_axis],
                         group.nghost[0],
                        LegacyHaloPatternKind::ParallelFace);
        for (std::size_t region_index = 0;
             region_index < pattern.regions.size();
             ++region_index)
        {
            const HaloRegion &region = pattern.regions[region_index];
            const int target_block = region.this_block;
            if (target_block < 0 ||
                target_block >= field.num_blocks())
                continue;

            const int source_axis =
                region.trans.perm[destination_axis];
            TOPO::IndexTransform inverse;
            for (int local_axis = 0; local_axis < 3;
                 ++local_axis)
          {
                const int peer_axis =
                    region.trans.perm[local_axis];
                inverse.perm[peer_axis] = local_axis;
                inverse.sign[peer_axis] =
                    region.trans.sign[local_axis];
            }
            const int sign =
                legacy_sign(group.form_kind, inverse, source_axis);
            const int source_rank = region.neighbor_rank;
            const int source_block = region.neighbor_block;
            const Box3 &box = region.recv_box;
            for (int i = box.lo.i; i < box.hi.i; ++i)
                for (int j = box.lo.j; j < box.hi.j; ++j)
                    for (int k = box.lo.k; k < box.hi.k; ++k)
                    {
                        int si = 0, sj = 0, sk = 0;
                        map_index(region.trans, i, j, k, si, sj, sk);
                        for (int component = 0;
                               component < group.ncomp[0]; ++component)
                        {
                            const double value =
                              sign * marker(source_rank, source_block,
                                        source_axis, si, sj, sk,
                                      component);
                            const Target target_key{
                                destination_axis, target_block,
                              i, j, k, component};
                          expected_state[target_key] = value;
                            expected_final[target_key] =
                                {value, region_index, false};
                        }
                    }
        }
    }

    ValueCoverage coverage;
    for (const auto &entry : expected_final)
    {
        const auto &[axis, block, i, j, k, component] = entry.first;
        const ExpectedValue &expected = entry.second;
        const int field_id = field.field_id(group.field_names[axis]);
        const FieldBlock &target = field.field(field_id, block);
        const double actual = target(i, j, k, component);
        check_equal(summary, rank, group.group_name, form,
                    expected.inner ? "inner-final-value"
                                  : "parallel-final-value",
                    expected.region, "transferred_value",
                    expected.value, actual);
        ++coverage.total;
        if (expected.inner)
            ++coverage.inner;
        else
            ++coverage.parallel;
    }
    return coverage;
}

void check_box(Summary &summary, int rank, const std::string &group,
               const char *kind, const char *role, std::size_t region,
               const char *prefix, const Box3 &expected,
               const core::LogicalBox &actual)
{
    const int e[6] = {expected.lo.i, expected.lo.j, expected.lo.k,
                      expected.hi.i, expected.hi.j, expected.hi.k};
    const int a[6] = {actual.lo.i, actual.lo.j, actual.lo.k,
                      actual.hi.i, actual.hi.j, actual.hi.k};
    const char *suffix[6] = {
        ".lo.i", ".lo.j", ".lo.k", ".hi.i", ".hi.j", ".hi.k"};
    for (int index = 0; index < 6; ++index)
        check_equal(summary, rank, group, kind, role, region,
                    std::string(prefix) + suffix[index],
                    e[index], a[index]);
}

void compare_transfer(
    Summary &summary, int rank,
    const LegacyOrientedTripletSnapshot &group,
    LegacyHaloPatternKind pattern_kind,
    int destination_axis,
    std::size_t region_index,
    std::size_t expected_order,
    const HaloRegion &source_region,
    const HaloRegion &destination_region,
    const parallel::OrientedFormTransfer &actual,
    long long &nonidentity,
    long long &negative,
    long long &inner,
    long long &cross_rank)
{
    const auto expected_kind = form_kind(group.form_kind);
    const char *form = kind_name(expected_kind);
    const bool is_parallel =
        pattern_kind == LegacyHaloPatternKind::ParallelFace;
    const char *role = is_parallel ? "parallel-send+recv" : "inner-copy";
    const int source_axis =
        legacy_source_axis(destination_region.trans, destination_axis);
    const int sign =
        legacy_sign(group.form_kind, destination_region.trans, source_axis);

    check_equal(summary, rank, group.group_name, form, role, region_index,
                "stage", (int)parallel::HaloStage::Face, (int)actual.stage);
    check_equal(summary, rank, group.group_name, form, role, region_index,
                "transfer_kind",
                (int)(is_parallel
                          ? parallel::OrientedTransferKind::Parallel
                          : parallel::OrientedTransferKind::Inner),
                (int)actual.kind);
    check_equal(summary, rank, group.group_name, form, role, region_index,
                "source.rank", source_region.this_rank,
                actual.source_block.rank);
    check_equal(summary, rank, group.group_name, form, role, region_index,
                "source.block", source_region.this_block,
                actual.source_block.local_block);
    check_equal(summary, rank, group.group_name, form, role, region_index,
                "destination.rank", destination_region.neighbor_rank,
                actual.destination_block.rank);
    check_equal(summary, rank, group.group_name, form, role, region_index,
                "destination.block", destination_region.neighbor_block,
                actual.destination_block.local_block);
    check_equal(summary, rank, group.group_name, form, role, region_index,
                "source_axis", source_axis, actual.source_axis);
    check_equal(summary, rank, group.group_name, form, role, region_index,
                "destination_axis", destination_axis,
                actual.destination_axis);
    check_equal(summary, rank, group.group_name, form, role, region_index,
                "source_quantity", group.field_names[source_axis],
                actual.source_quantity);
    check_equal(summary, rank, group.group_name, form, role, region_index,
                "destination_quantity", group.field_names[destination_axis],
                actual.destination_quantity);
    check_equal(summary, rank, group.group_name, form, role, region_index,
                "ncomp", group.ncomp[0], actual.ncomp);
    check_equal(summary, rank, group.group_name, form, role, region_index,
                "nghost", group.nghost[0], actual.ghost_layers);
    check_equal(summary, rank, group.group_name, form, role, region_index,
                "source_region_order", region_index,
                actual.source_region_order);
    check_equal(summary, rank, group.group_name, form, role, region_index,
                "destination_region_order", region_index,
                actual.destination_region_order);
    check_equal(summary, rank, group.group_name, form, role, region_index,
                "transfer_order", expected_order, actual.transfer_order);
    check_equal(summary, rank, group.group_name, form, role, region_index,
                "packing_order",
                (int)parallel::PackingOrder::TransformedDestinationIJKComponent,
                (int)actual.packing_order);

    check_box(summary, rank, group.group_name, form, role, region_index,
              "send_box", source_region.send_box, actual.send_box);
    check_box(summary, rank, group.group_name, form, role, region_index,
              "recv_box", destination_region.recv_box, actual.recv_box);
    for (int axis = 0; axis < 3; ++axis)
    {
        check_equal(summary, rank, group.group_name, form, role, region_index,
                    "transform.permutation[" + std::to_string(axis) + "]",
                    destination_region.trans.perm[axis],
                    actual.transform.permutation[axis]);
        check_equal(summary, rank, group.group_name, form, role, region_index,
                    "transform.sign[" + std::to_string(axis) + "]",
                    destination_region.trans.sign[axis],
                    actual.transform.sign[axis]);
        if (destination_region.trans.perm[axis] != axis)
            ++nonidentity;
    }
    check_equal(summary, rank, group.group_name, form, role, region_index,
                "transform.offset.i", destination_region.trans.offset.i,
                actual.transform.offset.i);
    check_equal(summary, rank, group.group_name, form, role, region_index,
                "transform.offset.j", destination_region.trans.offset.j,
                actual.transform.offset.j);
    check_equal(summary, rank, group.group_name, form, role, region_index,
                "transform.offset.k", destination_region.trans.offset.k,
                actual.transform.offset.k);
    check_equal(summary, rank, group.group_name, form, role, region_index,
                "orientation_sign", sign, actual.orientation_sign);

    const long long source_marker =
        100000LL * (source_axis + 1) +
        100LL * (region_index + 1) + 7;
    const long long actual_marker =
        (100000LL * (actual.source_axis + 1) +
         100LL * (actual.source_region_order + 1) + 7) *
        actual.orientation_sign;
    check_equal(summary, rank, group.group_name, form, role, region_index,
                "marked_transfer_value", source_marker * sign,
                actual_marker);

    check_equal(summary, rank, group.group_name, form, role, region_index,
                "send_tag.has_value", is_parallel,
                actual.send_tag.has_value());
    check_equal(summary, rank, group.group_name, form, role, region_index,
                "recv_tag.has_value", is_parallel,
                actual.recv_tag.has_value());
    if (actual.send_tag)
        check_equal(summary, rank, group.group_name, form, role, region_index,
                    "send_tag", destination_region.send_flag,
                    *actual.send_tag);
    if (actual.recv_tag)
        check_equal(summary, rank, group.group_name, form, role, region_index,
                    "recv_tag", destination_region.recv_flag,
                    *actual.recv_tag);

    if (sign < 0)
        ++negative;
    if (is_parallel)
        ++cross_rank;
    else
        ++inner;
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

void run_contract_tests(Summary &summary, int rank)
{
    using parallel::FieldLocation;
    using parallel::FormQuantityIdentity;
    using parallel::OrientedFormKind;
    using parallel::OrientedFormTransfer;
    using parallel::OrientedFormTransferPlan;

    std::array<FormQuantityIdentity, 3> quantities = {{
        {"q_xi", OrientedFormKind::Edge1Form, FieldLocation::EdgeXi, 0},
        {"q_eta", OrientedFormKind::Edge1Form, FieldLocation::EdgeEta, 1},
        {"q_zeta", OrientedFormKind::Edge1Form, FieldLocation::EdgeZeta, 2}}};

    OrientedFormTransfer transfer;
    transfer.source_block = {rank, 0};
    transfer.destination_block = {rank, 0};
    transfer.source_axis = 0;
    transfer.destination_axis = 0;
    transfer.source_quantity = "q_xi";
    transfer.destination_quantity = "q_xi";
    transfer.ncomp = 1;
    transfer.ghost_layers = 1;
    transfer.transform.permutation = {0, 1, 2};
    transfer.transform.sign = {1, 1, 1};
    transfer.orientation_sign = 1;

    auto missing_member = quantities;
    missing_member[1].name.clear();
    check_equal(summary, rank, "contract", "Edge1Form", "negative", 0,
                "missing_member_rejected", true,
                rejected([&]
                {
                    OrientedFormTransferPlan(
                        "q", OrientedFormKind::Edge1Form,
                        missing_member, 1, 1, {});
                }));

    auto mixed_kind = quantities;
    mixed_kind[2].form_kind = OrientedFormKind::Face2Form;
    check_equal(summary, rank, "contract", "Edge1Form", "negative", 1,
                "mixed_value_kind_rejected", true,
                rejected([&]
                {
                    OrientedFormTransferPlan(
                        "q", OrientedFormKind::Edge1Form,
                        mixed_kind, 1, 1, {});
                }));

    auto inconsistent_ghost = transfer;
    inconsistent_ghost.ghost_layers = 2;
    check_equal(summary, rank, "contract", "Edge1Form", "negative", 2,
                "inconsistent_ghost_rejected", true,
                rejected([&]
                {
                    OrientedFormTransferPlan(
                        "q", OrientedFormKind::Edge1Form,
                        quantities, 1, 1, {inconsistent_ghost});
                }));

    auto invalid_axis = transfer;
    invalid_axis.source_axis = 3;
    check_equal(summary, rank, "contract", "Edge1Form", "negative", 3,
                "invalid_axis_rejected", true,
                rejected([&]
                {
                    OrientedFormTransferPlan(
                        "q", OrientedFormKind::Edge1Form,
                        quantities, 1, 1, {invalid_axis});
                }));

    auto wrong_sign = transfer;
    wrong_sign.orientation_sign = -1;
    check_equal(summary, rank, "contract", "Edge1Form", "negative", 4,
                "wrong_orientation_sign_rejected", true,
                rejected([&]
                {
                    OrientedFormTransferPlan(
                        "q", OrientedFormKind::Edge1Form,
                        quantities, 1, 1, {wrong_sign});
                }));
}

void run_pure_mapping_tests(Summary &summary, int rank)
{
    mesh::IndexTransform transform;
    transform.permutation = {1, 2, 0};
    transform.sign = {-1, 1, 1};

    const int expected_source[3] = {2, 0, 1};
    const int expected_edge_sign[3] = {1, -1, 1};
    const int expected_face_sign[3] = {-1, 1, -1};
    for (int destination = 0; destination < 3; ++destination)
    {
        const int source =
            parallel::inverse_permuted_axis(transform, destination);
        check_equal(summary, rank, "synthetic", "mapping",
                    "nonidentity", destination, "source_axis",
                    expected_source[destination], source);
        check_equal(summary, rank, "synthetic", "Edge1Form",
                    "signed-marker", destination, "mapped_value",
                    expected_edge_sign[destination] *
                        (1000 + expected_source[destination]),
                    parallel::oriented_form_sign(
                        parallel::OrientedFormKind::Edge1Form,
                        transform, source) *
                        (1000 + source));
        check_equal(summary, rank, "synthetic", "Face2Form",
                    "signed-marker", destination, "mapped_value",
                    expected_face_sign[destination] *
                        (2000 + expected_source[destination]),
                    parallel::oriented_form_sign(
                        parallel::OrientedFormKind::Face2Form,
                        transform, source) *
                        (2000 + source));
    }
}
} // namespace

namespace Z0_TEST
{
bool RunOrientedFormPlanTests(const Grid &grid,
                              const TOPO::Topology &topology,
                              Field &field,
                              Halo &halo,
                              int rank)
{
    int size = 0;
    PARALLEL::mpi_size(&size);
    const auto partition =
        compat::LegacyMeshPartitionAdapter::build(
            grid, topology, rank, size);

    Summary summary;
    long long plans = 0;
    long long transfers = 0;
    long long nonidentity = 0;
    long long negative = 0;
    long long inner = 0;
    long long cross_rank = 0;

    const auto requests = field.halo_requests();
    auto groups = halo.debug_oriented_form_triplets();
    groups.erase(
        std::remove_if(
            groups.begin(), groups.end(),
            [&](const LegacyOrientedTripletSnapshot &group)
            {
                return std::none_of(
                    requests.begin(), requests.end(),
                    [&](const FieldHaloRequest &request)
                    {
                        const std::string request_group =
                            request.sync_group.empty()
                                ? request.field_name
                                : request.sync_group;
                        return request.do_halo &&
                               request.orientation_aware &&
                               request_group == group.group_name;
                    });
            }),
        groups.end());
    std::sort(groups.begin(), groups.end(),
              [](const auto &left, const auto &right)
              {
                  return left.group_name < right.group_name;
              });
    for (const auto &group : groups)
    {
        const auto request = std::find_if(
            requests.begin(), requests.end(),
            [&](const FieldHaloRequest &candidate)
            {
                const std::string candidate_group =
                    candidate.sync_group.empty()
                        ? candidate.field_name
                        : candidate.sync_group;
                return candidate.do_halo &&
                       candidate.orientation_aware &&
                       candidate_group == group.group_name;
            });
        if (request == requests.end())
            throw std::logic_error(
                "oriented group has no orientation-aware source request");
        const auto plan =
            compat::LegacyOrientedFormPlanAdapter::build(
                partition, halo, *request);
        ++plans;

        const auto expected_kind = form_kind(group.form_kind);
        const char *form = kind_name(expected_kind);
        check_equal(summary, rank, group.group_name, form, "plan", 0,
                    "group", group.group_name, plan.sync_group());
        check_equal(summary, rank, group.group_name, form, "plan", 0,
                    "form_kind", (int)expected_kind,
                    (int)plan.form_kind());
        check_equal(summary, rank, group.group_name, form, "plan", 0,
                    "ncomp", group.ncomp[0], plan.ncomp());
        check_equal(summary, rank, group.group_name, form, "plan", 0,
                    "nghost", group.nghost[0], plan.ghost_layers());
        for (int axis = 0; axis < 3; ++axis)
        {
            check_equal(summary, rank, group.group_name, form, "plan", axis,
                        "quantity.name", group.field_names[axis],
                        plan.quantities()[axis].name);
            check_equal(summary, rank, group.group_name, form, "plan", axis,
                        "quantity.form_kind", (int)expected_kind,
                        (int)plan.quantities()[axis].form_kind);
            check_equal(summary, rank, group.group_name, form, "plan", axis,
                        "quantity.location", (int)location(group.locations[axis]),
                        (int)plan.quantities()[axis].location);
            check_equal(summary, rank, group.group_name, form, "plan", axis,
                        "quantity.axis", axis,
                        plan.quantities()[axis].axis);
        }

        std::size_t order = 0;
        for (LegacyHaloPatternKind pattern_kind :
             {LegacyHaloPatternKind::InnerCopy,
              LegacyHaloPatternKind::ParallelFace})
        {
            for (int destination_axis = 0; destination_axis < 3;
                 ++destination_axis)
            {
                const auto destination_pattern =
                    face_pattern(halo, group.locations[destination_axis],
                                 group.nghost[0], pattern_kind);
                for (std::size_t region_index = 0;
                     region_index < destination_pattern.regions.size();
                     ++region_index, ++order)
                {
                    const HaloRegion &destination_region =
                        destination_pattern.regions[region_index];
                    const int source_axis =
                        legacy_source_axis(destination_region.trans,
                                           destination_axis);
                    const auto source_pattern =
                        face_pattern(halo, group.locations[source_axis],
                                     group.nghost[0], pattern_kind);
                    if (order >= plan.transfers().size() ||
                        region_index >= source_pattern.regions.size())
                    {
                        ++summary.checked;
                        ++summary.failed;
                        continue;
                    }
                    compare_transfer(
                        summary, rank, group, pattern_kind,
                        destination_axis, region_index, order,
                        source_pattern.regions[region_index],
                        destination_region, plan.transfers()[order],
                        nonidentity, negative, inner, cross_rank);
                    ++transfers;
                }
            }
        }
        check_equal(summary, rank, group.group_name, form, "plan", 0,
                    "transfer_count", order, plan.transfers().size());
    }

    long long value_checks = 0;
    long long inner_value_checks = 0;
    long long parallel_value_checks = 0;
    for (const auto &group : groups)
    {
        seed_triplet(field, group, rank);
        halo.sync_group(group.group_name, HaloLevel::Corner1D);
        const ValueCoverage coverage =
            check_transferred_values(summary, rank, field, halo, group);
        value_checks += coverage.total;
        inner_value_checks += coverage.inner;
        parallel_value_checks += coverage.parallel;
    }

    FieldHaloRequest component_request;
    component_request.field_name = "component-form";
    component_request.sync_group = "component-form-group";
    component_request.value_kind = FieldValueKind::EdgeCovariant1Form;
    component_request.do_halo = true;
    component_request.orientation_aware = false;
    check_equal(summary, rank, "component-form-group", "Edge1Form",
                "negative", 0, "non_oriented_request_rejected", true,
                rejected([&]
                {
                    compat::LegacyOrientedFormPlanAdapter::build(
                        partition, halo, component_request);
                }));

    FieldHaloRequest mixed_request = component_request;
    mixed_request.orientation_aware = true;
    mixed_request.value_kind = FieldValueKind::CartesianVector;
    check_equal(summary, rank, "component-form-group", "mixed",
                "negative", 1, "unsupported_value_kind_rejected", true,
                rejected([&]
                {
                    compat::LegacyOrientedFormPlanAdapter::build(
                        partition, halo, mixed_request);
                }));

    FieldHaloRequest missing_request = component_request;
    missing_request.field_name = "missing-oriented-field";
    missing_request.sync_group = "missing-oriented-group";
    missing_request.orientation_aware = true;
    check_equal(summary, rank, "missing", "unknown", "plan", 0,
                "missing_group_rejected", true,
                rejected([&]
                {
                    compat::LegacyOrientedFormPlanAdapter::build(
                        partition, halo, missing_request);
                }));

    run_pure_mapping_tests(summary, rank);
    run_contract_tests(summary, rank);

    const long long global_plans = global_sum(plans);
    const long long global_transfers = global_sum(transfers);
    const long long global_checked = global_sum(summary.checked);
    const long long global_failed = global_sum(summary.failed);
    const long long global_nonidentity = global_sum(nonidentity);
    const long long global_negative = global_sum(negative);
    const long long global_inner = global_sum(inner);
    const long long global_cross_rank = global_sum(cross_rank);
    const long long global_value_checks = global_sum(value_checks);
    const long long global_inner_value_checks =
        global_sum(inner_value_checks);
    const long long global_parallel_value_checks =
        global_sum(parallel_value_checks);

    std::ostringstream detail;
    detail << "plans=" << global_plans
           << " transfers=" << global_transfers
           << " checks=" << global_checked
           << " failed=" << global_failed
           << " inner=" << global_inner
           << " parallel=" << global_cross_rank
           << " nonidentity_axis_entries=" << global_nonidentity
           << " negative_signs=" << global_negative
           << " value_checks=" << global_value_checks
           << " inner_value_checks=" << global_inner_value_checks
           << " parallel_value_checks=" << global_parallel_value_checks
           << " synthetic_nontrivial_cases=18";
    return print_result("Face oriented-form plan parity",
                        summary.failed == 0, detail.str());
}
} // namespace Z0_TEST
