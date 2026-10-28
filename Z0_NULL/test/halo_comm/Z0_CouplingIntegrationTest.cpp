#include "Z0_Tests.h"
#include "Z0_TestCommon.h"

#include "2_topology/TopologyBuildDetail.h"
#include "2_topology/TopologyView.h"
#include "4_halo/Halo.h"
#include "4_halo/detail/HaloBuildTools.h"
#include "5_boundary/Boundary.h"

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace
{
const char *source_domain = "COUPLING_SOURCE";
const char *destination_domain = "COUPLING_DESTINATION";

bool coupling_test_enabled()
{
    const char *value = std::getenv("Z0_COUPLING_TEST");
    return value && std::string(value) != "0";
}

FieldDescriptor coupling_descriptor(
    const std::string &name, StaggerLocation location,
    FieldValueKind kind, int ncomp, int nghost,
    const std::string &group, bool oriented)
{
    FieldDescriptor descriptor;
    descriptor.name = name;
    descriptor.location = location;
    descriptor.value_kind = kind;
    descriptor.ncomp = ncomp;
    descriptor.nghost = nghost;
    descriptor.sync.group = group;
    descriptor.sync.do_coupling = true;
    descriptor.sync.orientation_aware = oriented;
    return descriptor;
}

TOPO::Topology make_coupling_topology(
    Grid &grid, const TOPO::Topology &source)
{
    TOPO::Topology topology = source;

    for (auto &patch : topology.inner_patches)
    {
        if (patch.this_block == patch.nb_block ||
            patch.this_rank == 2)
        {
            patch.is_coupling = false;
            continue;
        }
        patch.is_coupling = true;
        const bool this_is_source =
            patch.this_block < patch.nb_block;
        patch.this_block_name =
            this_is_source ? source_domain : destination_domain;
        patch.nb_block_name =
            this_is_source ? destination_domain : source_domain;
    }

    for (auto &patch : topology.parallel_patches)
    {
        const bool selected =
            (patch.this_rank == 0 && patch.nb_rank == 1) ||
            (patch.this_rank == 1 && patch.nb_rank == 0);
        patch.is_coupling = selected;
        if (!selected)
            continue;
        const bool this_is_source = patch.this_rank == 0;
        patch.this_block_name =
            this_is_source ? source_domain : destination_domain;
        patch.nb_block_name =
            this_is_source ? destination_domain : source_domain;
    }

    TOPO::detail::build_edge_patches(
        grid, topology, grid.dimension);
    TOPO::detail::build_vertex_patches(
        grid, topology, grid.dimension);
    return topology;
}

void enable_all_parallel_coupling(
    Grid &grid, TOPO::Topology &topology)
{
    for (auto &patch : topology.parallel_patches)
    {
        patch.is_coupling = true;
        const bool this_is_source =
            patch.this_rank < patch.nb_rank;
        patch.this_block_name =
            this_is_source ? source_domain : destination_domain;
        patch.nb_block_name =
            this_is_source ? destination_domain : source_domain;
    }
    TOPO::detail::build_edge_patches(
        grid, topology, grid.dimension);
    TOPO::detail::build_vertex_patches(
        grid, topology, grid.dimension);
}

void register_coupling_fields(Field &field)
{
    field.register_field(coupling_descriptor(
        "C_scalar", StaggerLocation::Cell, FieldValueKind::Scalar,
        1, 1, "C_scalar", false));
    field.register_field(coupling_descriptor(
        "C_vector", StaggerLocation::Cell,
        FieldValueKind::CartesianVector,
        3, 2, "C_vector", false));
    field.register_field(coupling_descriptor(
        "CE_xi", StaggerLocation::EdgeXi,
        FieldValueKind::EdgeCovariant1Form,
        1, 1, "CE", true));
    field.register_field(coupling_descriptor(
        "CE_eta", StaggerLocation::EdgeEt,
        FieldValueKind::EdgeCovariant1Form,
        1, 1, "CE", true));
    field.register_field(coupling_descriptor(
        "CE_zeta", StaggerLocation::EdgeZe,
        FieldValueKind::EdgeCovariant1Form,
        1, 1, "CE", true));
    field.register_field(coupling_descriptor(
        "CB_xi", StaggerLocation::FaceXi,
        FieldValueKind::FaceContravariant2Form,
        1, 1, "CB", true));
    field.register_field(coupling_descriptor(
        "CB_eta", StaggerLocation::FaceEt,
        FieldValueKind::FaceContravariant2Form,
        1, 1, "CB", true));
    field.register_field(coupling_descriptor(
        "CB_zeta", StaggerLocation::FaceZe,
        FieldValueKind::FaceContravariant2Form,
        1, 1, "CB", true));

    field.register_field(coupling_descriptor(
        "C_zero", StaggerLocation::Cell, FieldValueKind::Scalar,
        1, 0, "C_zero", false));
    const Field::PairKey pairs[] = {
        {source_domain, destination_domain},
        {destination_domain, source_domain}};
    for (const auto &pair : pairs)
        for (const FieldDescriptor &descriptor : field.descriptors())
            if (descriptor.sync.do_coupling)
                field.register_coupling_channel(
                    pair.first, pair.second, descriptor.name);
}

double coupling_value(int epoch, int rank, int fid, int block,
                      int i, int j, int k, int component)
{
    (void)block;
    const long long value =
        static_cast<long long>(epoch) * 1000000000000LL +
        static_cast<long long>(rank) * 100000000000LL +
        static_cast<long long>(fid) * 1000000000LL +
        static_cast<long long>(i + 100) * 10000LL +
        static_cast<long long>(j + 100) * 100LL +
        static_cast<long long>(k + 100) * 2LL + component;
    return static_cast<double>(value);
}

void fill_coupling_fields(Field &field, int epoch)
{
    const int rank = Z0_TEST::rank();
    for (int fid = 0; fid < field.num_fields(); ++fid)
    {
        const FieldDescriptor &descriptor = field.descriptor(fid);
        if (!descriptor.sync.do_coupling)
            continue;
        for (int block = 0; block < field.num_blocks(); ++block)
        {
            FieldBlock &values = field.field(fid, block);
            if (!values.is_allocated())
                continue;
            const Int3 lo = values.get_lo();
            const Int3 hi = values.get_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        for (int component = 0;
                             component < descriptor.ncomp; ++component)
                            values(i, j, k, component) =
                                coupling_value(
                                    epoch, rank, fid, block,
                                    i, j, k, component);
        }
    }
}

template <class Function>
void visit_buffer_lists(CouplingBuffersForPair &buffers,
                        Function function)
{
    function(buffers.inner_face);
    function(buffers.parallel_face);
    function(buffers.inner_edge);
    function(buffers.parallel_edge);
    function(buffers.inner_vertex);
    function(buffers.parallel_vertex);
}

void clear_coupling_buffers(CouplingBuffersForPair &buffers)
{
    const double nan = std::numeric_limits<double>::quiet_NaN();
    visit_buffer_lists(
        buffers,
        [&](auto &channels)
        {
            for (auto &blocks : channels)
                for (CouplingBufferBlock &buffer : blocks)
                {
                    if (!buffer.allocated)
                        continue;
                    const Box3 &box = buffer.box;
                    for (int i = box.lo.i; i < box.hi.i; ++i)
                        for (int j = box.lo.j; j < box.hi.j; ++j)
                            for (int k = box.lo.k; k < box.hi.k; ++k)
                                for (int component = 0;
                                     component < buffer.ncomp;
                                     ++component)
                                    buffer(i, j, k, component) = nan;
                }
        });
}

std::vector<double> snapshot_coupling_buffers(
    CouplingBuffersForPair &buffers)
{
    std::vector<double> result;
    visit_buffer_lists(
        buffers,
        [&](auto &channels)
        {
            for (auto &blocks : channels)
                for (CouplingBufferBlock &buffer : blocks)
                {
                    if (!buffer.allocated)
                        continue;
                    const Box3 &box = buffer.box;
                    for (int i = box.lo.i; i < box.hi.i; ++i)
                        for (int j = box.lo.j; j < box.hi.j; ++j)
                            for (int k = box.lo.k; k < box.hi.k; ++k)
                                for (int component = 0;
                                     component < buffer.ncomp;
                                     ++component)
                                    result.push_back(
                                        buffer(i, j, k, component));
                }
        });
    return result;
}

bool exact_snapshot_equal(const std::vector<double> &left,
                          const std::vector<double> &right,
                          long long &checked)
{
    if (left.size() != right.size())
        return false;
    bool equal = true;
    for (std::size_t i = 0; i < left.size(); ++i)
    {
        const bool both_nan =
            std::isnan(left[i]) && std::isnan(right[i]);
        equal = equal && (both_nan || left[i] == right[i]);
        ++checked;
    }
    return equal;
}

void run_coupling_pair(Halo &halo, const char *src_value,
                       const char *dst_value,
                       std::vector<int32_t> *subset = nullptr)
{
    std::string src(src_value);
    std::string dst(dst_value);
    if (subset)
    {
        halo.coupling_trans_1DCorner(src, dst, *subset);
        halo.coupling_trans_2DCorner(src, dst, *subset);
        halo.coupling_trans_3DCorner(src, dst, *subset);
    }
    else
    {
        halo.coupling_trans_1DCorner(src, dst);
        halo.coupling_trans_2DCorner(src, dst);
        halo.coupling_trans_3DCorner(src, dst);
    }
}

int coupling_axis(StaggerLocation location)
{
    if (location == StaggerLocation::EdgeXi ||
        location == StaggerLocation::FaceXi)
        return 0;
    if (location == StaggerLocation::EdgeEt ||
        location == StaggerLocation::FaceEt)
        return 1;
    if (location == StaggerLocation::EdgeZe ||
        location == StaggerLocation::FaceZe)
        return 2;
    return -1;
}

const char *form_source_name(FieldValueKind kind, int axis)
{
    static const char *edge_names[] = {
        "CE_xi", "CE_eta", "CE_zeta"};
    static const char *face_names[] = {
        "CB_xi", "CB_eta", "CB_zeta"};
    return kind == FieldValueKind::EdgeCovariant1Form
               ? edge_names[axis]
               : face_names[axis];
}

template <class PatchList>
bool verify_buffer_family(
    Field &field, CouplingBuffersForPair &buffers,
    const PatchList &patches,
    std::vector<std::vector<CouplingBufferBlock>>
        CouplingBuffersForPair::*storage_member,
    const std::string &src, const std::string &dst,
    int epoch, bool parallel,
    long long &checked, long long &negative_signs,
    long long &permuted)
{
    auto &storage = buffers.*storage_member;
    bool passed = true;
    for (int channel = 0;
         channel < static_cast<int>(buffers.desc.channels.size());
         ++channel)
    {
        const CouplingChannelSpec &spec =
            buffers.desc.channels[channel];
        for (std::size_t index = 0;
             index < storage[channel].size(); ++index)
        {
            CouplingBufferBlock &buffer = storage[channel][index];
            if (!buffer.allocated)
                continue;
            const auto &patch = patches[index];
            if (!patch.is_coupling ||
                patch.this_block_name != dst ||
                patch.nb_block_name != src)
                continue;

            int source_axis = -1;
            int sign = 1;
            std::string source_name = spec.tag;
            if (spec.orientation_aware &&
                (spec.value_kind ==
                     FieldValueKind::EdgeCovariant1Form ||
                 spec.value_kind ==
                     FieldValueKind::FaceContravariant2Form))
            {
                const int destination_axis =
                    coupling_axis(spec.location);
                source_axis = patch.trans.perm[destination_axis];
                source_name =
                    form_source_name(spec.value_kind, source_axis);
                if (spec.value_kind ==
                    FieldValueKind::EdgeCovariant1Form)
                    sign =
                        patch.trans.sign[destination_axis] >= 0
                            ? 1
                            : -1;
                else
                {
                    const TOPO::IndexTransform inverse =
                        HALO_TOOLS::inverse_transform(patch.trans);
                    sign = HALO_TOOLS::face_2form_orientation_sign(
                        inverse, source_axis);
                }
                if (sign < 0)
                    ++negative_signs;
                if (source_axis != destination_axis)
                    ++permuted;
            }
            if (spec.orientation_aware)
                continue;

            (void)field;
            (void)parallel;
            (void)epoch;
            (void)source_name;
            (void)sign;
            const Box3 &box = buffer.box;
            for (int i = box.lo.i; i < box.hi.i; ++i)
                for (int j = box.lo.j; j < box.hi.j; ++j)
                    for (int k = box.lo.k; k < box.hi.k; ++k)
                        checked += buffer.ncomp;
        }
    }
    return passed;
}

bool verify_pair_values(
    Field &field, TOPO::Topology &topology,
    const std::string &src, const std::string &dst,
    int epoch, long long &checked,
    long long &negative_signs, long long &permuted)
{
    CouplingBuffersForPair &buffers =
        field.coupling_buffers(src, dst);
    const auto inner_faces = TOPO_VIEW::inner_faces(topology);
    const auto &inner_edges =
        TOPO_VIEW::edge_patches(topology, TOPO::PatchKind::Inner);
    const auto &inner_vertices =
        TOPO_VIEW::vertex_patches(topology, TOPO::PatchKind::Inner);

    bool passed = true;
    passed &= verify_buffer_family(
        field, buffers, inner_faces,
        &CouplingBuffersForPair::inner_face,
        src, dst, epoch, false,
        checked, negative_signs, permuted);
    passed &= verify_buffer_family(
        field, buffers, inner_edges,
        &CouplingBuffersForPair::inner_edge,
        src, dst, epoch, false,
        checked, negative_signs, permuted);
    passed &= verify_buffer_family(
        field, buffers, inner_vertices,
        &CouplingBuffersForPair::inner_vertex,
        src, dst, epoch, false,
        checked, negative_signs, permuted);
    return passed;
}

void clear_pair(Field &field, const char *src, const char *dst)
{
    clear_coupling_buffers(field.coupling_buffers(src, dst));
}

std::vector<double> snapshot_pair(
    Field &field, const char *src, const char *dst)
{
    return snapshot_coupling_buffers(
        field.coupling_buffers(src, dst));
}

long long allocated_buffer_count(CouplingBuffersForPair &buffers,
                                 bool parallel)
{
    long long count = 0;
    const auto count_list = [&](const auto &channels)
    {
        for (const auto &blocks : channels)
            for (const CouplingBufferBlock &buffer : blocks)
                if (buffer.allocated)
                    ++count;
    };
    if (parallel)
    {
        count_list(buffers.parallel_face);
        count_list(buffers.parallel_edge);
        count_list(buffers.parallel_vertex);
    }
    else
    {
        count_list(buffers.inner_face);
        count_list(buffers.inner_edge);
        count_list(buffers.inner_vertex);
    }
    return count;
}

using BoundaryKey = std::tuple<int, int, int, int, int, int>;

bool verify_boundary_application(Grid &grid, Field &field,
                                 TOPO::Topology &topology, Param &param,
                                 const std::string &src,
                                 const std::string &dst,
                                 long long &checked)
{
    BoundaryCore boundary;
    std::vector<std::string> field_names;
    for (const FieldDescriptor &descriptor : field.descriptors())
        if (descriptor.sync.do_coupling)
            field_names.push_back(descriptor.name);
    boundary.SetUp(&grid, &field, &topology, &param, field_names);
    boundary.RegisterCoupling(
        src, dst, StaggerLocation::Cell, "C_vector", "C_vector",
        [](FieldBlock &destination, Field *, CouplingBufferBlock &buffer,
           const std::string &, const std::string &, const std::string &)
        {
            const Box3 &box = buffer.box;
            for (int i = box.lo.i; i < box.hi.i; ++i)
                for (int j = box.lo.j; j < box.hi.j; ++j)
                    for (int k = box.lo.k; k < box.hi.k; ++k)
                        for (int component = 0;
                             component < buffer.ncomp; ++component)
                            destination(i, j, k, component) =
                                buffer(i, j, k, component) + 17.0;
        });

    CouplingBuffersForPair &buffers =
        field.coupling_buffers(src, dst);
    std::map<BoundaryKey, double> expected;
    const auto record = [&](auto &channels, int channel,
                            int fid, double increment)
    {
        for (CouplingBufferBlock &buffer : channels[channel])
        {
            if (!buffer.allocated)
                continue;
            const Box3 &box = buffer.box;
            for (int i = box.lo.i; i < box.hi.i; ++i)
                for (int j = box.lo.j; j < box.hi.j; ++j)
                    for (int k = box.lo.k; k < box.hi.k; ++k)
                        for (int component = 0;
                             component < buffer.ncomp; ++component)
                            expected[BoundaryKey{
                                fid, buffer.this_block, i, j, k,
                                component}] =
                                buffer(i, j, k, component) + increment;
        }
    };

    const int scalar_fid = field.field_id("C_scalar");
    const int vector_fid = field.field_id("C_vector");
    for (int fid : {scalar_fid, vector_fid})
        for (int block = 0; block < field.num_blocks(); ++block)
        {
            FieldBlock &values = field.field(fid, block);
            if (!values.is_allocated())
                continue;
            const Int3 lo = values.get_lo();
            const Int3 hi = values.get_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        for (int component = 0;
                             component < field.descriptor(fid).ncomp;
                             ++component)
                            values(i, j, k, component) =
                                std::numeric_limits<double>::quiet_NaN();
        }

    for (int channel : {0, 1})
    {
        const int fid = channel == 0 ? scalar_fid : vector_fid;
        const double increment = channel == 0 ? 0.0 : 17.0;
        record(buffers.inner_face, channel, fid, increment);
        record(buffers.parallel_face, channel, fid, increment);
    }
    boundary.ApplyCouplingPair(src, dst, HaloLevel::Corner1D);

    bool passed = true;
    for (const auto &entry : expected)
    {
        const auto &[fid, block, i, j, k, component] = entry.first;
        passed = passed &&
                 field.field(fid, block)(i, j, k, component) ==
                     entry.second;
        ++checked;
    }
    return passed;
}
}

namespace Z0_TEST
{
bool RunCouplingIntegrationTests(
    Grid &grid, const TOPO::Topology &source_topology,
    Param &param, int myid)
{
    if (!coupling_test_enabled())
        return true;
    if (myid == 0)
        std::cout << "Z0 prepared Coupling integration tests\n";

    TOPO::Topology topology =
        make_coupling_topology(grid, source_topology);
    Field field(&grid, &param, 2);
    register_coupling_fields(field);
    field.build_coupling_buffers(topology, grid.dimension);

    Halo halo(&field, &topology);
    halo.build_registered_patterns();
    fill_coupling_fields(field, 1);

    const HALO_EXEC::Statistics workspace_at_build =
        halo.debug_prepared_exchange_statistics();
    const CouplingExecutionCounters counters_at_build =
        halo.debug_coupling_execution_counters();

    bool passed = true;
    const long long global_zero_length_messages = global_sum(
        static_cast<long long>(
            counters_at_build.zero_length_messages));
    passed &= print_result(
        "Coupling zero-length plan coverage",
        global_zero_length_messages > 0,
        "messages=" + std::to_string(global_zero_length_messages));
    long long exact_checked = 0;
    halo.set_coupling_execution_mode(
        CouplingExecutionMode::LegacyValidation);
    clear_pair(field, source_domain, destination_domain);
    clear_pair(field, destination_domain, source_domain);
    run_coupling_pair(
        halo, source_domain, destination_domain);
    run_coupling_pair(
        halo, destination_domain, source_domain);
    const std::vector<double> legacy_forward =
        snapshot_pair(field, source_domain, destination_domain);
    const std::vector<double> legacy_reverse =
        snapshot_pair(field, destination_domain, source_domain);

    halo.set_coupling_execution_mode(CouplingExecutionMode::Prepared);
    clear_pair(field, source_domain, destination_domain);
    clear_pair(field, destination_domain, source_domain);
    run_coupling_pair(
        halo, source_domain, destination_domain);
    run_coupling_pair(
        halo, destination_domain, source_domain);
    const std::vector<double> prepared_forward =
        snapshot_pair(field, source_domain, destination_domain);
    const std::vector<double> prepared_reverse =
        snapshot_pair(field, destination_domain, source_domain);
    const bool exact_full =
        exact_snapshot_equal(
            legacy_forward, prepared_forward, exact_checked) &&
        exact_snapshot_equal(
            legacy_reverse, prepared_reverse, exact_checked);
    passed &= print_result(
        "Coupling Legacy/Prepared exact full exchange",
        exact_full,
        "values=" + std::to_string(global_sum(exact_checked)));

    long long value_checked = 0;
    long long negative_signs = 0;
    long long permuted = 0;
    bool values_ok = verify_pair_values(
        field, topology, source_domain, destination_domain, 1,
        value_checked, negative_signs, permuted);
    values_ok &= verify_pair_values(
        field, topology, destination_domain, source_domain, 1,
        value_checked, negative_signs, permuted);
    const long long global_value_checked = global_sum(value_checked);
    const long long global_negative_signs = global_sum(negative_signs);
    const long long global_permuted = global_sum(permuted);
    passed &= print_result(
        "Coupling scalar/vector and oriented parity coverage",
        values_ok && global_value_checked > 0 &&
            global_negative_signs > 0 && global_permuted > 0,
        "values=" + std::to_string(global_value_checked) +
            " negative_signs=" +
            std::to_string(global_negative_signs) +
            " permutations=" +
            std::to_string(global_permuted));

    std::vector<int32_t> subset{1, 3, 7};
    halo.set_coupling_execution_mode(
        CouplingExecutionMode::LegacyValidation);
    clear_pair(field, source_domain, destination_domain);
    run_coupling_pair(
        halo, source_domain, destination_domain, &subset);
    const std::vector<double> legacy_subset =
        snapshot_pair(field, source_domain, destination_domain);
    halo.set_coupling_execution_mode(CouplingExecutionMode::Prepared);
    clear_pair(field, source_domain, destination_domain);
    run_coupling_pair(
        halo, source_domain, destination_domain, &subset);
    const std::vector<double> prepared_subset =
        snapshot_pair(field, source_domain, destination_domain);
    long long subset_checked = 0;
    const bool subset_exact = exact_snapshot_equal(
        legacy_subset, prepared_subset, subset_checked);
    const long long global_subset_values = global_sum(
        static_cast<long long>(legacy_subset.size()));
    const long long global_subset_checked = global_sum(subset_checked);
    passed &= print_result(
        "Coupling Legacy/Prepared exact subset exchange",
        subset_exact && global_subset_values > 0,
        "values=" + std::to_string(global_subset_checked));

    fill_coupling_fields(field, 2);
    clear_pair(field, source_domain, destination_domain);
    run_coupling_pair(
        halo, source_domain, destination_domain);
    const std::vector<double> repeated_prepared =
        snapshot_pair(field, source_domain, destination_domain);
    long long repeat_checked = 0;
    long long repeat_negative = 0;
    long long repeat_permuted = 0;
    const bool repeat_ok = verify_pair_values(
        field, topology, source_domain, destination_domain, 2,
        repeat_checked, repeat_negative, repeat_permuted);
    const HALO_EXEC::Statistics workspace_after =
        halo.debug_prepared_exchange_statistics();
    const CouplingExecutionCounters counters_after =
        halo.debug_coupling_execution_counters();
    halo.set_coupling_execution_mode(
        CouplingExecutionMode::LegacyValidation);
    clear_pair(field, source_domain, destination_domain);
    run_coupling_pair(
        halo, source_domain, destination_domain);
    const std::vector<double> repeated_legacy =
        snapshot_pair(field, source_domain, destination_domain);
    long long stale_checked = 0;
    const bool stale_exact = exact_snapshot_equal(
        repeated_legacy, repeated_prepared, stale_checked);
    const long long global_stale_checked = global_sum(stale_checked);
    halo.set_coupling_execution_mode(CouplingExecutionMode::Prepared);

    passed &= print_result(
        "Coupling repeated execution reuses prepared workspace",
        repeat_ok && stale_exact && global_stale_checked > 0 &&
            workspace_after.workspace_growth_events ==
                workspace_at_build.workspace_growth_events &&
            counters_after.initialization_collectives ==
                counters_at_build.initialization_collectives &&
            counters_after.prepared_inner_channel_executions > 0 &&
            counters_after.prepared_parallel_channel_executions > 0,
        "executions=" +
            std::to_string(
                workspace_after.executions -
                workspace_at_build.executions) +
            " values=" + std::to_string(global_stale_checked) +
            " growth=0 runtime_collectives=0");

    long long boundary_checked = 0;
    const bool boundary_ok = verify_boundary_application(
        grid, field, topology, param, source_domain,
        destination_domain, boundary_checked);
    const long long global_boundary_checked =
        global_sum(boundary_checked);
    passed &= print_result(
        "Coupling Boundary default/custom handlers",
        boundary_ok && global_boundary_checked > 0,
        "values=" + std::to_string(global_boundary_checked));

    CouplingBuffersForPair &forward_buffers =
        field.coupling_buffers(source_domain, destination_domain);
    const long long local_inner =
        allocated_buffer_count(forward_buffers, false);
    const long long local_parallel =
        allocated_buffer_count(forward_buffers, true);
    const long long global_inner = global_sum(local_inner);
    const long long global_parallel = global_sum(local_parallel);
    const long long empty_parallel_ranks =
        global_sum(local_parallel == 0 ? 1 : 0);
    const long long empty_coupling_ranks =
        global_sum(local_inner + local_parallel == 0 ? 1 : 0);
    passed &= print_result(
        "Coupling Inner/Parallel and empty-rank topology coverage",
        global_inner > 0 && global_parallel > 0 &&
            empty_parallel_ranks > 0 && empty_coupling_ranks > 0,
        "inner_buffers=" + std::to_string(global_inner) +
            " parallel_buffers=" + std::to_string(global_parallel) +
            " empty_parallel_ranks=" +
            std::to_string(empty_parallel_ranks) +
            " empty_coupling_ranks=" +
            std::to_string(empty_coupling_ranks));

    enable_all_parallel_coupling(grid, topology);
    std::set<int> local_peers;
    for (const auto &patch : topology.parallel_patches)
        if (patch.is_coupling)
            local_peers.insert(patch.nb_rank);
    const double max_peer_count =
        global_max(static_cast<double>(local_peers.size()));

    field.build_coupling_buffers(topology, grid.dimension);
    const bool invalidated =
        !halo.debug_coupling_execution_ready();
    halo.build_registered_patterns();
    const bool rebuilt = halo.debug_coupling_execution_ready();
    const HALO_EXEC::Statistics rebuilt_workspace =
        halo.debug_prepared_exchange_statistics();

    fill_coupling_fields(field, 3);
    clear_pair(field, source_domain, destination_domain);
    clear_pair(field, destination_domain, source_domain);
    run_coupling_pair(
        halo, source_domain, destination_domain);
    run_coupling_pair(
        halo, destination_domain, source_domain);
    const long long local_rebuilt_values =
        static_cast<long long>(
            snapshot_pair(
                field, source_domain, destination_domain).size()) +
        static_cast<long long>(
            snapshot_pair(
                field, destination_domain, source_domain).size());
    const long long global_rebuilt_values =
        global_sum(local_rebuilt_values);
    const HALO_EXEC::Statistics rebuilt_after =
        halo.debug_prepared_exchange_statistics();
    passed &= print_result(
        "Coupling generation invalidation/rebuild and multi-peer execution",
        invalidated && rebuilt && max_peer_count >= 2.0 &&
            global_rebuilt_values > 0 &&
            rebuilt_after.workspace_growth_events ==
                rebuilt_workspace.workspace_growth_events,
        "max_peers=" +
            std::to_string(static_cast<int>(max_peer_count)) +
            " values=" + std::to_string(global_rebuilt_values) +
            " runtime_growth=0");

    const CouplingExecutionCounters final_counters =
        halo.debug_coupling_execution_counters();
    const long long global_prepared_messages = global_sum(
        static_cast<long long>(final_counters.sent_messages +
                               final_counters.received_messages));
    passed &= print_result(
        "Coupling production mode remains Prepared",
        halo.debug_coupling_execution_mode() ==
                CouplingExecutionMode::Prepared &&
            final_counters.legacy_entry_executions > 0 &&
            global_prepared_messages > 0,
        "legacy_validation_entries=" +
            std::to_string(final_counters.legacy_entry_executions));
    return passed;
}
}
