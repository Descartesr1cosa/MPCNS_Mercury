#include "Z0_Tests.h"
#include "Z0_TestCommon.h"

#include "0_basic/MPI_WRAPPER.h"
#include "1_grid/1_MPCNS_Grid.h"
#include "2_topology/Topology.h"
#include "3_field/Field.h"
#include "4_halo/Halo.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace
{
using StateKey = std::tuple<int, int, int, int, int, int>;
using State = std::map<StateKey, double>;

struct Summary
{
    long long checked = 0;
    long long failed = 0;
    long long reported = 0;
};

void check(Summary &summary, int rank, const std::string &label,
           bool condition)
{
    ++summary.checked;
    if (condition)
        return;
    ++summary.failed;
    if (summary.reported++ < 12)
        std::cerr << "[PreparedFace2FormExecution][rank " << rank
                  << "] " << label << " failed\n";
}

void seed_bface(
    Field &field, const std::array<std::string, 3> &fields,
    int rank, int ncomp)
{
    for (int axis = 0; axis < 3; ++axis)
    {
        const int fid = field.field_id(fields[axis]);
        for (int block = 0; block < field.num_blocks(); ++block)
        {
            FieldBlock &fb = field.field(fid, block);
            const Int3 lo = fb.get_lo();
            const Int3 hi = fb.get_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        for (int component = 0; component < ncomp;
                             ++component)
                            fb(i, j, k, component) =
                                Z0_TEST::unique_code(
                                    rank, block, i, j, k, component) +
                                static_cast<double>(axis) * 0.125;
        }
    }
}

State snapshot_bface(
    Field &field, const std::array<std::string, 3> &fields,
    int ncomp)
{
    State state;
    for (int axis = 0; axis < 3; ++axis)
    {
        const int fid = field.field_id(fields[axis]);
        for (int block = 0; block < field.num_blocks(); ++block)
        {
            const FieldBlock &fb = field.field(fid, block);
            const Int3 lo = fb.get_lo();
            const Int3 hi = fb.get_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        for (int component = 0; component < ncomp;
                             ++component)
                            state[{axis, block, i, j, k, component}] =
                                fb(i, j, k, component);
        }
    }
    return state;
}

void compare_state(
    Summary &summary, int rank, const std::string &label,
    const State &expected, const State &actual)
{
    check(summary, rank, label + ".size", expected.size() == actual.size());
    for (const auto &entry : expected)
    {
        const auto found = actual.find(entry.first);
        check(summary, rank, label + ".present",
              found != actual.end());
        if (found != actual.end())
            check(summary, rank, label + ".value",
                  found->second == entry.second);
    }
}

bool registered(
    const Halo &halo, LegacyOrientedFormKind kind,
    const char *group, HaloLevel stage)
{
    const auto triplets = halo.debug_oriented_form_triplets();
    return std::any_of(
        triplets.begin(), triplets.end(),
        [&](const LegacyOrientedTripletSnapshot &triplet)
        {
            return triplet.form_kind == kind &&
                   triplet.group_name == group &&
                   triplet.orientation_aware &&
                   static_cast<int>(triplet.level) >=
                       static_cast<int>(stage);
        });
}

template <class Counters>
bool same_data_per_sync(
    const Counters &legacy_before, const Counters &legacy_after,
    const Counters &prepared_before, const Counters &prepared_after,
    std::uint64_t prepared_runs)
{
    const std::uint64_t legacy_send =
        legacy_after.send_messages - legacy_before.send_messages;
    const std::uint64_t legacy_recv =
        legacy_after.recv_messages - legacy_before.recv_messages;
    const std::uint64_t legacy_payload =
        legacy_after.sent_doubles - legacy_before.sent_doubles;
    const std::uint64_t legacy_received =
        legacy_after.received_doubles - legacy_before.received_doubles;
    const bool same =
        prepared_after.send_messages - prepared_before.send_messages ==
            prepared_runs * legacy_send &&
        prepared_after.recv_messages - prepared_before.recv_messages ==
            prepared_runs * legacy_recv &&
        prepared_after.sent_doubles - prepared_before.sent_doubles ==
            prepared_runs * legacy_payload &&
        prepared_after.received_doubles - prepared_before.received_doubles ==
            prepared_runs * legacy_received;
    if (!same)
        std::cerr << "[PreparedFace2FormExecution] parity legacy="
                  << legacy_send << "/" << legacy_recv << "/"
                  << legacy_payload << "/" << legacy_received
                  << " prepared="
                  << (prepared_after.send_messages - prepared_before.send_messages)
                  << "/" << (prepared_after.recv_messages - prepared_before.recv_messages)
                  << "/" << (prepared_after.sent_doubles - prepared_before.sent_doubles)
                  << "/" << (prepared_after.received_doubles - prepared_before.received_doubles)
                  << "\n";
    return same;
}
} // namespace

namespace Z0_TEST
{
bool RunPreparedFace2FormExecutionTests(
    const Grid &grid, const TOPO::Topology &topology,
    Field &field, Halo &halo, int rank)
{
    (void)grid;
    (void)topology;
    int size = 1;
    PARALLEL::mpi_size(&size);
    Summary summary;
    const std::array<std::string, 3> bfields = {
        "B_xi", "B_eta", "B_zeta"};

    const auto requests = field.halo_requests();
    const auto request = std::find_if(
        requests.begin(), requests.end(),
        [](const FieldHaloRequest &request)
        {
            const std::string group = request.sync_group.empty()
                ? request.field_name : request.sync_group;
            return request.do_halo && request.orientation_aware &&
                   request.value_kind ==
                       FieldValueKind::FaceContravariant2Form &&
                   group == "Bface";
        });

    check(summary, rank, "face2_edge_default_legacy",
          halo.debug_face_2form_edge_execution_mode() ==
              Face2FormEdgeExecutionMode::Legacy);
    check(summary, rank, "face2_vertex_default_legacy",
          halo.debug_face_2form_vertex_execution_mode() ==
              Face2FormVertexExecutionMode::Legacy);
    check(summary, rank, "edge1_edge_default_legacy",
          halo.debug_edge_1form_edge_execution_mode() ==
              Edge1FormEdgeExecutionMode::Legacy);
    check(summary, rank, "edge1_vertex_default_legacy",
          halo.debug_edge_1form_vertex_execution_mode() ==
              Edge1FormVertexExecutionMode::Legacy);

    const bool face_edge_registered = registered(
        halo, LegacyOrientedFormKind::Face2Form, "Bface",
        HaloLevel::Corner2D);
    const bool face_vertex_registered = registered(
        halo, LegacyOrientedFormKind::Face2Form, "Bface",
        HaloLevel::Corner3D);
    const bool edge_edge_registered = registered(
        halo, LegacyOrientedFormKind::Edge1Form, "Eedge",
        HaloLevel::Corner2D);
    const bool edge_vertex_registered = registered(
        halo, LegacyOrientedFormKind::Edge1Form, "Eedge",
        HaloLevel::Corner3D);

    check(summary, rank, "face2_edge_explicit_legacy",
          halo.deactivate_prepared_face_2form_edge_execution());
    check(summary, rank, "face2_vertex_explicit_legacy",
          halo.deactivate_prepared_face_2form_vertex_execution());
    check(summary, rank, "edge1_edge_explicit_legacy",
          halo.deactivate_prepared_edge_1form_edge_execution());
    check(summary, rank, "edge1_vertex_explicit_legacy",
          halo.deactivate_prepared_edge_1form_vertex_execution());

    if (request == requests.end())
    {
        check(summary, rank, "orientation_disabled_face2_edge_rejected",
              !halo.activate_prepared_face_2form_edge_execution());
        check(summary, rank, "orientation_disabled_face2_vertex_rejected",
              !halo.activate_prepared_face_2form_vertex_execution());
        const long long total = global_sum(summary.checked);
        const long long failures = global_sum(summary.failed);
        std::ostringstream detail;
        detail << "orientation_aware=0 checks=" << total
               << " failed=" << failures;
        return print_result(
            "Prepared Face2Form Edge/Vertex execution",
            summary.failed == 0, detail.str());
    }

    check(summary, rank, "prepare_face2_edge",
          halo.prepare_edge_face_2form_metadata());
    check(summary, rank, "prepare_face2_vertex",
          halo.prepare_vertex_face_2form_metadata());
    check(summary, rank, "prepare_edge1_edge",
          halo.prepare_edge_1form_corner_metadata());
    check(summary, rank, "prepare_edge1_vertex",
          halo.prepare_vertex_edge_1form_metadata());

    const auto run_bface =
        [&](bool include_vertex, bool execute_edge1)
        {
            seed_bface(field, bfields, rank, request->ncomp);
            if (execute_edge1)
                halo.sync_group("Eedge", HaloLevel::Corner1D);
            halo.sync_group("Bface", HaloLevel::Corner1D);
            if (execute_edge1)
                halo.sync_group("Eedge", HaloLevel::Corner2D);
            halo.sync_group("Bface", HaloLevel::Corner2D);
            if (include_vertex)
            {
                if (execute_edge1)
                    halo.sync_group("Eedge", HaloLevel::Corner3D);
                halo.sync_group("Bface", HaloLevel::Corner3D);
            }
            return snapshot_bface(field, bfields, request->ncomp);
        };

    State legacy_edge;
    if (face_edge_registered)
    {
        const auto legacy_before =
            halo.debug_face_2form_edge_execution_counters();
        legacy_edge = run_bface(false, false);
        const auto legacy_after =
            halo.debug_face_2form_edge_execution_counters();

        if (size > 1)
        {
            const Face2FormEdgeExecutionMode mismatch =
                rank == 0 ? Face2FormEdgeExecutionMode::Prepared
                          : Face2FormEdgeExecutionMode::Legacy;
            check(summary, rank, "face2_edge_mixed_rank_rejected",
                  !halo.set_face_2form_edge_execution_mode(mismatch));
        }
        check(summary, rank, "face2_edge_activate",
              halo.activate_prepared_face_2form_edge_execution());
        const auto prepared_before =
            halo.debug_face_2form_edge_execution_counters();
        const std::uint64_t prepares_before =
            halo.debug_prepared_edge_face_2form_metadata_exchange_count();
        const State prepared = run_bface(false, false);
        const State repeated = run_bface(false, false);
        const auto prepared_after =
            halo.debug_face_2form_edge_execution_counters();
        compare_state(summary, rank, "face2_edge_prepared", legacy_edge,
                      prepared);
        compare_state(summary, rank, "face2_edge_repeat", legacy_edge,
                      repeated);
        check(summary, rank, "face2_edge_prepared_runs",
              prepared_after.prepared_executions ==
                  prepared_before.prepared_executions + 2);
        check(summary, rank, "face2_edge_no_runtime_alltoall",
              prepared_after.legacy_metadata_alltoall ==
                  prepared_before.legacy_metadata_alltoall);
        check(summary, rank, "face2_edge_no_runtime_alltoallv",
              prepared_after.legacy_metadata_alltoallv ==
                  prepared_before.legacy_metadata_alltoallv);
        check(summary, rank, "face2_edge_reuses_prepare",
              halo.debug_prepared_edge_face_2form_metadata_exchange_count() ==
                  prepares_before);
        check(summary, rank, "face2_edge_message_payload_parity",
              same_data_per_sync(
                  legacy_before, legacy_after, prepared_before,
                  prepared_after, 2));
        check(summary, rank, "face2_edge_deactivate",
              halo.deactivate_prepared_face_2form_edge_execution());
    }

    State legacy_vertex;
    if (face_vertex_registered)
    {
        const auto legacy_before =
            halo.debug_face_2form_vertex_execution_counters();
        legacy_vertex = run_bface(true, false);
        const auto legacy_after =
            halo.debug_face_2form_vertex_execution_counters();
        check(summary, rank, "face2_vertex_activate",
              halo.activate_prepared_face_2form_vertex_execution());
        const auto prepared_before =
            halo.debug_face_2form_vertex_execution_counters();
        const std::uint64_t prepares_before =
            halo.debug_prepared_vertex_face_2form_metadata_exchange_count();
        const State prepared = run_bface(true, false);
        const State repeated = run_bface(true, false);
        const auto prepared_after =
            halo.debug_face_2form_vertex_execution_counters();
        compare_state(summary, rank, "face2_vertex_prepared", legacy_vertex,
                      prepared);
        compare_state(summary, rank, "face2_vertex_repeat", legacy_vertex,
                      repeated);
        check(summary, rank, "face2_vertex_prepared_runs",
              prepared_after.prepared_executions ==
                  prepared_before.prepared_executions + 2);
        check(summary, rank, "face2_vertex_no_runtime_alltoall",
              prepared_after.legacy_metadata_alltoall ==
                  prepared_before.legacy_metadata_alltoall);
        check(summary, rank, "face2_vertex_no_runtime_alltoallv",
              prepared_after.legacy_metadata_alltoallv ==
                  prepared_before.legacy_metadata_alltoallv);
        check(summary, rank, "face2_vertex_reuses_prepare",
              halo.debug_prepared_vertex_face_2form_metadata_exchange_count() ==
                  prepares_before);
        check(summary, rank, "face2_vertex_message_payload_parity",
              same_data_per_sync(
                  legacy_before, legacy_after, prepared_before,
                  prepared_after, 2));
        check(summary, rank, "face2_vertex_deactivate",
              halo.deactivate_prepared_face_2form_vertex_execution());
    }

    if (face_vertex_registered && edge_edge_registered &&
        edge_vertex_registered)
    {
        const State all_legacy = run_bface(true, true);
        check(summary, rank, "face2_edge_only_activate",
              halo.activate_prepared_face_2form_edge_execution());
        const State face2_edge_only = run_bface(true, true);
        compare_state(summary, rank, "face2_edge_only_bface", all_legacy,
                      face2_edge_only);
        check(summary, rank, "face2_vertex_after_edge_activate",
              halo.activate_prepared_face_2form_vertex_execution());
        const State face2_only = run_bface(true, true);
        compare_state(summary, rank, "face2_only_bface", all_legacy,
                      face2_only);
        check(summary, rank, "all_edge1_edge_activate",
              halo.activate_prepared_edge_1form_edge_execution());
        check(summary, rank, "all_edge1_vertex_activate",
              halo.activate_prepared_edge_1form_vertex_execution());
        const auto e1e_before =
            halo.debug_edge_1form_edge_execution_counters();
        const auto e1v_before =
            halo.debug_edge_1form_vertex_execution_counters();
        const auto f2e_before =
            halo.debug_face_2form_edge_execution_counters();
        const auto f2v_before =
            halo.debug_face_2form_vertex_execution_counters();
        const State combined = run_bface(true, true);
        const auto e1e_after =
            halo.debug_edge_1form_edge_execution_counters();
        const auto e1v_after =
            halo.debug_edge_1form_vertex_execution_counters();
        const auto f2e_after =
            halo.debug_face_2form_edge_execution_counters();
        const auto f2v_after =
            halo.debug_face_2form_vertex_execution_counters();
        compare_state(summary, rank, "all_prepared_bface", all_legacy,
                      combined);
        check(summary, rank, "all_prepared_edge1_edge_executes",
              e1e_after.prepared_executions ==
                  e1e_before.prepared_executions + 1);
        check(summary, rank, "all_prepared_edge1_vertex_executes",
              e1v_after.prepared_executions ==
                  e1v_before.prepared_executions + 1);
        check(summary, rank, "all_prepared_face2_edge_executes",
              f2e_after.prepared_executions ==
                  f2e_before.prepared_executions + 1);
        check(summary, rank, "all_prepared_face2_vertex_executes",
              f2v_after.prepared_executions ==
                  f2v_before.prepared_executions + 1);
        check(summary, rank, "all_prepared_face2_collectives_zero",
              f2e_after.legacy_metadata_alltoall ==
                  f2e_before.legacy_metadata_alltoall &&
              f2v_after.legacy_metadata_alltoall ==
                  f2v_before.legacy_metadata_alltoall);
        halo.build_registered_patterns();
        check(summary, rank, "rebuild_invalidates_face2_edge_view",
              !halo.debug_prepared_face_2form_edge_execution_active());
        check(summary, rank, "rebuild_invalidates_face2_vertex_view",
              !halo.debug_prepared_face_2form_vertex_execution_active());
        check(summary, rank, "rebuild_invalidates_edge1_edge_view",
              !halo.debug_prepared_edge_1form_edge_execution_active());
        check(summary, rank, "rebuild_invalidates_edge1_vertex_view",
              !halo.debug_prepared_edge_1form_vertex_execution_active());
        check(summary, rank, "stale_face2_edge_activation_rejected",
              !halo.activate_prepared_face_2form_edge_execution());
        check(summary, rank, "stale_face2_vertex_activation_rejected",
              !halo.activate_prepared_face_2form_vertex_execution());
        check(summary, rank, "rebuild_prepare_face2_edge",
              halo.prepare_edge_face_2form_metadata());
        check(summary, rank, "rebuild_prepare_face2_vertex",
              halo.prepare_vertex_face_2form_metadata());
        check(summary, rank, "rebuild_reactivate_face2_edge",
              halo.activate_prepared_face_2form_edge_execution());
        check(summary, rank, "rebuild_reactivate_face2_vertex",
              halo.activate_prepared_face_2form_vertex_execution());
        check(summary, rank, "rebuild_deactivate_face2_edge",
              halo.deactivate_prepared_face_2form_edge_execution());
        check(summary, rank, "rebuild_deactivate_face2_vertex",
              halo.deactivate_prepared_face_2form_vertex_execution());

        check(summary, rank, "all_prepared_deactivate_edge1_edge",
              halo.deactivate_prepared_edge_1form_edge_execution());
        check(summary, rank, "all_prepared_deactivate_edge1_vertex",
              halo.deactivate_prepared_edge_1form_vertex_execution());
        check(summary, rank, "all_prepared_deactivate_face2_edge",
              halo.deactivate_prepared_face_2form_edge_execution());
        check(summary, rank, "all_prepared_deactivate_face2_vertex",
              halo.deactivate_prepared_face_2form_vertex_execution());
    }

    const long long total = global_sum(summary.checked);
    const long long failures = global_sum(summary.failed);
    const long long edge_send = global_sum(
        static_cast<long long>(
            halo.debug_face_2form_edge_execution_counters().send_messages));
    const long long edge_recv = global_sum(
        static_cast<long long>(
            halo.debug_face_2form_edge_execution_counters().recv_messages));
    const long long vertex_send = global_sum(
        static_cast<long long>(
            halo.debug_face_2form_vertex_execution_counters().send_messages));
    const long long vertex_recv = global_sum(
        static_cast<long long>(
            halo.debug_face_2form_vertex_execution_counters().recv_messages));
    const long long edge_sent_doubles = global_sum(
        static_cast<long long>(
            halo.debug_face_2form_edge_execution_counters().sent_doubles));
    const long long edge_received_doubles = global_sum(
        static_cast<long long>(
            halo.debug_face_2form_edge_execution_counters().received_doubles));
    const long long vertex_sent_doubles = global_sum(
        static_cast<long long>(
            halo.debug_face_2form_vertex_execution_counters().sent_doubles));
    const long long vertex_received_doubles = global_sum(
        static_cast<long long>(
            halo.debug_face_2form_vertex_execution_counters().received_doubles));
    std::ostringstream detail;
    detail << "checks=" << total << " failed=" << failures
           << " face2_edge_registered=" << face_edge_registered
           << " face2_vertex_registered=" << face_vertex_registered
           << " face2_edge_messages=" << edge_send << "/" << edge_recv
           << " face2_vertex_messages=" << vertex_send << "/"
           << vertex_recv
           << " face2_edge_doubles=" << edge_sent_doubles << "/"
           << edge_received_doubles
           << " face2_vertex_doubles=" << vertex_sent_doubles << "/"
           << vertex_received_doubles
           << " prepared_runtime_metadata_collectives=0";
    return print_result(
        "Prepared Face2Form Edge/Vertex execution",
        summary.failed == 0, detail.str());
}
} // namespace Z0_TEST

