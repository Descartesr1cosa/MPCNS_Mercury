#pragma once
#include <cstdint>
#include <vector>

#include "0_basic/MPI_WRAPPER.h"
#include "2_topology/2_MPCNS_Topology_Equiv.h"
#include "3_field/2_MPCNS_Field.h"
#include "4_halo/Halo_EdgeOwner_Type.h"

namespace HALO_OWNER
{
    // ============================================================
    // edge-owner sync pattern
    // ============================================================

    struct EdgeOwnerLocalAliasItem
    {
        TOPO::EdgeLocalID owner; // local owner rep
        TOPO::EdgeLocalID rep;   // local non-owner rep
        int8_t sign;             // for 1-form; ignored by vec copy
    };

    struct EdgeOwnerSendItem
    {
        int tar_id;              // target rank id
        TOPO::EdgeKey key;       //
        TOPO::EdgeLocalID owner; // local owner rep
    };

    struct EdgeOwnerRecvItem
    {
        int tar_id;            // remote owner rank id
        TOPO::EdgeKey key;     //
        TOPO::EdgeLocalID rep; // local non-owner rep
        int8_t sign;           // for 1-form; ignored by vec copy
    };

    struct EdgeOwnerSyncPattern
    {
        std::vector<EdgeOwnerLocalAliasItem> local_alias;

        std::vector<EdgeOwnerSendItem> send_items;
        std::vector<EdgeOwnerRecvItem> recv_items;

        std::vector<int> send_counts;
        std::vector<int> recv_counts;
        std::vector<int> send_displs;
        std::vector<int> recv_displs;

        // ------------------------------------------------------------
        // runtime caches (allocated once after pattern build)
        // ------------------------------------------------------------
        std::vector<std::vector<double>> send_buf_cache;
        std::vector<std::vector<double>> recv_buf_cache;

        std::vector<MPI_Request> req_recv_cache;
        std::vector<MPI_Status> stat_recv_cache;

        std::vector<MPI_Request> req_send_cache;
        std::vector<MPI_Status> stat_send_cache;

        void clear()
        {
            local_alias.clear();
            send_items.clear();
            recv_items.clear();

            send_counts.clear();
            recv_counts.clear();
            send_displs.clear();
            recv_displs.clear();

            send_buf_cache.clear();
            recv_buf_cache.clear();
            req_recv_cache.clear();
            stat_recv_cache.clear();
            req_send_cache.clear();
            stat_send_cache.clear();
        }
    };

    // ============================================================
    // build
    // ============================================================

    void build_edge_owner_sync_pattern(
        const TOPO::TopologyEquiv &equiv,
        EdgeOwnerSyncPattern &pattern);

    // ============================================================
    // runtime sync
    // ============================================================

    void gather_local_owner_edges_sorted(
        const TOPO::TopologyEquiv &equiv,
        std::vector<TOPO::EdgeLocalID> &owner_edges_sorted);

    // Treat every component on edge as a 1-form component:
    // rep(comp) = sign * owner(comp)
    void sync_edge_1form(
        Field &fld,
        const IdTriplet &field_id,
        EdgeOwnerSyncPattern &pattern);

    // Treat every component on edge as a Cartesian vector component:
    // rep(comp) = owner(comp)
    void sync_edge_vec(
        Field &fld,
        const IdTriplet &field_id,
        EdgeOwnerSyncPattern &pattern);

    void pack_owner_edge_1form_local(
        Field &fld,
        const IdTriplet &field_id,
        const TOPO::TopologyEquiv &equiv,
        const std::vector<TOPO::EdgeLocalID> &owner_edges_sorted,
        std::vector<double> &buf_local);

    void unpack_owner_edge_1form_local(
        const std::vector<double> &buf_local,
        Field &fld,
        const IdTriplet &field_id,
        const TOPO::TopologyEquiv &equiv,
        const std::vector<TOPO::EdgeLocalID> &owner_edges_sorted,
        EdgeOwnerSyncPattern &pattern);

} // namespace HALO_OWNER
