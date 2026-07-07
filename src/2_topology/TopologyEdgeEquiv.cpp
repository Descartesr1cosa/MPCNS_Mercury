#include "2_topology/TopologyEquivDetail.h"

#include "0_basic/MPI_WRAPPER.h"
#include "1_grid/1_MPCNS_Grid.h"

#include <algorithm>
#include <stdexcept>
#include <unordered_map>

namespace TOPO::detail
{
    std::vector<EntityKey> collect_all_local_edges(
        Grid &grid,
        int my_rank,
        int dimension)
    {
        std::vector<EntityKey> edges;

        for (int ib = 0; ib < grid.nblock; ++ib)
        {
        const auto &blk = grid.grids(ib);

        // Xi
        for (int i = 0; i < blk.mx; ++i)
            for (int j = 0; j <= blk.my; ++j)
            for (int k = 0; k <= blk.mz; ++k)
            {
                edges.push_back(make_edge(my_rank, ib, i, j, k, EntityAxis::Xi));
            }

        if (dimension >= 2)
        {
            // Eta
            for (int i = 0; i <= blk.mx; ++i)
            for (int j = 0; j < blk.my; ++j)
                for (int k = 0; k <= blk.mz; ++k)
                {
                edges.push_back(make_edge(my_rank, ib, i, j, k, EntityAxis::Eta));
                }
        }

        if (dimension >= 3)
        {
            // Zeta
            for (int i = 0; i <= blk.mx; ++i)
            for (int j = 0; j <= blk.my; ++j)
                for (int k = 0; k < blk.mz; ++k)
                {
                edges.push_back(make_edge(my_rank, ib, i, j, k, EntityAxis::Zeta));
                }
        }
        }
        return edges;
    }

    
    static int node_equiv_count(
        const std::unordered_map<EntityKey, int, EntityKey::Hash> &rep_count,
        const EntityKey &node)
    {
        auto it = rep_count.find(node);
        return (it != rep_count.end()) ? it->second : 1;
    }

    
    bool edge_key_is_shared(
        const EdgeKey &key,
        const std::unordered_map<EntityKey, int, EntityKey::Hash> &rep_count)
    {
        return node_equiv_count(rep_count, key.a) > 1 ||
           node_equiv_count(rep_count, key.b) > 1;
    }

    
    void build_edge_equivalence_from_nodes(
        const std::vector<EntityKey> &all_local_edges,
        Topology &equiv,
        EdgeBuildScratch &scratch)
    {
        equiv.edges.local_to_qkey.clear();
        equiv.edges.local_to_qsign.clear();
        scratch.qkey_to_local_members.clear();

        for (const auto &e : all_local_edges)
        {
        int8_t sign_to_canonical = 0;
        EdgeKey key = make_edge_key(e, equiv.nodes.local_to_rep, sign_to_canonical);

        equiv.edges.local_to_qkey[e] = key;
        equiv.edges.local_to_qsign[e] = sign_to_canonical;
        scratch.qkey_to_local_members[key].push_back(e);
        }
    }

    
    void pack_edge_local(std::vector<int> &buf, const EntityKey &e)
    {
        buf.push_back(e.rank);
        buf.push_back(e.block);
        buf.push_back(e.i);
        buf.push_back(e.j);
        buf.push_back(e.k);
        buf.push_back(axis_number(e.axis));
    }

    
    EntityKey unpack_edge_local(const int *p)
    {
        return make_edge(p[0], p[1], p[2], p[3], p[4], entity_axis(p[5]));
    }

    
    void pack_edge_key(std::vector<int> &buf, const EdgeKey &key)
    {
        pack_node(buf, key.a);
        pack_node(buf, key.b);
    }

    
    EdgeKey unpack_edge_key(const int *p)
    {
        EntityKey a = unpack_node(p + 0);
        EntityKey b = unpack_node(p + 5);
        return EdgeKey{a, b};
    }

    
    void build_edge_entity_ids(Topology &equiv)
    {
        std::vector<EdgeKey> local_keys;
        local_keys.reserve(equiv.edges.local_to_qkey.size());
        for (const auto &[edge, key] : equiv.edges.local_to_qkey)
        {
        (void)edge;
        local_keys.push_back(key);
        }
        std::sort(local_keys.begin(), local_keys.end());
        local_keys.erase(std::unique(local_keys.begin(), local_keys.end()), local_keys.end());

        std::vector<int> send_buf;
        send_buf.reserve(local_keys.size() * 10);
        for (const EdgeKey &key : local_keys)
        pack_edge_key(send_buf, key);

        const std::vector<int> recv_buf =
        allgather_packed_records(send_buf, 10, "build_topology edge EntityId");
        std::vector<EdgeKey> global_keys;
        global_keys.reserve(recv_buf.size() / 10);
        for (std::size_t n = 0; n < recv_buf.size(); n += 10)
        global_keys.push_back(unpack_edge_key(recv_buf.data() + n));

        std::sort(global_keys.begin(), global_keys.end());
        global_keys.erase(std::unique(global_keys.begin(), global_keys.end()), global_keys.end());
        equiv.edges.qkey_to_qid.clear();
        for (std::size_t n = 0; n < global_keys.size(); ++n)
        equiv.edges.qkey_to_qid[global_keys[n]] = static_cast<int>(n);
    }

    
    void pack_edge_owner_candidate(
        std::vector<int> &buf,
        const EdgeKey &key,
        const EntityKey &e,
        int8_t sign)
    {
        // 10 ints for EdgeKey + 6 ints for EntityKey + 1 sign = 17 ints
        pack_edge_key(buf, key);
        pack_edge_local(buf, e);
        buf.push_back(static_cast<int>(sign));
    }

    
    void select_edge_owner_parallel(
        Topology &equiv,
        const EdgeBuildScratch &scratch,
        const std::unordered_map<EntityKey, int, EntityKey::Hash> &rep_count)
    {
        // ------------------------------------------------------------
        // 1) 打包本 rank 的所有 local edge candidates
        //    每个 local member 都作为一个候选发出去
        // ------------------------------------------------------------
        std::vector<int> send_buf;
        send_buf.reserve(1024);

        for (const auto &[key, members] : scratch.qkey_to_local_members)
        {
        if (!edge_key_is_shared(key, rep_count) && members.size() <= 1)
            continue;

        for (const auto &e : members)
        {
            const auto sign_it = equiv.edges.local_to_qsign.find(e);
            const int8_t sign = (sign_it != equiv.edges.local_to_qsign.end()) ? sign_it->second : int8_t{+1};
            pack_edge_owner_candidate(send_buf, key, e, sign);
        }
        }

        // ------------------------------------------------------------
        // 2) gather counts + bcast counts + allgatherv data
        // ------------------------------------------------------------
        int nrank = 1;
        PARALLEL::mpi_size(&nrank);

        int send_count = static_cast<int>(send_buf.size());

        std::vector<int> recv_counts(nrank, 0);
        std::vector<int> displs(nrank, 0);

        PARALLEL::mpi_gather(&send_count, 1, recv_counts.data(), 1, 0);
        PARALLEL::mpi_bcast(recv_counts.data(), nrank, 0);

        int total_recv = 0;
        for (int r = 0; r < nrank; ++r)
        {
        displs[r] = total_recv;
        total_recv += recv_counts[r];
        }

        std::vector<int> recv_buf(total_recv, 0);

        PARALLEL::mpi_allgatherv(
        send_count > 0 ? send_buf.data() : nullptr,
        send_count,
        total_recv > 0 ? recv_buf.data() : nullptr,
        recv_counts.data(),
        displs.data());

        if (total_recv % 17 != 0)
        {
        throw std::runtime_error(
            "build_topology: gathered edge-candidate int count is not multiple of 17.");
        }

        // ------------------------------------------------------------
        // 3) 在每个 rank 本地重建同一份全局 owner 选择结果
        // ------------------------------------------------------------
        std::unordered_map<EdgeKey, EntityKey, EdgeKey::Hash> global_owner;
        std::unordered_map<EdgeKey, std::vector<EntityKey>, EdgeKey::Hash> global_members;

        const int ncand = total_recv / 17;
        for (int c = 0; c < ncand; ++c)
        {
        const int *base = recv_buf.data() + 17 * c;

        EdgeKey key = unpack_edge_key(base + 0);
        EntityKey e = unpack_edge_local(base + 10);
        int8_t sign = static_cast<int8_t>(base[16]);

        equiv.edges.local_to_qkey[e] = key;
        equiv.edges.local_to_qsign[e] = sign;
        global_members[key].push_back(e);

        auto it = global_owner.find(key);
        if (it == global_owner.end() || e < it->second)
        {
            global_owner[key] = e;
        }
        }

        // ------------------------------------------------------------
        // 4) 只回填本 rank 本地涉及到的 key / local members
        // ------------------------------------------------------------
        equiv.edges.qkey_to_owner.clear();
        equiv.edges.local_is_owner.clear();
        equiv.edges.qkey_to_members.clear();

        for (auto &[key, members] : global_members)
        {
        std::sort(members.begin(), members.end());

        auto it = global_owner.find(key);
        if (it == global_owner.end())
        {
            throw std::runtime_error(
            "build_topology: shared edge key missing in global_owner.");
        }

        const EntityKey &owner = it->second;
        equiv.edges.qkey_to_owner[key] = owner;
        equiv.edges.qkey_to_members[key] = members;

        for (const auto &e : members)
        {
            equiv.edges.local_is_owner[e] = (e == owner);
        }
        }
    }

    
    void build_edge_owner_gid(int my_rank, Topology &equiv)
    {
        std::vector<EntityKey> local_owner_edges;
        local_owner_edges.reserve(equiv.edges.local_is_owner.size());

        for (const auto &[e, is_owner] : equiv.edges.local_is_owner)
        {
        if (is_owner && e.rank == my_rank)
            local_owner_edges.push_back(e);
        }

        std::sort(local_owner_edges.begin(), local_owner_edges.end());

        equiv.edges.n_local_owner = static_cast<int>(local_owner_edges.size());

        int nrank = 1;
        PARALLEL::mpi_size(&nrank);

        std::vector<int> counts(nrank, 0);
        PARALLEL::mpi_gather(&equiv.edges.n_local_owner, 1, counts.data(), 1, 0);
        PARALLEL::mpi_bcast(counts.data(), nrank, 0);

        equiv.edges.owner_gid_begin = 0;
        for (int r = 0; r < my_rank; ++r)
        {
        equiv.edges.owner_gid_begin += counts[r];
        }

        equiv.edges.owner_gid_end = equiv.edges.owner_gid_begin + equiv.edges.n_local_owner;

        equiv.edges.n_global_owner = 0;
        for (int r = 0; r < nrank; ++r)
        {
        equiv.edges.n_global_owner += counts[r];
        }

        equiv.edges.owner_to_gid.clear();
        equiv.edges.gid_to_owner.clear();

        std::vector<int> send_buf;
        send_buf.reserve(local_owner_edges.size() * 7);

        for (int n = 0; n < equiv.edges.n_local_owner; ++n)
        {
        const EntityKey &e = local_owner_edges[n];
        int gid = equiv.edges.owner_gid_begin + n;

        pack_edge_local(send_buf, e);
        send_buf.push_back(gid);
        }

        std::vector<int> recv_counts(nrank, 0);
        std::vector<int> displs(nrank, 0);
        int send_count = static_cast<int>(send_buf.size());

        PARALLEL::mpi_gather(&send_count, 1, recv_counts.data(), 1, 0);
        PARALLEL::mpi_bcast(recv_counts.data(), nrank, 0);

        int total_recv = 0;
        for (int r = 0; r < nrank; ++r)
        {
        displs[r] = total_recv;
        total_recv += recv_counts[r];
        }

        std::vector<int> recv_buf(total_recv, 0);
        PARALLEL::mpi_allgatherv(
        send_count > 0 ? send_buf.data() : nullptr,
        send_count,
        total_recv > 0 ? recv_buf.data() : nullptr,
        recv_counts.data(),
        displs.data());

        if (total_recv % 7 != 0)
        {
        throw std::runtime_error(
            "build_topology: gathered edge-owner-gid int count is not multiple of 7.");
        }

        const int nowners = total_recv / 7;
        for (int n = 0; n < nowners; ++n)
        {
        const int *base = recv_buf.data() + 7 * n;
        EntityKey e = unpack_edge_local(base);
        int gid = base[6];

        equiv.edges.owner_to_gid[e] = gid;
        equiv.edges.gid_to_owner[gid] = e;
        }
    }

    } // namespace TOPO::detail
