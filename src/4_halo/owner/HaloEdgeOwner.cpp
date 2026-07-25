#include "4_halo/HaloEdgeOwner.h"
#include "0_basic/MPI_WRAPPER.h"

#include <algorithm>
#include <stdexcept>
#include <vector>

namespace HALO_OWNER
{
    namespace
    {
        inline int8_t factor_from_signs(int8_t s_rep, int8_t s_owner)
        {
            // factor = s_rep / s_owner = s_rep * s_owner, because s_owner = ±1
            return static_cast<int8_t>(s_rep * s_owner);
        }
    } // namespace

    void build_edge_owner_sync_pattern(
        const TOPO::Topology &equiv,
        EdgeOwnerSyncPattern &pattern)
    {
        pattern.clear();

        int my_rank = 0;
        int nrank = 1;
        PARALLEL::mpi_rank(&my_rank);
        PARALLEL::mpi_size(&nrank);

        // ------------------------------------------------------------
        // 1) 直接使用 topology 已经收集好的全局共享边成员表。
        //
        // select_edge_owner_parallel() 已在每个 rank 上生成相同的
        // qkey_to_members/qkey_to_owner/local_to_qsign。这里若再次
        // Allgatherv，会把每个成员重复 nrank 次，令内存及同步条目
        // 都产生 O(nrank) 的额外放大。
        // ------------------------------------------------------------
        for (const auto &[key, members_all] : equiv.edges.qkey_to_members)
        {
            auto it_owner = equiv.edges.qkey_to_owner.find(key);
            if (it_owner == equiv.edges.qkey_to_owner.end())
            {
                throw std::runtime_error(
                    "build_edge_owner_sync_pattern: shared edge key missing owner.");
            }

            const TOPO::EntityKey &owner = it_owner->second;
            if (owner.rank < 0 || owner.rank >= nrank)
            {
                throw std::runtime_error(
                    "build_edge_owner_sync_pattern: owner rank is out of range.");
            }

            auto it_owner_sign = equiv.edges.local_to_qsign.find(owner);
            if (it_owner_sign == equiv.edges.local_to_qsign.end())
            {
                throw std::runtime_error(
                    "build_edge_owner_sync_pattern: owner sign missing in topology.");
            }
            const int8_t owner_sign = it_owner_sign->second;

            for (const auto &rep : members_all)
            {
                if (rep.rank < 0 || rep.rank >= nrank)
                {
                    throw std::runtime_error(
                        "build_edge_owner_sync_pattern: member rank is out of range.");
                }

                if (rep == owner)
                    continue;

                auto it_rep_sign = equiv.edges.local_to_qsign.find(rep);
                if (it_rep_sign == equiv.edges.local_to_qsign.end())
                {
                    throw std::runtime_error(
                        "build_edge_owner_sync_pattern: member sign missing in topology.");
                }
                const int8_t sign = factor_from_signs(it_rep_sign->second, owner_sign);

                if (owner.rank == my_rank)
                {
                    if (rep.rank == my_rank)
                    {
                        pattern.local_alias.push_back(
                            EdgeOwnerLocalAliasItem{owner, rep, sign});
                    }
                    else
                    {
                        pattern.send_items.push_back(
                            EdgeOwnerSendItem{rep.rank, key, owner});
                    }
                }
                else if (rep.rank == my_rank)
                {
                    pattern.recv_items.push_back(
                        EdgeOwnerRecvItem{owner.rank, key, rep, sign});
                }
            }
        }

        // ------------------------------------------------------------
        // 2) 为 send/recv items 按 tar_id 排序，并生成 counts/displs
        // ------------------------------------------------------------
        std::sort(pattern.send_items.begin(), pattern.send_items.end(),
                  [](const EdgeOwnerSendItem &a, const EdgeOwnerSendItem &b)
                  {
                      if (a.tar_id != b.tar_id)
                          return a.tar_id < b.tar_id;
                      return a.key < b.key;
                  });

        std::sort(pattern.recv_items.begin(), pattern.recv_items.end(),
                  [](const EdgeOwnerRecvItem &a, const EdgeOwnerRecvItem &b)
                  {
                      if (a.tar_id != b.tar_id)
                          return a.tar_id < b.tar_id;
                      return a.key < b.key;
                  });

        pattern.send_counts.assign(nrank, 0);
        pattern.recv_counts.assign(nrank, 0);
        pattern.send_displs.assign(nrank, 0);
        pattern.recv_displs.assign(nrank, 0);

        for (const auto &it : pattern.send_items)
            ++pattern.send_counts[it.tar_id];
        for (const auto &it : pattern.recv_items)
            ++pattern.recv_counts[it.tar_id];

        for (int r = 1; r < nrank; ++r)
        {
            pattern.send_displs[r] = pattern.send_displs[r - 1] + pattern.send_counts[r - 1];
            pattern.recv_displs[r] = pattern.recv_displs[r - 1] + pattern.recv_counts[r - 1];
        }

        // ------------------------------------------------------------
        // 3) runtime caches: allocate once
        // ------------------------------------------------------------
        pattern.send_buf_cache.clear();
        pattern.recv_buf_cache.clear();
        pattern.send_buf_cache.resize(nrank);
        pattern.recv_buf_cache.resize(nrank);

        int nrecv_active = 0;
        int nsend_active = 0;

        for (int r = 0; r < nrank; ++r)
        {
            const int send_len = pattern.send_counts[r];
            const int recv_len = pattern.recv_counts[r];

            if (send_len > 0)
            {
                pattern.send_buf_cache[r].resize(
                    static_cast<std::size_t>(send_len) * 1u); // ncomp unknown here; later resize if needed
                ++nsend_active;
            }

            if (recv_len > 0)
            {
                pattern.recv_buf_cache[r].resize(
                    static_cast<std::size_t>(recv_len) * 1u); // ncomp unknown here; later resize if needed
                ++nrecv_active;
            }
        }

        pattern.req_recv_cache.resize(nrecv_active);
        pattern.stat_recv_cache.resize(nrecv_active);
        pattern.req_send_cache.resize(nsend_active);
        pattern.stat_send_cache.resize(nsend_active);
    }

} // namespace HALO_OWNER
