#include "4_halo/HaloEdgeOwner.h"
#include "0_basic/MPI_WRAPPER.h"

#include <algorithm>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace HALO_OWNER
{
    namespace
    {
        struct EdgeMemberCandidate
        {
            TOPO::EdgeKey key;
            TOPO::EdgeLocalID rep;
            int8_t sign_to_canonical; // ±1
        };

        inline void pack_node(std::vector<int> &buf, const TOPO::NodeEqID &x)
        {
            buf.push_back(x.rank);
            buf.push_back(x.gblock);
            buf.push_back(x.i);
            buf.push_back(x.j);
            buf.push_back(x.k);
        }

        inline TOPO::NodeEqID unpack_node_eq(const int *p)
        {
            return TOPO::NodeEqID{p[0], p[1], p[2], p[3], p[4]};
        }

        inline void pack_edge_key(std::vector<int> &buf, const TOPO::EdgeKey &key)
        {
            pack_node(buf, key.a);
            pack_node(buf, key.b);
        }

        inline TOPO::EdgeKey unpack_edge_key(const int *p)
        {
            return TOPO::EdgeKey{
                unpack_node_eq(p + 0),
                unpack_node_eq(p + 5)};
        }

        inline void pack_edge_local(std::vector<int> &buf, const TOPO::EdgeLocalID &e)
        {
            buf.push_back(e.rank);
            buf.push_back(e.gblock);
            buf.push_back(e.i);
            buf.push_back(e.j);
            buf.push_back(e.k);
            buf.push_back(e.dir);
        }

        inline TOPO::EdgeLocalID unpack_edge_local(const int *p)
        {
            return TOPO::EdgeLocalID{p[0], p[1], p[2], p[3], p[4], p[5]};
        }

        inline void pack_edge_member_candidate(std::vector<int> &buf,
                                               const EdgeMemberCandidate &c)
        {
            pack_edge_key(buf, c.key);
            pack_edge_local(buf, c.rep);
            buf.push_back((int)c.sign_to_canonical);
        }

        inline EdgeMemberCandidate unpack_edge_member_candidate(const int *p)
        {
            EdgeMemberCandidate c;
            c.key = unpack_edge_key(p + 0);
            c.rep = unpack_edge_local(p + 10);
            c.sign_to_canonical = static_cast<int8_t>(p[16]);
            return c;
        }

        inline int8_t factor_from_signs(int8_t s_rep, int8_t s_owner)
        {
            // factor = s_rep / s_owner = s_rep * s_owner, because s_owner = ±1
            return static_cast<int8_t>(s_rep * s_owner);
        }
    } // namespace

    void build_edge_owner_sync_pattern(
        const TOPO::TopologyEquiv &equiv,
        EdgeOwnerSyncPattern &pattern)
    {
        pattern.clear();

        int my_rank = 0;
        int nrank = 1;
        PARALLEL::mpi_rank(&my_rank);
        PARALLEL::mpi_size(&nrank);

        // ------------------------------------------------------------
        // 1) 打包本 rank 所有 local edge members
        // ------------------------------------------------------------
        std::vector<int> send_buf;
        send_buf.reserve(1024);

        for (const auto &[key, members] : equiv.edge_members)
        {
            for (const auto &rep : members)
            {
                auto it_sign = equiv.edge2sign.find(rep);
                if (it_sign == equiv.edge2sign.end())
                {
                    throw std::runtime_error(
                        "build_edge_owner_sync_pattern: local rep missing in equiv.edge2sign.");
                }

                EdgeMemberCandidate c;
                c.key = key;
                c.rep = rep;
                c.sign_to_canonical = it_sign->second;
                pack_edge_member_candidate(send_buf, c);
            }
        }

        // ------------------------------------------------------------
        // 2) gather counts + bcast counts + allgatherv data
        // ------------------------------------------------------------
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
                "build_edge_owner_sync_pattern: gathered member int count is not multiple of 17.");
        }

        // ------------------------------------------------------------
        // 3) 重建全局成员表：EdgeKey -> all global members
        // ------------------------------------------------------------
        std::unordered_map<TOPO::EdgeKey,
                           std::vector<EdgeMemberCandidate>,
                           TOPO::EdgeKey::Hash>
            global_members;

        const int ncand = total_recv / 17;
        for (int c = 0; c < ncand; ++c)
        {
            const int *base = recv_buf.data() + 17 * c;
            EdgeMemberCandidate cand = unpack_edge_member_candidate(base);
            global_members[cand.key].push_back(cand);
        }

        // ------------------------------------------------------------
        // 4) 针对“本 rank 涉及到的 key”，构造 local_alias / send / recv
        // ------------------------------------------------------------
        // ------------------------------------------------------------
        // 4) 针对“本 rank 涉及到的 key”，构造 local_alias / send / recv
        //    关键：不能只遍历 equiv.edge_members，
        //         因为 non-owner rank 可能没有 local edge_members[key]，
        //         但一定有 local edge2key[rep]。
        // ------------------------------------------------------------
        {
            std::unordered_set<TOPO::EdgeKey, TOPO::EdgeKey::Hash> local_keys;
            local_keys.reserve(equiv.edge2key.size());

            for (const auto &[rep, key] : equiv.edge2key)
                local_keys.insert(key);

            for (const auto &key : local_keys)
            {
                auto it_owner = equiv.edge_owner.find(key);
                if (it_owner == equiv.edge_owner.end())
                {
                    throw std::runtime_error(
                        "build_edge_owner_sync_pattern: local key missing in equiv.edge_owner.");
                }

                const TOPO::EdgeLocalID &owner = it_owner->second;

                auto it_global = global_members.find(key);
                if (it_global == global_members.end())
                {
                    throw std::runtime_error(
                        "build_edge_owner_sync_pattern: local key missing in global member table.");
                }

                const auto &members_all = it_global->second;

                // 找 owner 对应的 sign_to_canonical
                bool owner_sign_found = false;
                int8_t owner_sign = 0;

                for (const auto &m : members_all)
                {
                    if (m.rep == owner)
                    {
                        owner_sign = m.sign_to_canonical;
                        owner_sign_found = true;
                        break;
                    }
                }

                if (!owner_sign_found)
                {
                    throw std::runtime_error(
                        "build_edge_owner_sync_pattern: owner sign not found in global members.");
                }

                for (const auto &m : members_all)
                {
                    if (m.rep == owner)
                        continue;

                    const int8_t sign = factor_from_signs(m.sign_to_canonical, owner_sign);

                    if (owner.rank == my_rank)
                    {
                        if (m.rep.rank == my_rank)
                        {
                            // owner 和 rep 都在本 rank：本地 alias
                            pattern.local_alias.push_back(
                                EdgeOwnerLocalAliasItem{owner, m.rep, sign});
                        }
                        else
                        {
                            // 本 rank 是 owner，给远端 rep 所在 rank 发送
                            pattern.send_items.push_back(
                                EdgeOwnerSendItem{m.rep.rank, key, owner});
                        }
                    }
                    else
                    {
                        if (m.rep.rank == my_rank)
                        {
                            // owner 在远端，本 rank 持有本地 non-owner rep：需要 recv
                            pattern.recv_items.push_back(
                                EdgeOwnerRecvItem{owner.rank, key, m.rep, sign});
                        }
                    }
                }
            }
        }

        // ------------------------------------------------------------
        // 5) 为 send/recv items 按 tar_id 排序，并生成 counts/displs
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
        // 6) runtime caches: allocate once
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