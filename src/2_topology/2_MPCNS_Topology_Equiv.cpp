#include "2_topology/2_MPCNS_Topology_Equiv.h"
#include "2_topology/TopologyOps.h"

#include "0_basic/MPI_WRAPPER.h"
#include "0_basic/BoxOps.h"

#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace TOPO
{
    namespace
    {
        // ============================================================
        // small internal helpers
        // ============================================================

        struct UnionFind
        {
            std::vector<int> parent;
            std::vector<int> rank;

            int add()
            {
                int id = static_cast<int>(parent.size());
                parent.push_back(id);
                rank.push_back(0);
                return id;
            }

            int find(int x)
            {
                if (parent[x] != x)
                    parent[x] = find(parent[x]);
                return parent[x];
            }

            void unite(int a, int b)
            {
                a = find(a);
                b = find(b);
                if (a == b)
                    return;

                if (rank[a] < rank[b])
                    std::swap(a, b);
                parent[b] = a;
                if (rank[a] == rank[b])
                    ++rank[a];
            }
        };

        inline void pack_node(std::vector<int> &buf, const LocalNodeID &x)
        {
            buf.push_back(x.rank);
            buf.push_back(x.gblock);
            buf.push_back(x.i);
            buf.push_back(x.j);
            buf.push_back(x.k);
        }

        inline LocalNodeID unpack_node(const int *p)
        {
            return LocalNodeID{p[0], p[1], p[2], p[3], p[4]};
        }

        inline void pack_pair(std::vector<int> &buf,
                              const LocalNodeID &a,
                              const LocalNodeID &b)
        {
            pack_node(buf, a);
            pack_node(buf, b);
        }

        // ============================================================
        // local collection
        // ============================================================

        std::vector<LocalNodeID> collect_all_local_nodes_impl(
            Grid &grid,
            int my_rank)
        {
            std::vector<LocalNodeID> nodes;

            for (int ib = 0; ib < grid.nblock; ++ib)
            {
                const auto &blk = grid.grids(ib);

                for (int i = 0; i <= blk.mx; ++i)
                    for (int j = 0; j <= blk.my; ++j)
                        for (int k = 0; k <= blk.mz; ++k)
                        {
                            nodes.push_back(LocalNodeID{my_rank, ib, i, j, k});
                        }
            }
            return nodes;
        }

        std::vector<EdgeLocalID> collect_all_local_edges_impl(
            Grid &grid,
            int my_rank,
            int dimension)
        {
            std::vector<EdgeLocalID> edges;

            for (int ib = 0; ib < grid.nblock; ++ib)
            {
                const auto &blk = grid.grids(ib);

                // Xi
                for (int i = 0; i < blk.mx; ++i)
                    for (int j = 0; j <= blk.my; ++j)
                        for (int k = 0; k <= blk.mz; ++k)
                        {
                            edges.push_back(EdgeLocalID{my_rank, ib, i, j, k, 1});
                        }

                if (dimension >= 2)
                {
                    // Eta
                    for (int i = 0; i <= blk.mx; ++i)
                        for (int j = 0; j < blk.my; ++j)
                            for (int k = 0; k <= blk.mz; ++k)
                            {
                                edges.push_back(EdgeLocalID{my_rank, ib, i, j, k, 2});
                            }
                }

                if (dimension >= 3)
                {
                    // Zeta
                    for (int i = 0; i <= blk.mx; ++i)
                        for (int j = 0; j <= blk.my; ++j)
                            for (int k = 0; k < blk.mz; ++k)
                            {
                                edges.push_back(EdgeLocalID{my_rank, ib, i, j, k, 3});
                            }
                }
            }
            return edges;
        }

        // ============================================================
        // internal build phases
        // ============================================================

        void reconcile_node_equivalence_parallel_impl(
            const Topology &topo,
            int my_rank,
            const std::vector<LocalNodeID> &all_local_nodes,
            TopologyEquiv &equiv)
        {
            std::vector<int> send_buf;
            send_buf.reserve(1024);

            auto collect_patch_pairs = [&](const InterfacePatch &patch)
            {
                for (int i = patch.this_box_node.lo.i; i < patch.this_box_node.hi.i; ++i)
                    for (int j = patch.this_box_node.lo.j; j < patch.this_box_node.hi.j; ++j)
                        for (int k = patch.this_box_node.lo.k; k < patch.this_box_node.hi.k; ++k)
                        {
                            Int3 p_this{i, j, k};
                            Int3 p_nb = map_node_point(p_this, patch.trans);

                            if (!BOX::contains_point(patch.nb_box_node, p_nb))
                            {
                                std::ostringstream oss;
                                oss << "build_topology_equiv: mapped node out of nb_box_node. "
                                    << "this_rank=" << patch.this_rank
                                    << " nb_rank=" << patch.nb_rank
                                    << " this_block=" << patch.this_block
                                    << " nb_block=" << patch.nb_block
                                    << " p_this=(" << p_this.i << "," << p_this.j << "," << p_this.k << ")"
                                    << " p_nb=(" << p_nb.i << "," << p_nb.j << "," << p_nb.k << ")";
                                throw std::runtime_error(oss.str());
                            }

                            LocalNodeID a{patch.this_rank, patch.this_block, i, j, k};
                            LocalNodeID b{patch.nb_rank, patch.nb_block, p_nb.i, p_nb.j, p_nb.k};
                            pack_pair(send_buf, a, b);
                        }
            };

            for (const auto &patch : topo.inner_patches)
                collect_patch_pairs(patch);
            for (const auto &patch : topo.parallel_patches)
                collect_patch_pairs(patch);

            int nrank = 1;
            PARALLEL::mpi_size(&nrank);

            int send_count = static_cast<int>(send_buf.size());

            std::vector<int> recv_counts(nrank, 0);
            std::vector<int> displs(nrank, 0);

            // 先 gather 到 root=0
            PARALLEL::mpi_gather(&send_count, 1, recv_counts.data(), 1, 0);

            // 再广播给所有 rank
            PARALLEL::mpi_bcast(recv_counts.data(), nrank, 0);

            int total_recv = 0;
            for (int r = 0; r < nrank; ++r)
            {
                displs[r] = total_recv;
                total_recv += recv_counts[r];
            }

            std::vector<int> recv_buf(total_recv, 0);

            // allgatherv 真正收所有 pair 数据
            PARALLEL::mpi_allgatherv(
                send_count > 0 ? send_buf.data() : nullptr,
                send_count,
                total_recv > 0 ? recv_buf.data() : nullptr,
                recv_counts.data(),
                displs.data());

            if (total_recv % 10 != 0)
            {
                throw std::runtime_error("build_topology_equiv: gathered node-pair int count is not multiple of 10.");
            }

            std::unordered_map<LocalNodeID, int, LocalNodeID::Hash> node2idx;
            std::vector<LocalNodeID> idx2node;
            UnionFind uf;

            auto get_or_add_node = [&](const LocalNodeID &nid) -> int
            {
                auto it = node2idx.find(nid);
                if (it != node2idx.end())
                    return it->second;

                int idx = uf.add();
                node2idx.emplace(nid, idx);
                idx2node.push_back(nid);
                return idx;
            };

            for (const auto &nid : all_local_nodes)
            {
                get_or_add_node(nid);
            }

            const int npair = total_recv / 10;
            for (int p = 0; p < npair; ++p)
            {
                const int *base = recv_buf.data() + 10 * p;
                LocalNodeID a = unpack_node(base + 0);
                LocalNodeID b = unpack_node(base + 5);

                int ia = get_or_add_node(a);
                int ib = get_or_add_node(b);
                uf.unite(ia, ib);
            }

            std::unordered_map<int, LocalNodeID> root_min_node;
            for (int idx = 0; idx < static_cast<int>(idx2node.size()); ++idx)
            {
                int root = uf.find(idx);
                const LocalNodeID &nid = idx2node[idx];

                auto it = root_min_node.find(root);
                if (it == root_min_node.end() || nid < it->second)
                {
                    root_min_node[root] = nid;
                }
            }

            equiv.node2eq.clear();
            for (const auto &nid : all_local_nodes)
            {
                auto it = node2idx.find(nid);
                if (it == node2idx.end())
                {
                    std::ostringstream oss;
                    oss << "build_topology_equiv: local node missing in node2idx. "
                        << "(" << nid.rank << "," << nid.gblock << ","
                        << nid.i << "," << nid.j << "," << nid.k << ")";
                    throw std::runtime_error(oss.str());
                }

                int root = uf.find(it->second);
                const LocalNodeID &rep = root_min_node.at(root);
                equiv.node2eq[nid] = to_node_eq_id(rep);
            }
        }

        void build_edge_equivalence_from_nodes_impl(
            const std::vector<EdgeLocalID> &all_local_edges,
            TopologyEquiv &equiv)
        {
            equiv.edge2key.clear();
            equiv.edge2sign.clear();
            equiv.edge_members.clear();

            for (const auto &e : all_local_edges)
            {
                int8_t sign_to_canonical = 0;
                EdgeKey key = make_edge_key(e, equiv.node2eq, sign_to_canonical);

                equiv.edge2key[e] = key;
                equiv.edge2sign[e] = sign_to_canonical;
                equiv.edge_members[key].push_back(e);
            }
        }

        inline void pack_edge_local(std::vector<int> &buf, const EdgeLocalID &e)
        {
            buf.push_back(e.rank);
            buf.push_back(e.gblock);
            buf.push_back(e.i);
            buf.push_back(e.j);
            buf.push_back(e.k);
            buf.push_back(e.dir);
        }

        inline EdgeLocalID unpack_edge_local(const int *p)
        {
            return EdgeLocalID{p[0], p[1], p[2], p[3], p[4], p[5]};
        }

        inline void pack_edge_key(std::vector<int> &buf, const EdgeKey &key)
        {
            pack_node(buf, to_local_node_id(key.a));
            pack_node(buf, to_local_node_id(key.b));
        }

        inline EdgeKey unpack_edge_key(const int *p)
        {
            LocalNodeID a = unpack_node(p + 0);
            LocalNodeID b = unpack_node(p + 5);
            return EdgeKey{to_node_eq_id(a), to_node_eq_id(b)};
        }

        inline void pack_edge_owner_candidate(
            std::vector<int> &buf,
            const EdgeKey &key,
            const EdgeLocalID &e)
        {
            // 10 ints for EdgeKey + 6 ints for EdgeLocalID = 16 ints
            pack_edge_key(buf, key);
            pack_edge_local(buf, e);
        }

        inline void select_edge_owner_parallel_impl(TopologyEquiv &equiv)
        {
            // ------------------------------------------------------------
            // 1) 打包本 rank 的所有 local edge candidates
            //    每个 local member 都作为一个候选发出去
            // ------------------------------------------------------------
            std::vector<int> send_buf;
            send_buf.reserve(1024);

            for (const auto &[key, members] : equiv.edge_members)
            {
                for (const auto &e : members)
                {
                    pack_edge_owner_candidate(send_buf, key, e);
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

            if (total_recv % 16 != 0)
            {
                throw std::runtime_error(
                    "build_topology_equiv: gathered edge-candidate int count is not multiple of 16.");
            }

            // ------------------------------------------------------------
            // 3) 在每个 rank 本地重建同一份全局 owner 选择结果
            // ------------------------------------------------------------
            std::unordered_map<EdgeKey, EdgeLocalID, EdgeKey::Hash> global_owner;

            const int ncand = total_recv / 16;
            for (int c = 0; c < ncand; ++c)
            {
                const int *base = recv_buf.data() + 16 * c;

                EdgeKey key = unpack_edge_key(base + 0);
                EdgeLocalID e = unpack_edge_local(base + 10);

                auto it = global_owner.find(key);
                if (it == global_owner.end() || e < it->second)
                {
                    global_owner[key] = e;
                }
            }

            // ------------------------------------------------------------
            // 4) 只回填本 rank 本地涉及到的 key / local members
            // ------------------------------------------------------------
            equiv.edge_owner.clear();
            equiv.edge_is_owner.clear();

            for (const auto &[key, members] : equiv.edge_members)
            {
                auto it = global_owner.find(key);
                if (it == global_owner.end())
                {
                    throw std::runtime_error(
                        "build_topology_equiv: local edge key missing in global_owner.");
                }

                const EdgeLocalID &owner = it->second;
                equiv.edge_owner[key] = owner;

                for (const auto &e : members)
                {
                    equiv.edge_is_owner[e] = (e == owner);
                }
            }
        }

        inline void build_edge_owner_gid_impl(int my_rank, TopologyEquiv &equiv)
        {
            std::vector<EdgeLocalID> local_owner_edges;
            local_owner_edges.reserve(equiv.edge_is_owner.size());

            for (const auto &[e, is_owner] : equiv.edge_is_owner)
            {
                if (is_owner)
                    local_owner_edges.push_back(e);
            }

            std::sort(local_owner_edges.begin(), local_owner_edges.end());

            equiv.n_local_edge_owner = static_cast<int>(local_owner_edges.size());

            int nrank = 1;
            PARALLEL::mpi_size(&nrank);

            std::vector<int> counts(nrank, 0);
            PARALLEL::mpi_gather(&equiv.n_local_edge_owner, 1, counts.data(), 1, 0);
            PARALLEL::mpi_bcast(counts.data(), nrank, 0);

            equiv.edge_owner_gid_begin = 0;
            for (int r = 0; r < my_rank; ++r)
            {
                equiv.edge_owner_gid_begin += counts[r];
            }

            equiv.edge_owner_gid_end = equiv.edge_owner_gid_begin + equiv.n_local_edge_owner;

            equiv.n_global_edge_owner = 0;
            for (int r = 0; r < nrank; ++r)
            {
                equiv.n_global_edge_owner += counts[r];
            }

            equiv.edge_owner_gid.clear();
            equiv.gid2edge_owner.clear();

            for (int n = 0; n < equiv.n_local_edge_owner; ++n)
            {
                const EdgeLocalID &e = local_owner_edges[n];
                int gid = equiv.edge_owner_gid_begin + n;

                equiv.edge_owner_gid[e] = gid;
                equiv.gid2edge_owner[gid] = e;
            }
        }
    } // anonymous namespace

    std::pair<LocalNodeID, LocalNodeID> endpoints(const EdgeLocalID &e)
    {
        LocalNodeID n0{e.rank, e.gblock, e.i, e.j, e.k};
        LocalNodeID n1 = n0;

        switch (e.dir)
        {
        case 1:
            ++n1.i;
            break; // Xi
        case 2:
            ++n1.j;
            break; // Eta
        case 3:
            ++n1.k;
            break; // Zeta
        default:
        {
            std::ostringstream oss;
            oss << "TOPO::endpoints: invalid edge dir = " << e.dir;
            throw std::runtime_error(oss.str());
        }
        }
        return {n0, n1};
    }

    EdgeKey make_edge_key(
        const EdgeLocalID &e,
        const std::unordered_map<LocalNodeID, NodeEqID, LocalNodeID::Hash> &node2eq,
        int8_t &sign_to_canonical)
    {
        const auto [ln0, ln1] = endpoints(e);

        auto it0 = node2eq.find(ln0);
        auto it1 = node2eq.find(ln1);

        if (it0 == node2eq.end() || it1 == node2eq.end())
        {
            throw std::runtime_error("TOPO::make_edge_key: endpoint not found in node2eq.");
        }

        const NodeEqID &g0 = it0->second;
        const NodeEqID &g1 = it1->second;

        if (g0 == g1)
        {
            throw std::runtime_error("TOPO::make_edge_key: two endpoints collapse to the same NodeEqID.");
        }

        if (g0 < g1)
        {
            sign_to_canonical = +1;
            return EdgeKey{g0, g1};
        }
        else
        {
            sign_to_canonical = -1;
            return EdgeKey{g1, g0};
        }
    }

    void build_topology_equiv(
        const Topology &topo,
        Grid &grid,
        int my_rank,
        int dimension,
        TopologyEquiv &equiv)
    {
        equiv.clear();

        auto all_local_nodes = collect_all_local_nodes_impl(grid, my_rank);
        reconcile_node_equivalence_parallel_impl(topo, my_rank, all_local_nodes, equiv);

        auto all_local_edges = collect_all_local_edges_impl(grid, my_rank, dimension);
        build_edge_equivalence_from_nodes_impl(all_local_edges, equiv);

        select_edge_owner_parallel_impl(equiv);
        build_edge_owner_gid_impl(my_rank, equiv);
    }

} // namespace TOPO