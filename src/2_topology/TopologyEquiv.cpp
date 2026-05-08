#include "2_topology/TopologyEquiv.h"
#include "2_topology/TopologyOps.h"

#include "0_basic/MPI_WRAPPER.h"
#include "0_basic/BoxOps.h"
#include "0_basic/Error.h"
#include "1_grid/1_MPCNS_Grid.h"

#include <algorithm>
#include <array>
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

        std::vector<FaceLocalID> collect_all_local_faces_impl(
            Grid &grid,
            int my_rank,
            int dimension)
        {
            std::vector<FaceLocalID> faces;

            for (int ib = 0; ib < grid.nblock; ++ib)
            {
                const auto &blk = grid.grids(ib);

                if (dimension >= 3)
                {
                    // FaceXi
                    for (int i = 0; i <= blk.mx; ++i)
                        for (int j = 0; j < blk.my; ++j)
                            for (int k = 0; k < blk.mz; ++k)
                            {
                                faces.push_back(FaceLocalID{my_rank, ib, i, j, k, 1});
                            }

                    // FaceEt
                    for (int i = 0; i < blk.mx; ++i)
                        for (int j = 0; j <= blk.my; ++j)
                            for (int k = 0; k < blk.mz; ++k)
                            {
                                faces.push_back(FaceLocalID{my_rank, ib, i, j, k, 2});
                            }

                    // FaceZe
                    for (int i = 0; i < blk.mx; ++i)
                        for (int j = 0; j < blk.my; ++j)
                            for (int k = 0; k <= blk.mz; ++k)
                            {
                                faces.push_back(FaceLocalID{my_rank, ib, i, j, k, 3});
                            }
                }
                else if (dimension >= 2)
                {
                    // 2D uses xi/eta face locations as line faces in the single k plane.
                    const int k = 0;

                    for (int i = 0; i <= blk.mx; ++i)
                        for (int j = 0; j < blk.my; ++j)
                        {
                            faces.push_back(FaceLocalID{my_rank, ib, i, j, k, 1});
                        }

                    for (int i = 0; i < blk.mx; ++i)
                        for (int j = 0; j <= blk.my; ++j)
                        {
                            faces.push_back(FaceLocalID{my_rank, ib, i, j, k, 2});
                        }
                }
            }
            return faces;
        }

        // ============================================================
        // internal build phases
        // ============================================================

        void reconcile_node_equivalence_parallel_impl(
            const Topology &topo,
            int my_rank,
            const std::vector<LocalNodeID> &all_local_nodes,
            TopologyEquiv &equiv,
            std::unordered_map<NodeEqID, int, NodeEqID::Hash> &node_eq_counts)
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
            std::unordered_map<int, int> root_node_count;
            for (int idx = 0; idx < static_cast<int>(idx2node.size()); ++idx)
            {
                int root = uf.find(idx);
                const LocalNodeID &nid = idx2node[idx];

                ++root_node_count[root];

                auto it = root_min_node.find(root);
                if (it == root_min_node.end() || nid < it->second)
                {
                    root_min_node[root] = nid;
                }
            }

            equiv.node2eq.clear();
            node_eq_counts.clear();
            for (const auto &[root, rep] : root_min_node)
            {
                node_eq_counts[to_node_eq_id(rep)] = root_node_count[root];
            }

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

        inline int node_equiv_count(
            const std::unordered_map<NodeEqID, int, NodeEqID::Hash> &node_eq_counts,
            const NodeEqID &node)
        {
            auto it = node_eq_counts.find(node);
            return (it != node_eq_counts.end()) ? it->second : 1;
        }

        inline bool edge_key_is_shared(
            const EdgeKey &key,
            const std::unordered_map<NodeEqID, int, NodeEqID::Hash> &node_eq_counts)
        {
            return node_equiv_count(node_eq_counts, key.a) > 1 ||
                   node_equiv_count(node_eq_counts, key.b) > 1;
        }

        inline bool face_key_is_shared(
            const FaceKey &key,
            const std::unordered_map<NodeEqID, int, NodeEqID::Hash> &node_eq_counts)
        {
            return node_equiv_count(node_eq_counts, key.a) > 1 ||
                   node_equiv_count(node_eq_counts, key.b) > 1 ||
                   node_equiv_count(node_eq_counts, key.c) > 1 ||
                   node_equiv_count(node_eq_counts, key.d) > 1;
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

        inline void select_edge_owner_parallel_impl(
            TopologyEquiv &equiv,
            const std::unordered_map<NodeEqID, int, NodeEqID::Hash> &node_eq_counts)
        {
            // ------------------------------------------------------------
            // 1) 打包本 rank 的所有 local edge candidates
            //    每个 local member 都作为一个候选发出去
            // ------------------------------------------------------------
            std::vector<int> send_buf;
            send_buf.reserve(1024);

            for (const auto &[key, members] : equiv.edge_members)
            {
                if (!edge_key_is_shared(key, node_eq_counts) && members.size() <= 1)
                    continue;

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
            std::unordered_map<EdgeKey, std::vector<EdgeLocalID>, EdgeKey::Hash> global_members;

            const int ncand = total_recv / 16;
            for (int c = 0; c < ncand; ++c)
            {
                const int *base = recv_buf.data() + 16 * c;

                EdgeKey key = unpack_edge_key(base + 0);
                EdgeLocalID e = unpack_edge_local(base + 10);

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
            equiv.edge_owner.clear();
            equiv.edge_is_owner.clear();
            equiv.edge_members.clear();

            for (auto &[key, members] : global_members)
            {
                std::sort(members.begin(), members.end());

                auto it = global_owner.find(key);
                if (it == global_owner.end())
                {
                    throw std::runtime_error(
                        "build_topology_equiv: shared edge key missing in global_owner.");
                }

                const EdgeLocalID &owner = it->second;
                equiv.edge_owner[key] = owner;
                equiv.edge_members[key] = members;

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
                if (is_owner && e.rank == my_rank)
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

            std::vector<int> send_buf;
            send_buf.reserve(local_owner_edges.size() * 7);

            for (int n = 0; n < equiv.n_local_edge_owner; ++n)
            {
                const EdgeLocalID &e = local_owner_edges[n];
                int gid = equiv.edge_owner_gid_begin + n;

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
                    "build_topology_equiv: gathered edge-owner-gid int count is not multiple of 7.");
            }

            const int nowners = total_recv / 7;
            for (int n = 0; n < nowners; ++n)
            {
                const int *base = recv_buf.data() + 7 * n;
                EdgeLocalID e = unpack_edge_local(base);
                int gid = base[6];

                equiv.edge_owner_gid[e] = gid;
                equiv.gid2edge_owner[gid] = e;
            }
        }

        void build_face_equivalence_from_nodes_impl(
            const std::vector<FaceLocalID> &all_local_faces,
            TopologyEquiv &equiv)
        {
            equiv.face2key.clear();
            equiv.face2sign.clear();
            equiv.face_members.clear();

            for (const auto &f : all_local_faces)
            {
                int8_t sign_to_canonical = 0;
                FaceKey key = make_face_key(f, equiv.node2eq, sign_to_canonical);

                equiv.face2key[f] = key;
                equiv.face2sign[f] = sign_to_canonical;
                equiv.face_members[key].push_back(f);
            }
        }

        inline void pack_face_local(std::vector<int> &buf, const FaceLocalID &f)
        {
            buf.push_back(f.rank);
            buf.push_back(f.gblock);
            buf.push_back(f.i);
            buf.push_back(f.j);
            buf.push_back(f.k);
            buf.push_back(f.dir);
        }

        inline FaceLocalID unpack_face_local(const int *p)
        {
            return FaceLocalID{p[0], p[1], p[2], p[3], p[4], p[5]};
        }

        inline void pack_face_key(std::vector<int> &buf, const FaceKey &key)
        {
            pack_node(buf, to_local_node_id(key.a));
            pack_node(buf, to_local_node_id(key.b));
            pack_node(buf, to_local_node_id(key.c));
            pack_node(buf, to_local_node_id(key.d));
        }

        inline FaceKey unpack_face_key(const int *p)
        {
            LocalNodeID a = unpack_node(p + 0);
            LocalNodeID b = unpack_node(p + 5);
            LocalNodeID c = unpack_node(p + 10);
            LocalNodeID d = unpack_node(p + 15);
            return FaceKey{to_node_eq_id(a), to_node_eq_id(b), to_node_eq_id(c), to_node_eq_id(d)};
        }

        inline void pack_face_owner_candidate(
            std::vector<int> &buf,
            const FaceKey &key,
            const FaceLocalID &f,
            int8_t sign)
        {
            // 20 ints for FaceKey + 6 ints for FaceLocalID + 1 sign = 27 ints
            pack_face_key(buf, key);
            pack_face_local(buf, f);
            buf.push_back(static_cast<int>(sign));
        }

        inline void select_face_owner_parallel_impl(
            TopologyEquiv &equiv,
            const std::unordered_map<NodeEqID, int, NodeEqID::Hash> &node_eq_counts)
        {
            std::vector<int> send_buf;
            send_buf.reserve(1024);

            for (const auto &[key, members] : equiv.face_members)
            {
                if (!face_key_is_shared(key, node_eq_counts) && members.size() <= 1)
                    continue;

                for (const auto &f : members)
                {
                    auto sign_it = equiv.face2sign.find(f);
                    const int8_t sign = (sign_it != equiv.face2sign.end()) ? sign_it->second : int8_t{+1};
                    pack_face_owner_candidate(send_buf, key, f, sign);
                }
            }

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

            if (total_recv % 27 != 0)
            {
                throw std::runtime_error(
                    "build_topology_equiv: gathered face-candidate int count is not multiple of 27.");
            }

            std::unordered_map<FaceKey, FaceLocalID, FaceKey::Hash> global_owner;
            std::unordered_map<FaceKey, std::vector<FaceLocalID>, FaceKey::Hash> global_members;

            const int ncand = total_recv / 27;
            for (int c = 0; c < ncand; ++c)
            {
                const int *base = recv_buf.data() + 27 * c;

                FaceKey key = unpack_face_key(base + 0);
                FaceLocalID f = unpack_face_local(base + 20);
                int8_t sign = static_cast<int8_t>(base[26]);

                equiv.face2sign[f] = sign;
                equiv.face2key[f] = key;
                global_members[key].push_back(f);

                auto it = global_owner.find(key);
                if (it == global_owner.end() || f < it->second)
                {
                    global_owner[key] = f;
                }
            }

            equiv.face_owner.clear();
            equiv.face_is_owner.clear();
            equiv.face_members.clear();

            for (auto &[key, members] : global_members)
            {
                std::sort(members.begin(), members.end());

                auto it = global_owner.find(key);
                if (it == global_owner.end())
                {
                    throw std::runtime_error(
                        "build_topology_equiv: local face key missing in global_owner.");
                }

                const FaceLocalID &owner = it->second;
                equiv.face_owner[key] = owner;
                equiv.face_members[key] = members;

                for (const auto &f : members)
                {
                    equiv.face_is_owner[f] = (f == owner);
                }
            }

        }

        inline void build_face_owner_gid_impl(int my_rank, TopologyEquiv &equiv)
        {
            std::vector<FaceLocalID> local_owner_faces;
            local_owner_faces.reserve(equiv.face_is_owner.size());

            for (const auto &[f, is_owner] : equiv.face_is_owner)
            {
                if (is_owner && f.rank == my_rank)
                    local_owner_faces.push_back(f);
            }

            std::sort(local_owner_faces.begin(), local_owner_faces.end());

            equiv.n_local_face_owner = static_cast<int>(local_owner_faces.size());

            int nrank = 1;
            PARALLEL::mpi_size(&nrank);

            std::vector<int> counts(nrank, 0);
            PARALLEL::mpi_gather(&equiv.n_local_face_owner, 1, counts.data(), 1, 0);
            PARALLEL::mpi_bcast(counts.data(), nrank, 0);

            equiv.face_owner_gid_begin = 0;
            for (int r = 0; r < my_rank; ++r)
            {
                equiv.face_owner_gid_begin += counts[r];
            }

            equiv.face_owner_gid_end = equiv.face_owner_gid_begin + equiv.n_local_face_owner;

            equiv.n_global_face_owner = 0;
            for (int r = 0; r < nrank; ++r)
            {
                equiv.n_global_face_owner += counts[r];
            }

            equiv.face_owner_gid.clear();
            equiv.gid2face_owner.clear();

            std::vector<int> send_buf;
            send_buf.reserve(local_owner_faces.size() * 7);

            for (int n = 0; n < equiv.n_local_face_owner; ++n)
            {
                const FaceLocalID &f = local_owner_faces[n];
                int gid = equiv.face_owner_gid_begin + n;

                pack_face_local(send_buf, f);
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
                    "build_topology_equiv: gathered face-owner-gid int count is not multiple of 7.");
            }

            const int nowners = total_recv / 7;
            for (int n = 0; n < nowners; ++n)
            {
                const int *base = recv_buf.data() + 7 * n;
                FaceLocalID f = unpack_face_local(base);
                int gid = base[6];

                equiv.face_owner_gid[f] = gid;
                equiv.gid2face_owner[gid] = f;
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

    std::array<LocalNodeID, 4> corners(const FaceLocalID &f)
    {
        LocalNodeID n0{f.rank, f.gblock, f.i, f.j, f.k};
        LocalNodeID n1 = n0;
        LocalNodeID n2 = n0;
        LocalNodeID n3 = n0;

        switch (f.dir)
        {
        case 1:
            ++n1.j;
            ++n2.k;
            ++n3.j;
            ++n3.k;
            break; // FaceXi
        case 2:
            ++n1.i;
            ++n2.k;
            ++n3.i;
            ++n3.k;
            break; // FaceEt
        case 3:
            ++n1.i;
            ++n2.j;
            ++n3.i;
            ++n3.j;
            break; // FaceZe
        default:
        {
            std::ostringstream oss;
            oss << "TOPO::corners: invalid face dir = " << f.dir;
            throw std::runtime_error(oss.str());
        }
        }
        return {n0, n1, n2, n3};
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

    FaceKey make_face_key(
        const FaceLocalID &f,
        const std::unordered_map<LocalNodeID, NodeEqID, LocalNodeID::Hash> &node2eq,
        int8_t &sign_to_canonical)
    {
        auto local_corners = corners(f);
        std::array<NodeEqID, 4> corner_eq{};

        auto lookup = [&](const LocalNodeID &nid, NodeEqID &out) -> bool
        {
            auto it = node2eq.find(nid);
            if (it == node2eq.end())
                return false;
            out = it->second;
            return true;
        };

        bool all_found = true;
        for (int n = 0; n < 4; ++n)
        {
            all_found = lookup(local_corners[n], corner_eq[n]) && all_found;
        }

        if (!all_found && (f.dir == 1 || f.dir == 2))
        {
            // 2D grids have a single k plane. Represent xi/eta line faces as
            // degenerate four-corner faces with repeated endpoints.
            LocalNodeID a{f.rank, f.gblock, f.i, f.j, f.k};
            LocalNodeID b = a;
            if (f.dir == 1)
                ++b.j;
            else
                ++b.i;

            NodeEqID ga{};
            NodeEqID gb{};
            if (lookup(a, ga) && lookup(b, gb))
            {
                corner_eq = {ga, gb, ga, gb};
                all_found = true;
            }
        }

        if (!all_found)
        {
            throw std::runtime_error("TOPO::make_face_key: corner not found in node2eq.");
        }

        auto sorted = corner_eq;
        std::sort(sorted.begin(), sorted.end());

        const bool degenerate =
            (sorted[0] == sorted[1]) ||
            (sorted[1] == sorted[2]) ||
            (sorted[2] == sorted[3]);

        if (degenerate)
        {
            sign_to_canonical = +1;
        }
        else
        {
            std::array<int, 4> perm{};
            for (int n = 0; n < 4; ++n)
            {
                auto it = std::find(sorted.begin(), sorted.end(), corner_eq[n]);
                if (it == sorted.end())
                    throw std::runtime_error("TOPO::make_face_key: internal corner permutation error.");
                perm[n] = static_cast<int>(it - sorted.begin());
            }

            int inversions = 0;
            for (int a = 0; a < 4; ++a)
                for (int b = a + 1; b < 4; ++b)
                    if (perm[a] > perm[b])
                        ++inversions;

            sign_to_canonical = (inversions % 2 == 0) ? +1 : -1;
        }

        return FaceKey{sorted[0], sorted[1], sorted[2], sorted[3]};
    }

    namespace
    {
        StaggerLocation edge_location_from_dir(int dir)
        {
            switch (dir)
            {
            case 1:
                return StaggerLocation::EdgeXi;
            case 2:
                return StaggerLocation::EdgeEt;
            case 3:
                return StaggerLocation::EdgeZe;
            default:
                ERROR::Abort("TopologyEquiv: invalid edge dir");
                return StaggerLocation::Cell;
            }
        }

        EquivMember make_edge_member(const EdgeLocalID &e, int orient_sign, bool is_owner)
        {
            EquivMember m;
            m.rank = e.rank;
            m.block = e.gblock;
            m.location = edge_location_from_dir(e.dir);
            m.i = e.i;
            m.j = e.j;
            m.k = e.k;
            m.orient_sign = orient_sign;
            m.is_owner = is_owner;
            return m;
        }

        StaggerLocation face_location_from_dir(int dir)
        {
            switch (dir)
            {
            case 1:
                return StaggerLocation::FaceXi;
            case 2:
                return StaggerLocation::FaceEt;
            case 3:
                return StaggerLocation::FaceZe;
            default:
                ERROR::Abort("TopologyEquiv: invalid face dir");
                return StaggerLocation::Cell;
            }
        }

        EquivMember make_face_member(const FaceLocalID &f, int orient_sign, bool is_owner)
        {
            EquivMember m;
            m.rank = f.rank;
            m.block = f.gblock;
            m.location = face_location_from_dir(f.dir);
            m.i = f.i;
            m.j = f.j;
            m.k = f.k;
            m.orient_sign = orient_sign;
            m.is_owner = is_owner;
            return m;
        }
    }

    const std::vector<EquivClass> &TopologyEquiv::classes(EquivDofKind kind) const
    {
        switch (kind)
        {
        case EquivDofKind::Node:
            return node_classes;
        case EquivDofKind::Edge:
            return edge_classes_general;
        case EquivDofKind::Face:
            return face_classes;
        }

        ERROR::Abort("TopologyEquiv::classes: invalid EquivDofKind");
        return node_classes;
    }

    void TopologyEquiv::mirror_legacy_edge_equiv_to_general()
    {
        edge_classes_general.clear();
        edge_classes_general.reserve(edge_members.size());

        for (const auto &[key, members] : edge_members)
        {
            EquivClass cls;
            cls.kind = EquivDofKind::Edge;

            auto owner_it = edge_owner.find(key);
            const bool has_owner = (owner_it != edge_owner.end());
            EdgeLocalID owner{};

            if (has_owner)
            {
                owner = owner_it->second;
                auto gid_it = edge_owner_gid.find(owner);
                if (gid_it != edge_owner_gid.end())
                    cls.global_id = gid_it->second;

                int owner_sign = +1;
                auto sign_it = edge2sign.find(owner);
                if (sign_it != edge2sign.end())
                    owner_sign = static_cast<int>(sign_it->second);

                cls.owner = make_edge_member(owner, owner_sign, true);
            }

            cls.members.reserve(members.size());
            for (const auto &e : members)
            {
                int orient_sign = +1;
                auto sign_it = edge2sign.find(e);
                if (sign_it != edge2sign.end())
                    orient_sign = static_cast<int>(sign_it->second);

                bool is_owner = false;
                auto owner_flag_it = edge_is_owner.find(e);
                if (owner_flag_it != edge_is_owner.end())
                    is_owner = owner_flag_it->second;
                else if (has_owner)
                    is_owner = (e == owner);

                cls.members.push_back(make_edge_member(e, orient_sign, is_owner));
            }

            edge_classes_general.push_back(cls);
        }
    }

    void TopologyEquiv::mirror_legacy_face_equiv_to_general()
    {
        face_classes.clear();
        face_classes.reserve(face_members.size());

        for (const auto &[key, members] : face_members)
        {
            EquivClass cls;
            cls.kind = EquivDofKind::Face;

            auto owner_it = face_owner.find(key);
            const bool has_owner = (owner_it != face_owner.end());
            FaceLocalID owner{};

            if (has_owner)
            {
                owner = owner_it->second;
                auto gid_it = face_owner_gid.find(owner);
                if (gid_it != face_owner_gid.end())
                    cls.global_id = gid_it->second;

                int owner_sign = +1;
                auto sign_it = face2sign.find(owner);
                if (sign_it != face2sign.end())
                    owner_sign = static_cast<int>(sign_it->second);

                cls.owner = make_face_member(owner, owner_sign, true);
            }

            cls.members.reserve(members.size());
            for (const auto &f : members)
            {
                int orient_sign = +1;
                auto sign_it = face2sign.find(f);
                if (sign_it != face2sign.end())
                    orient_sign = static_cast<int>(sign_it->second);

                bool is_owner = false;
                auto owner_flag_it = face_is_owner.find(f);
                if (owner_flag_it != face_is_owner.end())
                    is_owner = owner_flag_it->second;
                else if (has_owner)
                    is_owner = (f == owner);

                cls.members.push_back(make_face_member(f, orient_sign, is_owner));
            }

            face_classes.push_back(cls);
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
        std::unordered_map<NodeEqID, int, NodeEqID::Hash> node_eq_counts;

        auto all_local_nodes = collect_all_local_nodes_impl(grid, my_rank);
        reconcile_node_equivalence_parallel_impl(topo, my_rank, all_local_nodes, equiv, node_eq_counts);

        auto all_local_edges = collect_all_local_edges_impl(grid, my_rank, dimension);
        build_edge_equivalence_from_nodes_impl(all_local_edges, equiv);

        select_edge_owner_parallel_impl(equiv, node_eq_counts);
        build_edge_owner_gid_impl(my_rank, equiv);
        equiv.mirror_legacy_edge_equiv_to_general();

        auto all_local_faces = collect_all_local_faces_impl(grid, my_rank, dimension);
        build_face_equivalence_from_nodes_impl(all_local_faces, equiv);

        select_face_owner_parallel_impl(equiv, node_eq_counts);
        build_face_owner_gid_impl(my_rank, equiv);
        equiv.mirror_legacy_face_equiv_to_general();
    }

    void build_node_equivalence(
        const Topology &topo,
        Grid &grid,
        int my_rank,
        int dimension,
        TopologyEquiv &equiv)
    {
        (void)topo;
        (void)grid;
        (void)my_rank;
        (void)dimension;
        (void)equiv;
        // TODO:
        // Build node equivalence classes for node-based owner sync.
        // This is reserved for Halo OwnerSyncPolicy::NodeOwner.
        // Current behavior unchanged.
    }

    void build_face_equivalence(
        const Topology &topo,
        Grid &grid,
        int my_rank,
        int dimension,
        TopologyEquiv &equiv)
    {
        if (equiv.node2eq.empty())
        {
            std::unordered_map<NodeEqID, int, NodeEqID::Hash> node_eq_counts;
            auto all_local_nodes = collect_all_local_nodes_impl(grid, my_rank);
            reconcile_node_equivalence_parallel_impl(topo, my_rank, all_local_nodes, equiv, node_eq_counts);
        }
        std::unordered_map<NodeEqID, int, NodeEqID::Hash> node_eq_counts;
        for (const auto &[node, eq] : equiv.node2eq)
        {
            (void)node;
            ++node_eq_counts[eq];
        }

        auto all_local_faces = collect_all_local_faces_impl(grid, my_rank, dimension);
        build_face_equivalence_from_nodes_impl(all_local_faces, equiv);

        select_face_owner_parallel_impl(equiv, node_eq_counts);
        build_face_owner_gid_impl(my_rank, equiv);
        equiv.mirror_legacy_face_equiv_to_general();
    }

} // namespace TOPO
