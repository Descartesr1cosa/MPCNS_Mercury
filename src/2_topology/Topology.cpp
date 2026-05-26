#include "2_topology/Topology.h"
#include "2_topology/LocalIncidence.h"
#include "2_topology/TopologyOps.h"

#include "0_basic/MPI_WRAPPER.h"
#include "0_basic/BoxOps.h"
#include "0_basic/Error.h"
#include "1_grid/1_MPCNS_Grid.h"

#include <algorithm>
#include <array>
#include <map>
#include <ostream>
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

        inline void pack_node(std::vector<int> &buf, const EntityKey &x)
        {
            buf.push_back(x.rank);
            buf.push_back(x.block);
            buf.push_back(x.i);
            buf.push_back(x.j);
            buf.push_back(x.k);
        }

        inline EntityKey unpack_node(const int *p)
        {
            return make_node(p[0], p[1], p[2], p[3], p[4]);
        }

        inline void pack_pair(std::vector<int> &buf,
                              const EntityKey &a,
                              const EntityKey &b)
        {
            pack_node(buf, a);
            pack_node(buf, b);
        }

        std::vector<int> allgather_packed_records(std::vector<int> send_buf,
                                                  int record_size,
                                                  const char *context)
        {
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

            if (total_recv % record_size != 0)
            {
                std::ostringstream oss;
                oss << context << ": gathered int count is not multiple of " << record_size << ".";
                throw std::runtime_error(oss.str());
            }
            return recv_buf;
        }

        // ============================================================
        // local collection
        // ============================================================

        std::vector<EntityKey> collect_all_local_nodes_impl(
            Grid &grid,
            int my_rank)
        {
            std::vector<EntityKey> nodes;

            for (int ib = 0; ib < grid.nblock; ++ib)
            {
                const auto &blk = grid.grids(ib);

                for (int i = 0; i <= blk.mx; ++i)
                    for (int j = 0; j <= blk.my; ++j)
                        for (int k = 0; k <= blk.mz; ++k)
                        {
                            nodes.push_back(make_node(my_rank, ib, i, j, k));
                        }
            }
            return nodes;
        }

        std::vector<EntityKey> collect_all_local_edges_impl(
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

        std::vector<EntityKey> collect_all_local_faces_impl(
            Grid &grid,
            int my_rank,
            int dimension)
        {
            std::vector<EntityKey> faces;

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
                                faces.push_back(make_face(my_rank, ib, i, j, k, EntityAxis::Xi));
                            }

                    // FaceEt
                    for (int i = 0; i < blk.mx; ++i)
                        for (int j = 0; j <= blk.my; ++j)
                            for (int k = 0; k < blk.mz; ++k)
                            {
                                faces.push_back(make_face(my_rank, ib, i, j, k, EntityAxis::Eta));
                            }

                    // FaceZe
                    for (int i = 0; i < blk.mx; ++i)
                        for (int j = 0; j < blk.my; ++j)
                            for (int k = 0; k <= blk.mz; ++k)
                            {
                                faces.push_back(make_face(my_rank, ib, i, j, k, EntityAxis::Zeta));
                            }
                }
                else if (dimension >= 2)
                {
                    // 2D uses xi/eta face locations as line faces in the single k plane.
                    const int k = 0;

                    for (int i = 0; i <= blk.mx; ++i)
                        for (int j = 0; j < blk.my; ++j)
                        {
                            faces.push_back(make_face(my_rank, ib, i, j, k, EntityAxis::Xi));
                        }

                    for (int i = 0; i < blk.mx; ++i)
                        for (int j = 0; j <= blk.my; ++j)
                        {
                            faces.push_back(make_face(my_rank, ib, i, j, k, EntityAxis::Eta));
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
            const std::vector<EntityKey> &all_local_nodes,
            Topology &equiv,
            std::unordered_map<EntityKey, int, EntityKey::Hash> &node_eq_counts)
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
                                oss << "build_topology_equivalence: mapped node out of nb_box_node. "
                                    << "this_rank=" << patch.this_rank
                                    << " nb_rank=" << patch.nb_rank
                                    << " this_block=" << patch.this_block
                                    << " nb_block=" << patch.nb_block
                                    << " p_this=(" << p_this.i << "," << p_this.j << "," << p_this.k << ")"
                                    << " p_nb=(" << p_nb.i << "," << p_nb.j << "," << p_nb.k << ")";
                                throw std::runtime_error(oss.str());
                            }

                            EntityKey a = make_node(patch.this_rank, patch.this_block, i, j, k);
                            EntityKey b = make_node(patch.nb_rank, patch.nb_block, p_nb.i, p_nb.j, p_nb.k);
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
                throw std::runtime_error("build_topology_equivalence: gathered node-pair int count is not multiple of 10.");
            }

            std::unordered_map<EntityKey, int, EntityKey::Hash> node2idx;
            std::vector<EntityKey> idx2node;
            UnionFind uf;

            auto get_or_add_node = [&](const EntityKey &nid) -> int
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
                EntityKey a = unpack_node(base + 0);
                EntityKey b = unpack_node(base + 5);

                int ia = get_or_add_node(a);
                int ib = get_or_add_node(b);
                uf.unite(ia, ib);
            }

            std::unordered_map<int, EntityKey> root_min_node;
            std::unordered_map<int, int> root_node_count;
            for (int idx = 0; idx < static_cast<int>(idx2node.size()); ++idx)
            {
                int root = uf.find(idx);
                const EntityKey &nid = idx2node[idx];

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
                node_eq_counts[rep] = root_node_count[root];
            }

            for (const auto &nid : all_local_nodes)
            {
                auto it = node2idx.find(nid);
                if (it == node2idx.end())
                {
                    std::ostringstream oss;
                    oss << "build_topology_equivalence: local node missing in node2idx. "
                        << "(" << nid.rank << "," << nid.block << ","
                        << nid.i << "," << nid.j << "," << nid.k << ")";
                    throw std::runtime_error(oss.str());
                }

                int root = uf.find(it->second);
                const EntityKey &rep = root_min_node.at(root);
                equiv.node2eq[nid] = rep;
            }
        }

        void build_node_entity_ids_impl(Topology &equiv)
        {
            std::vector<EntityKey> local_keys;
            local_keys.reserve(equiv.node2eq.size());
            for (const auto &[node, eq] : equiv.node2eq)
            {
                (void)node;
                local_keys.push_back(eq);
            }
            std::sort(local_keys.begin(), local_keys.end());
            local_keys.erase(std::unique(local_keys.begin(), local_keys.end()), local_keys.end());

            std::vector<int> send_buf;
            send_buf.reserve(local_keys.size() * 5);
            for (const EntityKey &eq : local_keys)
                pack_node(send_buf, eq);

            const std::vector<int> recv_buf =
                allgather_packed_records(send_buf, 5, "build_topology_equivalence node EntityId");
            std::vector<EntityKey> global_keys;
            global_keys.reserve(recv_buf.size() / 5);
            for (std::size_t n = 0; n < recv_buf.size(); n += 5)
                global_keys.push_back(unpack_node(recv_buf.data() + n));

            std::sort(global_keys.begin(), global_keys.end());
            global_keys.erase(std::unique(global_keys.begin(), global_keys.end()), global_keys.end());
            equiv.node_eq_to_id.clear();
            for (std::size_t n = 0; n < global_keys.size(); ++n)
                equiv.node_eq_to_id[global_keys[n]] = static_cast<int>(n);
        }

        inline int node_equiv_count(
            const std::unordered_map<EntityKey, int, EntityKey::Hash> &node_eq_counts,
            const EntityKey &node)
        {
            auto it = node_eq_counts.find(node);
            return (it != node_eq_counts.end()) ? it->second : 1;
        }

        inline bool edge_key_is_shared(
            const EdgeKey &key,
            const std::unordered_map<EntityKey, int, EntityKey::Hash> &node_eq_counts)
        {
            return node_equiv_count(node_eq_counts, key.a) > 1 ||
                   node_equiv_count(node_eq_counts, key.b) > 1;
        }

        inline bool face_key_is_shared(
            const FaceKey &key,
            const std::unordered_map<EntityKey, int, EntityKey::Hash> &node_eq_counts)
        {
            return node_equiv_count(node_eq_counts, key.a) > 1 ||
                   node_equiv_count(node_eq_counts, key.b) > 1 ||
                   node_equiv_count(node_eq_counts, key.c) > 1 ||
                   node_equiv_count(node_eq_counts, key.d) > 1;
        }

        void build_edge_equivalence_from_nodes_impl(
            const std::vector<EntityKey> &all_local_edges,
            Topology &equiv)
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

        inline void pack_edge_local(std::vector<int> &buf, const EntityKey &e)
        {
            buf.push_back(e.rank);
            buf.push_back(e.block);
            buf.push_back(e.i);
            buf.push_back(e.j);
            buf.push_back(e.k);
            buf.push_back(axis_number(e.axis));
        }

        inline EntityKey unpack_edge_local(const int *p)
        {
            return make_edge(p[0], p[1], p[2], p[3], p[4], entity_axis(p[5]));
        }

        inline void pack_edge_key(std::vector<int> &buf, const EdgeKey &key)
        {
            pack_node(buf, key.a);
            pack_node(buf, key.b);
        }

        inline EdgeKey unpack_edge_key(const int *p)
        {
            EntityKey a = unpack_node(p + 0);
            EntityKey b = unpack_node(p + 5);
            return EdgeKey{a, b};
        }

        void build_edge_entity_ids_impl(Topology &equiv)
        {
            std::vector<EdgeKey> local_keys;
            local_keys.reserve(equiv.edge2key.size());
            for (const auto &[edge, key] : equiv.edge2key)
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
                allgather_packed_records(send_buf, 10, "build_topology_equivalence edge EntityId");
            std::vector<EdgeKey> global_keys;
            global_keys.reserve(recv_buf.size() / 10);
            for (std::size_t n = 0; n < recv_buf.size(); n += 10)
                global_keys.push_back(unpack_edge_key(recv_buf.data() + n));

            std::sort(global_keys.begin(), global_keys.end());
            global_keys.erase(std::unique(global_keys.begin(), global_keys.end()), global_keys.end());
            equiv.edge_key_to_id.clear();
            for (std::size_t n = 0; n < global_keys.size(); ++n)
                equiv.edge_key_to_id[global_keys[n]] = static_cast<int>(n);
        }

        inline void pack_edge_owner_candidate(
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

        inline void select_edge_owner_parallel_impl(
            Topology &equiv,
            const std::unordered_map<EntityKey, int, EntityKey::Hash> &node_eq_counts)
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
                    const auto sign_it = equiv.edge2sign.find(e);
                    const int8_t sign = (sign_it != equiv.edge2sign.end()) ? sign_it->second : int8_t{+1};
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
                    "build_topology_equivalence: gathered edge-candidate int count is not multiple of 17.");
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

                equiv.edge2key[e] = key;
                equiv.edge2sign[e] = sign;
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
                        "build_topology_equivalence: shared edge key missing in global_owner.");
                }

                const EntityKey &owner = it->second;
                equiv.edge_owner[key] = owner;
                equiv.edge_members[key] = members;

                for (const auto &e : members)
                {
                    equiv.edge_is_owner[e] = (e == owner);
                }
            }
        }

        inline void build_edge_owner_gid_impl(int my_rank, Topology &equiv)
        {
            std::vector<EntityKey> local_owner_edges;
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
                const EntityKey &e = local_owner_edges[n];
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
                    "build_topology_equivalence: gathered edge-owner-gid int count is not multiple of 7.");
            }

            const int nowners = total_recv / 7;
            for (int n = 0; n < nowners; ++n)
            {
                const int *base = recv_buf.data() + 7 * n;
                EntityKey e = unpack_edge_local(base);
                int gid = base[6];

                equiv.edge_owner_gid[e] = gid;
                equiv.gid2edge_owner[gid] = e;
            }
        }

        void build_face_equivalence_from_nodes_impl(
            const std::vector<EntityKey> &all_local_faces,
            Topology &equiv)
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

        inline void pack_face_local(std::vector<int> &buf, const EntityKey &f)
        {
            buf.push_back(f.rank);
            buf.push_back(f.block);
            buf.push_back(f.i);
            buf.push_back(f.j);
            buf.push_back(f.k);
            buf.push_back(axis_number(f.axis));
        }

        inline EntityKey unpack_face_local(const int *p)
        {
            return make_face(p[0], p[1], p[2], p[3], p[4], entity_axis(p[5]));
        }

        inline void pack_face_key(std::vector<int> &buf, const FaceKey &key)
        {
            pack_node(buf, key.a);
            pack_node(buf, key.b);
            pack_node(buf, key.c);
            pack_node(buf, key.d);
        }

        inline FaceKey unpack_face_key(const int *p)
        {
            EntityKey a = unpack_node(p + 0);
            EntityKey b = unpack_node(p + 5);
            EntityKey c = unpack_node(p + 10);
            EntityKey d = unpack_node(p + 15);
            return FaceKey{a, b, c, d};
        }

        void build_face_entity_ids_impl(Topology &equiv)
        {
            std::vector<FaceKey> local_keys;
            local_keys.reserve(equiv.face2key.size());
            for (const auto &[face, key] : equiv.face2key)
            {
                (void)face;
                local_keys.push_back(key);
            }
            std::sort(local_keys.begin(), local_keys.end());
            local_keys.erase(std::unique(local_keys.begin(), local_keys.end()), local_keys.end());

            std::vector<int> send_buf;
            send_buf.reserve(local_keys.size() * 20);
            for (const FaceKey &key : local_keys)
                pack_face_key(send_buf, key);

            const std::vector<int> recv_buf =
                allgather_packed_records(send_buf, 20, "build_topology_equivalence face EntityId");
            std::vector<FaceKey> global_keys;
            global_keys.reserve(recv_buf.size() / 20);
            for (std::size_t n = 0; n < recv_buf.size(); n += 20)
                global_keys.push_back(unpack_face_key(recv_buf.data() + n));

            std::sort(global_keys.begin(), global_keys.end());
            global_keys.erase(std::unique(global_keys.begin(), global_keys.end()), global_keys.end());
            equiv.face_key_to_id.clear();
            for (std::size_t n = 0; n < global_keys.size(); ++n)
                equiv.face_key_to_id[global_keys[n]] = static_cast<int>(n);
        }

        inline void pack_face_owner_candidate(
            std::vector<int> &buf,
            const FaceKey &key,
            const EntityKey &f,
            int8_t sign)
        {
            // 20 ints for FaceKey + 6 ints for EntityKey + 1 sign = 27 ints
            pack_face_key(buf, key);
            pack_face_local(buf, f);
            buf.push_back(static_cast<int>(sign));
        }

        inline void select_face_owner_parallel_impl(
            Topology &equiv,
            const std::unordered_map<EntityKey, int, EntityKey::Hash> &node_eq_counts)
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
                    "build_topology_equivalence: gathered face-candidate int count is not multiple of 27.");
            }

            std::unordered_map<FaceKey, EntityKey, FaceKey::Hash> global_owner;
            std::unordered_map<FaceKey, std::vector<EntityKey>, FaceKey::Hash> global_members;

            const int ncand = total_recv / 27;
            for (int c = 0; c < ncand; ++c)
            {
                const int *base = recv_buf.data() + 27 * c;

                FaceKey key = unpack_face_key(base + 0);
                EntityKey f = unpack_face_local(base + 20);
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
                        "build_topology_equivalence: local face key missing in global_owner.");
                }

                const EntityKey &owner = it->second;
                equiv.face_owner[key] = owner;
                equiv.face_members[key] = members;

                for (const auto &f : members)
                {
                    equiv.face_is_owner[f] = (f == owner);
                }
            }

        }

        inline void build_face_owner_gid_impl(int my_rank, Topology &equiv)
        {
            std::vector<EntityKey> local_owner_faces;
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
                const EntityKey &f = local_owner_faces[n];
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
                    "build_topology_equivalence: gathered face-owner-gid int count is not multiple of 7.");
            }

            const int nowners = total_recv / 7;
            for (int n = 0; n < nowners; ++n)
            {
                const int *base = recv_buf.data() + 7 * n;
                EntityKey f = unpack_face_local(base);
                int gid = base[6];

                equiv.face_owner_gid[f] = gid;
                equiv.gid2face_owner[gid] = f;
            }
        }
    } // anonymous namespace

    std::pair<EntityKey, EntityKey> endpoints(const EntityKey &e)
    {
        EntityKey n0 = make_node(e.rank, e.block, e.i, e.j, e.k);
        EntityKey n1 = n0;

        switch (e.axis)
        {
        case EntityAxis::Xi:
            ++n1.i;
            break; // Xi
        case EntityAxis::Eta:
            ++n1.j;
            break; // Eta
        case EntityAxis::Zeta:
            ++n1.k;
            break; // Zeta
        default:
        {
            std::ostringstream oss;
            oss << "TOPO::endpoints: invalid edge axis";
            throw std::runtime_error(oss.str());
        }
        }
        return {n0, n1};
    }

    std::array<EntityKey, 4> corners(const EntityKey &f)
    {
        EntityKey n0 = make_node(f.rank, f.block, f.i, f.j, f.k);
        EntityKey n1 = n0;
        EntityKey n2 = n0;
        EntityKey n3 = n0;

        switch (f.axis)
        {
        case EntityAxis::Xi:
            ++n1.j;
            ++n2.k;
            ++n3.j;
            ++n3.k;
            break; // FaceXi
        case EntityAxis::Eta:
            ++n1.i;
            ++n2.k;
            ++n3.i;
            ++n3.k;
            break; // FaceEt
        case EntityAxis::Zeta:
            ++n1.i;
            ++n2.j;
            ++n3.i;
            ++n3.j;
            break; // FaceZe
        default:
        {
            std::ostringstream oss;
            oss << "TOPO::corners: invalid face axis";
            throw std::runtime_error(oss.str());
        }
        }
        return {n0, n1, n2, n3};
    }

    EdgeKey make_edge_key(
        const EntityKey &e,
        const std::unordered_map<EntityKey, EntityKey, EntityKey::Hash> &node2eq,
        int8_t &sign_to_canonical)
    {
        const auto [ln0, ln1] = endpoints(e);

        auto it0 = node2eq.find(ln0);
        auto it1 = node2eq.find(ln1);

        if (it0 == node2eq.end() || it1 == node2eq.end())
        {
            throw std::runtime_error("TOPO::make_edge_key: endpoint not found in node2eq.");
        }

        const EntityKey &g0 = it0->second;
        const EntityKey &g1 = it1->second;

        if (g0 == g1)
        {
            throw std::runtime_error("TOPO::make_edge_key: two endpoints collapse to the same EntityKey.");
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
        const EntityKey &f,
        const std::unordered_map<EntityKey, EntityKey, EntityKey::Hash> &node2eq,
        int8_t &sign_to_canonical)
    {
        // Retain the established convention: the sorted corner EntityKeys are
        // the physical FaceKey, and permutation parity supplies local sign
        // relative to that canonical ordering.  Strict DEC validation is
        // supplied below rather than changing this rule speculatively.
        auto local_corners = corners(f);
        std::array<EntityKey, 4> corner_eq{};

        auto lookup = [&](const EntityKey &nid, EntityKey &out) -> bool
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

        if (!all_found && (f.axis == EntityAxis::Xi || f.axis == EntityAxis::Eta))
        {
            // 2D grids have a single k plane. Represent xi/eta line faces as
            // degenerate four-corner faces with repeated endpoints.
            EntityKey a = make_node(f.rank, f.block, f.i, f.j, f.k);
            EntityKey b = a;
            if (f.axis == EntityAxis::Xi)
                ++b.j;
            else
                ++b.i;

            EntityKey ga{};
            EntityKey gb{};
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
        EquivMember make_edge_member(const EntityKey &e, int orient_sign, bool is_owner)
        {
            EquivMember m;
            m.entity = e;
            m.orient_sign = orient_sign;
            m.is_owner = is_owner;
            return m;
        }

        EquivMember make_face_member(const EntityKey &f, int orient_sign, bool is_owner)
        {
            EquivMember m;
            m.entity = f;
            m.orient_sign = orient_sign;
            m.is_owner = is_owner;
            return m;
        }
    }

    EntityId Topology::id_of(const EntityKey &key) const
    {
        switch (key.dim)
        {
        case EntityDim::Node:
        {
            const EntityKey local = key;
            const auto eq_it = node2eq.find(local);
            if (eq_it == node2eq.end())
                throw std::runtime_error("Topology::id_of: node is not present in node2eq.");
            const auto id_it = node_eq_to_id.find(eq_it->second);
            if (id_it == node_eq_to_id.end())
                throw std::runtime_error("Topology::id_of: node quotient id is not built.");
            return EntityId{EntityDim::Node, id_it->second};
        }
        case EntityDim::Edge:
        {
            const EntityKey local = key;
            const auto key_it = edge2key.find(local);
            if (key_it == edge2key.end())
                throw std::runtime_error("Topology::id_of: edge is not present in edge2key.");
            const auto id_it = edge_key_to_id.find(key_it->second);
            if (id_it == edge_key_to_id.end())
                throw std::runtime_error("Topology::id_of: edge quotient id is not built.");
            return EntityId{EntityDim::Edge, id_it->second};
        }
        case EntityDim::Face:
        {
            const EntityKey local = key;
            const auto key_it = face2key.find(local);
            if (key_it == face2key.end())
                throw std::runtime_error("Topology::id_of: face is not present in face2key.");
            const auto id_it = face_key_to_id.find(key_it->second);
            if (id_it == face_key_to_id.end())
                throw std::runtime_error("Topology::id_of: face quotient id is not built.");
            return EntityId{EntityDim::Face, id_it->second};
        }
        case EntityDim::Cell:
            throw std::runtime_error(
                "Topology::id_of: cell quotient ids are not implemented in this phase.");
        }
        throw std::runtime_error("Topology::id_of: invalid entity dimension.");
    }

    EntityKey Topology::owner_of(const EntityKey &key) const
    {
        switch (key.dim)
        {
        case EntityDim::Node:
        {
            const auto it = node2eq.find(key);
            if (it == node2eq.end())
                throw std::runtime_error("Topology::owner_of: node is not present in node2eq.");
            return make_node(it->second.rank, it->second.block,
                             it->second.i, it->second.j, it->second.k);
        }
        case EntityDim::Edge:
        {
            const EntityKey local = key;
            const auto key_it = edge2key.find(local);
            if (key_it == edge2key.end())
                throw std::runtime_error("Topology::owner_of: edge is not present in edge2key.");
            const auto owner_it = edge_owner.find(key_it->second);
            return owner_it == edge_owner.end() ? key : owner_it->second;
        }
        case EntityDim::Face:
        {
            const EntityKey local = key;
            const auto key_it = face2key.find(local);
            if (key_it == face2key.end())
                throw std::runtime_error("Topology::owner_of: face is not present in face2key.");
            const auto owner_it = face_owner.find(key_it->second);
            return owner_it == face_owner.end() ? key : owner_it->second;
        }
        case EntityDim::Cell:
            throw std::runtime_error(
                "Topology::owner_of: cell owner-alias is not implemented in this phase.");
        }
        throw std::runtime_error("Topology::owner_of: invalid entity dimension.");
    }

    int Topology::sign_to_owner(const EntityKey &key) const
    {
        // Existing maps store member signs relative to canonical key
        // orientation.  Convert that convention to member relative to the
        // selected owner: member_sign / owner_sign == member_sign*owner_sign.
        switch (key.dim)
        {
        case EntityDim::Node:
            (void)owner_of(key);
            return +1;
        case EntityDim::Edge:
        {
            const EntityKey local = key;
            const auto sign_it = edge2sign.find(local);
            if (sign_it == edge2sign.end())
                throw std::runtime_error("Topology::sign_to_owner: edge sign is not present.");
            const auto key_it = edge2key.find(local);
            if (key_it == edge2key.end())
                throw std::runtime_error("Topology::sign_to_owner: edge key is not present.");
            const auto owner_it = edge_owner.find(key_it->second);
            if (owner_it == edge_owner.end())
                return +1;
            const auto owner_sign_it = edge2sign.find(owner_it->second);
            if (owner_sign_it == edge2sign.end())
                throw std::runtime_error("Topology::sign_to_owner: owner edge sign is not present.");
            return static_cast<int>(sign_it->second) *
                   static_cast<int>(owner_sign_it->second);
        }
        case EntityDim::Face:
        {
            const EntityKey local = key;
            const auto sign_it = face2sign.find(local);
            if (sign_it == face2sign.end())
                throw std::runtime_error("Topology::sign_to_owner: face sign is not present.");
            const auto key_it = face2key.find(local);
            if (key_it == face2key.end())
                throw std::runtime_error("Topology::sign_to_owner: face key is not present.");
            const auto owner_it = face_owner.find(key_it->second);
            if (owner_it == face_owner.end())
                return +1;
            const auto owner_sign_it = face2sign.find(owner_it->second);
            if (owner_sign_it == face2sign.end())
                throw std::runtime_error("Topology::sign_to_owner: owner face sign is not present.");
            return static_cast<int>(sign_it->second) *
                   static_cast<int>(owner_sign_it->second);
        }
        case EntityDim::Cell:
            throw std::runtime_error(
                "Topology::sign_to_owner: cell owner-alias is not implemented in this phase.");
        }
        throw std::runtime_error("Topology::sign_to_owner: invalid entity dimension.");
    }

    bool Topology::is_owner(const EntityKey &key) const
    {
        return owner_of(key) == key;
    }

    namespace
    {
        void print_face_key(std::ostream &os, const FaceKey &key)
        {
            const EntityKey corners[4] = {key.a, key.b, key.c, key.d};
            os << "[";
            for (int n = 0; n < 4; ++n)
            {
                if (n != 0)
                    os << ", ";
                os << "(" << corners[n].rank << "," << corners[n].block << ","
                   << corners[n].i << "," << corners[n].j << "," << corners[n].k << ")";
            }
            os << "]";
        }

        void print_edge_stencil(std::ostream &os, const std::map<EntityId, int> &stencil)
        {
            os << "{";
            bool first = true;
            for (const auto &[edge, coefficient] : stencil)
            {
                if (!first)
                    os << ", ";
                first = false;
                os << "edge#" << edge.id << ":" << coefficient;
            }
            os << "}";
        }
    }

    bool validate_face_orientation_stencils(const Topology &equiv,
                                            std::ostream &diagnostics)
    {
        bool valid = true;
        for (const auto &[key, members] : equiv.face_members)
        {
            bool have_reference = false;
            std::map<EntityId, int> reference;
            std::map<EntityId, int> reference_raw;
            EntityKey reference_member{};
            int reference_sign = +1;

            for (const EntityKey &member : members)
            {
                std::map<EntityId, int> raw_stencil;
                std::map<EntityId, int> normalized_stencil;
                int face_sign = +1;
                try
                {
                    const EntityKey local_face = member;
                    face_sign = equiv.sign_to_owner(local_face);
                    for (const IncidenceEntry &local_edge : boundary_of_face(local_face))
                    {
                        const EntityId edge_id = equiv.id_of(local_edge.entity);
                        raw_stencil[edge_id] +=
                            local_edge.sign * equiv.sign_to_owner(local_edge.entity);
                    }

                    normalized_stencil = raw_stencil;
                    for (auto &entry : normalized_stencil)
                        entry.second *= face_sign;
                }
                catch (const std::exception &error)
                {
                    valid = false;
                    diagnostics << "Topology face orientation stencil validation unavailable: FaceKey=";
                    print_face_key(diagnostics, key);
                    diagnostics << "\n  member=(" << member.rank << "," << member.block << ","
                                << member.i << "," << member.j << "," << member.k << ",axis="
                                << axis_number(member.axis) << ") error=" << error.what() << "\n";
                    continue;
                }

                if (!have_reference)
                {
                    have_reference = true;
                    reference = normalized_stencil;
                    reference_raw = raw_stencil;
                    reference_member = member;
                    reference_sign = face_sign;
                    continue;
                }
                if (normalized_stencil == reference)
                    continue;

                valid = false;
                diagnostics << "Topology face orientation stencil mismatch: FaceKey=";
                print_face_key(diagnostics, key);
                diagnostics << "\n  reference member=(" << reference_member.rank << ","
                            << reference_member.block << "," << reference_member.i << ","
                            << reference_member.j << "," << reference_member.k << ",axis="
                            << axis_number(reference_member.axis) << ") local_sign=" << reference_sign << " raw=";
                print_edge_stencil(diagnostics, reference_raw);
                diagnostics << " normalized=";
                print_edge_stencil(diagnostics, reference);
                diagnostics << "\n  member=(" << member.rank << "," << member.block << ","
                            << member.i << "," << member.j << "," << member.k << ",axis="
                            << axis_number(member.axis) << ") local_sign=" << face_sign << " raw=";
                print_edge_stencil(diagnostics, raw_stencil);
                diagnostics << " normalized=";
                print_edge_stencil(diagnostics, normalized_stencil);
                diagnostics << "\n";
            }
        }
        return valid;
    }

    const std::vector<EquivClass> &Topology::classes(EquivDofKind kind) const
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

        ERROR::Abort("Topology::classes: invalid EquivDofKind");
        return node_classes;
    }

    void Topology::mirror_legacy_edge_equiv_to_general()
    {
        edge_classes_general.clear();
        edge_classes_general.reserve(edge_members.size());

        for (const auto &[key, members] : edge_members)
        {
            EquivClass cls;
            cls.kind = EquivDofKind::Edge;

            auto owner_it = edge_owner.find(key);
            const bool has_owner = (owner_it != edge_owner.end());
            EntityKey owner{};

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

    void Topology::mirror_legacy_face_equiv_to_general()
    {
        face_classes.clear();
        face_classes.reserve(face_members.size());

        for (const auto &[key, members] : face_members)
        {
            EquivClass cls;
            cls.kind = EquivDofKind::Face;

            auto owner_it = face_owner.find(key);
            const bool has_owner = (owner_it != face_owner.end());
            EntityKey owner{};

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

    void build_topology_equivalence(
        Topology &equiv,
        Grid &grid,
        int my_rank,
        int dimension)
    {
        equiv.clear_equivalence();
        std::unordered_map<EntityKey, int, EntityKey::Hash> node_eq_counts;

        auto all_local_nodes = collect_all_local_nodes_impl(grid, my_rank);
        reconcile_node_equivalence_parallel_impl(equiv, my_rank, all_local_nodes, equiv, node_eq_counts);
        build_node_entity_ids_impl(equiv);

        auto all_local_edges = collect_all_local_edges_impl(grid, my_rank, dimension);
        build_edge_equivalence_from_nodes_impl(all_local_edges, equiv);
        build_edge_entity_ids_impl(equiv);

        select_edge_owner_parallel_impl(equiv, node_eq_counts);
        build_edge_owner_gid_impl(my_rank, equiv);
        equiv.mirror_legacy_edge_equiv_to_general();

        auto all_local_faces = collect_all_local_faces_impl(grid, my_rank, dimension);
        build_face_equivalence_from_nodes_impl(all_local_faces, equiv);
        build_face_entity_ids_impl(equiv);

        select_face_owner_parallel_impl(equiv, node_eq_counts);
        build_face_owner_gid_impl(my_rank, equiv);
        equiv.mirror_legacy_face_equiv_to_general();
    }

    void build_node_equivalence(
        Topology &equiv,
        Grid &grid,
        int my_rank,
        int dimension)
    {
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
        Topology &equiv,
        Grid &grid,
        int my_rank,
        int dimension)
    {
        if (equiv.node2eq.empty())
        {
            std::unordered_map<EntityKey, int, EntityKey::Hash> node_eq_counts;
            auto all_local_nodes = collect_all_local_nodes_impl(grid, my_rank);
            reconcile_node_equivalence_parallel_impl(equiv, my_rank, all_local_nodes, equiv, node_eq_counts);
        }
        build_node_entity_ids_impl(equiv);
        std::unordered_map<EntityKey, int, EntityKey::Hash> node_eq_counts;
        for (const auto &[node, eq] : equiv.node2eq)
        {
            (void)node;
            ++node_eq_counts[eq];
        }

        auto all_local_faces = collect_all_local_faces_impl(grid, my_rank, dimension);
        build_face_equivalence_from_nodes_impl(all_local_faces, equiv);
        build_face_entity_ids_impl(equiv);

        select_face_owner_parallel_impl(equiv, node_eq_counts);
        build_face_owner_gid_impl(my_rank, equiv);
        equiv.mirror_legacy_face_equiv_to_general();
    }

} // namespace TOPO
