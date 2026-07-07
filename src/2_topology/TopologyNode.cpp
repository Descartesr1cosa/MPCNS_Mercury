#include "2_topology/TopologyEquivDetail.h"
#include "2_topology/TopologyOps.h"

#include "0_basic/BoxOps.h"
#include "0_basic/MPI_WRAPPER.h"
#include "1_grid/1_MPCNS_Grid.h"

#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace TOPO::detail
{
    namespace
    {
        void pack_pair(std::vector<int> &buf,
                       const EntityKey &a,
                       const EntityKey &b)
        {
            pack_node(buf, a);
            pack_node(buf, b);
        }
    }

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

    
    std::vector<EntityKey> collect_all_local_nodes(
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

    
    void reconcile_node_equivalence_parallel(
        const Topology &topo,
        int my_rank,
        const std::vector<EntityKey> &all_local_nodes,
        Topology &equiv,
        std::unordered_map<EntityKey, int, EntityKey::Hash> &rep_count)
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
                oss << "build_topology: mapped node out of nb_box_node. "
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
        throw std::runtime_error("build_topology: gathered node-pair int count is not multiple of 10.");
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

        equiv.nodes.local_to_rep.clear();
        rep_count.clear();
        for (const auto &[root, rep] : root_min_node)
        {
        rep_count[rep] = root_node_count[root];
        }

        for (const auto &nid : all_local_nodes)
        {
        auto it = node2idx.find(nid);
        if (it == node2idx.end())
        {
            std::ostringstream oss;
            oss << "build_topology: local node missing in node2idx. "
            << "(" << nid.rank << "," << nid.block << ","
            << nid.i << "," << nid.j << "," << nid.k << ")";
            throw std::runtime_error(oss.str());
        }

        int root = uf.find(it->second);
        const EntityKey &rep = root_min_node.at(root);
        equiv.nodes.local_to_rep[nid] = rep;
        }
    }

    
    void build_node_entity_ids(Topology &equiv)
    {
        std::vector<EntityKey> local_keys;
        local_keys.reserve(equiv.nodes.local_to_rep.size());
        for (const auto &[node, eq] : equiv.nodes.local_to_rep)
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
        allgather_packed_records(send_buf, 5, "build_topology node EntityId");
        std::vector<EntityKey> global_keys;
        global_keys.reserve(recv_buf.size() / 5);
        for (std::size_t n = 0; n < recv_buf.size(); n += 5)
        global_keys.push_back(unpack_node(recv_buf.data() + n));

        std::sort(global_keys.begin(), global_keys.end());
        global_keys.erase(std::unique(global_keys.begin(), global_keys.end()), global_keys.end());
        equiv.nodes.rep_to_qid.clear();
        for (std::size_t n = 0; n < global_keys.size(); ++n)
        equiv.nodes.rep_to_qid[global_keys[n]] = static_cast<int>(n);
    }

    } // namespace TOPO::detail
