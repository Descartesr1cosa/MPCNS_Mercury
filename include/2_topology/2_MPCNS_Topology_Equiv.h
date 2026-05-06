#pragma once
#include <compare>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <unordered_map>
#include <utility>
#include <vector>

#include "2_topology/2_MPCNS_Topology.h"

namespace TOPO
{
    // ============================================================
    // hash helpers
    // ============================================================

    inline void hash_combine_inplace(std::size_t &seed, std::size_t v)
    {
        seed ^= v + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
    }

    template <class T>
    inline void hash_combine_inplace(std::size_t &seed, const T &v)
    {
        hash_combine_inplace(seed, std::hash<T>{}(v));
    }

    // ============================================================
    // basic ids / keys
    // ============================================================

    // local unique node id on a rank/block
    // note: gblock here can be the block id on this rank;
    // (rank, gblock) together are globally unique
    struct LocalNodeID
    {
        int rank;
        int gblock;
        int i, j, k;

        auto operator<=>(const LocalNodeID &) const = default;

        struct Hash
        {
            std::size_t operator()(const LocalNodeID &x) const
            {
                std::size_t h = 0;
                hash_combine_inplace(h, x.rank);
                hash_combine_inplace(h, x.gblock);
                hash_combine_inplace(h, x.i);
                hash_combine_inplace(h, x.j);
                hash_combine_inplace(h, x.k);
                return h;
            }
        };
    };

    // canonical representative of a node equivalence class
    struct NodeEqID
    {
        int rank;
        int gblock;
        int i, j, k;

        auto operator<=>(const NodeEqID &) const = default;

        struct Hash
        {
            std::size_t operator()(const NodeEqID &x) const
            {
                std::size_t h = 0;
                hash_combine_inplace(h, x.rank);
                hash_combine_inplace(h, x.gblock);
                hash_combine_inplace(h, x.i);
                hash_combine_inplace(h, x.j);
                hash_combine_inplace(h, x.k);
                return h;
            }
        };
    };

    // local unique edge id on a rank/block
    // dir convention: 1->Xi, 2->Eta, 3->Zeta
    struct EdgeLocalID
    {
        int rank;
        int gblock;
        int i, j, k;
        int dir;

        auto operator<=>(const EdgeLocalID &) const = default;

        struct Hash
        {
            std::size_t operator()(const EdgeLocalID &x) const
            {
                std::size_t h = 0;
                hash_combine_inplace(h, x.rank);
                hash_combine_inplace(h, x.gblock);
                hash_combine_inplace(h, x.i);
                hash_combine_inplace(h, x.j);
                hash_combine_inplace(h, x.k);
                hash_combine_inplace(h, x.dir);
                return h;
            }
        };
    };

    // canonical physical edge key, with ordered endpoints a < b
    struct EdgeKey
    {
        NodeEqID a;
        NodeEqID b;

        auto operator<=>(const EdgeKey &) const = default;

        struct Hash
        {
            std::size_t operator()(const EdgeKey &x) const
            {
                std::size_t h = 0;
                hash_combine_inplace(h, NodeEqID::Hash{}(x.a));
                hash_combine_inplace(h, NodeEqID::Hash{}(x.b));
                return h;
            }
        };
    };

    // ============================================================
    // equivalence-class containers
    // ============================================================

    struct TopologyEquiv
    {
        // local node -> canonical node equivalence id
        std::unordered_map<LocalNodeID, NodeEqID, LocalNodeID::Hash> node2eq;

        // local edge -> canonical edge key
        std::unordered_map<EdgeLocalID, EdgeKey, EdgeLocalID::Hash> edge2key;

        // local edge -> sign to canonical edge direction
        // +1 : local direction == key.a -> key.b
        // -1 : local direction == key.b -> key.a
        std::unordered_map<EdgeLocalID, int8_t, EdgeLocalID::Hash> edge2sign;

        // canonical edge key -> all local members on this rank
        std::unordered_map<EdgeKey, std::vector<EdgeLocalID>, EdgeKey::Hash> edge_members;

        // canonical edge key -> chosen owner rep
        std::unordered_map<EdgeKey, EdgeLocalID, EdgeKey::Hash> edge_owner;

        // local edge -> whether this rep is owner
        std::unordered_map<EdgeLocalID, bool, EdgeLocalID::Hash> edge_is_owner;

        std::unordered_map<EdgeLocalID, int, EdgeLocalID::Hash> edge_owner_gid;
        std::unordered_map<int, EdgeLocalID> gid2edge_owner;
        int n_local_edge_owner = 0;
        int n_global_edge_owner = 0;
        int edge_owner_gid_begin = 0;
        int edge_owner_gid_end = 0; // half-open: [begin, end)

        void clear()
        {
            node2eq.clear();
            edge2key.clear();
            edge2sign.clear();
            edge_members.clear();
            edge_owner.clear();
            edge_is_owner.clear();

            edge_owner_gid.clear();
            gid2edge_owner.clear();

            n_local_edge_owner = 0;
            n_global_edge_owner = 0;
            edge_owner_gid_begin = 0;
            edge_owner_gid_end = 0; // half-open: [begin, end)
        }
    };

    // ============================================================
    // small helpers
    // ============================================================

    inline NodeEqID to_node_eq_id(const LocalNodeID &x)
    {
        return NodeEqID{x.rank, x.gblock, x.i, x.j, x.k};
    }

    inline LocalNodeID to_local_node_id(const NodeEqID &x)
    {
        return LocalNodeID{x.rank, x.gblock, x.i, x.j, x.k};
    }

    // local edge endpoints
    // dir=1 (Xi)   : (i,j,k) -> (i+1,j,k)
    // dir=2 (Eta)  : (i,j,k) -> (i,j+1,k)
    // dir=3 (Zeta) : (i,j,k) -> (i,j,k+1)
    std::pair<LocalNodeID, LocalNodeID> endpoints(const EdgeLocalID &e);

    // build canonical EdgeKey from a local edge using node2eq
    // returns sign_to_canonical:
    //   +1 if local edge direction matches key.a -> key.b
    //   -1 otherwise
    EdgeKey make_edge_key(
        const EdgeLocalID &e,
        const std::unordered_map<LocalNodeID, NodeEqID, LocalNodeID::Hash> &node2eq,
        int8_t &sign_to_canonical);

    // ============================================================
    // single public build entry
    // ============================================================

    void build_topology_equiv(
        const Topology &topo,
        Grid &grid,
        int my_rank,
        int dimension,
        TopologyEquiv &equiv);

} // namespace TOPO