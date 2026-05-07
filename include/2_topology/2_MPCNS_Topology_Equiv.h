#pragma once
#include <compare>
#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <unordered_map>
#include <utility>
#include <vector>

#include "0_basic/StaggerLocation.h"
#include "2_topology/TopologyTypes.h"

class Grid;

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

    // local unique face id on a rank/block
    // dir convention: 1->FaceXi, 2->FaceEt, 3->FaceZe
    struct FaceLocalID
    {
        int rank;
        int gblock;
        int i, j, k;
        int dir;

        auto operator<=>(const FaceLocalID &) const = default;

        struct Hash
        {
            std::size_t operator()(const FaceLocalID &x) const
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

    // canonical physical face key, with sorted canonical corner nodes
    struct FaceKey
    {
        NodeEqID a;
        NodeEqID b;
        NodeEqID c;
        NodeEqID d;

        auto operator<=>(const FaceKey &) const = default;

        struct Hash
        {
            std::size_t operator()(const FaceKey &x) const
            {
                std::size_t h = 0;
                hash_combine_inplace(h, NodeEqID::Hash{}(x.a));
                hash_combine_inplace(h, NodeEqID::Hash{}(x.b));
                hash_combine_inplace(h, NodeEqID::Hash{}(x.c));
                hash_combine_inplace(h, NodeEqID::Hash{}(x.d));
                return h;
            }
        };
    };

    // ============================================================
    // equivalence-class containers
    // ============================================================

    enum class EquivDofKind
    {
        Node,
        Edge,
        Face
    };

    struct EquivMember
    {
        int rank = 0;
        int block = -1;

        StaggerLocation location = StaggerLocation::Cell;

        int i = 0;
        int j = 0;
        int k = 0;

        // member 相对 canonical orientation 的符号。
        int orient_sign = +1;

        bool is_owner = false;
    };

    struct EquivClass
    {
        EquivDofKind kind = EquivDofKind::Node;

        int global_id = -1;

        EquivMember owner;
        std::vector<EquivMember> members;
    };

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

        std::unordered_map<FaceLocalID, FaceKey, FaceLocalID::Hash> face2key;
        std::unordered_map<FaceLocalID, int8_t, FaceLocalID::Hash> face2sign;
        std::unordered_map<FaceKey, std::vector<FaceLocalID>, FaceKey::Hash> face_members;
        std::unordered_map<FaceKey, FaceLocalID, FaceKey::Hash> face_owner;
        std::unordered_map<FaceLocalID, bool, FaceLocalID::Hash> face_is_owner;

        std::unordered_map<FaceLocalID, int, FaceLocalID::Hash> face_owner_gid;
        std::unordered_map<int, FaceLocalID> gid2face_owner;
        int n_local_face_owner = 0;
        int n_global_face_owner = 0;
        int face_owner_gid_begin = 0;
        int face_owner_gid_end = 0; // half-open: [begin, end)

        std::vector<EquivClass> node_classes;
        std::vector<EquivClass> edge_classes_general;
        std::vector<EquivClass> face_classes;

        bool has_node_equiv() const { return !node_classes.empty(); }
        bool has_edge_equiv() const
        {
            return !edge_classes_general.empty() || !edge_owner.empty();
        }
        bool has_face_equiv() const { return !face_classes.empty(); }

        const std::vector<EquivClass> &classes(EquivDofKind kind) const;

        void mirror_legacy_edge_equiv_to_general();
        void mirror_legacy_face_equiv_to_general();

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

            face2key.clear();
            face2sign.clear();
            face_members.clear();
            face_owner.clear();
            face_is_owner.clear();

            face_owner_gid.clear();
            gid2face_owner.clear();

            n_local_face_owner = 0;
            n_global_face_owner = 0;
            face_owner_gid_begin = 0;
            face_owner_gid_end = 0; // half-open: [begin, end)

            node_classes.clear();
            edge_classes_general.clear();
            face_classes.clear();
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

    // local face corners in oriented local order.
    // dir=1 (FaceXi) : (i,j,k), (i,j+1,k), (i,j,k+1), (i,j+1,k+1)
    // dir=2 (FaceEt) : (i,j,k), (i+1,j,k), (i,j,k+1), (i+1,j,k+1)
    // dir=3 (FaceZe) : (i,j,k), (i+1,j,k), (i,j+1,k), (i+1,j+1,k)
    std::array<LocalNodeID, 4> corners(const FaceLocalID &f);

    // build canonical EdgeKey from a local edge using node2eq
    // returns sign_to_canonical:
    //   +1 if local edge direction matches key.a -> key.b
    //   -1 otherwise
    EdgeKey make_edge_key(
        const EdgeLocalID &e,
        const std::unordered_map<LocalNodeID, NodeEqID, LocalNodeID::Hash> &node2eq,
        int8_t &sign_to_canonical);

    // build canonical FaceKey from a local face using node2eq
    // returns sign_to_canonical based on local corner ordering parity against
    // the sorted canonical corner ordering. Degenerate 2D faces use +1.
    FaceKey make_face_key(
        const FaceLocalID &f,
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

    void build_node_equivalence(
        const Topology &topo,
        Grid &grid,
        int my_rank,
        int dimension,
        TopologyEquiv &equiv);

    void build_face_equivalence(
        const Topology &topo,
        Grid &grid,
        int my_rank,
        int dimension,
        TopologyEquiv &equiv);

} // namespace TOPO
