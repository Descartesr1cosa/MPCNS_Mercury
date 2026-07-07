#pragma once
#include <compare>
#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iosfwd>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "0_basic/StaggerLocation.h"
#include "2_topology/Entity.h"
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

    // EntityKey is used for local entities and canonical node representatives.

    // canonical physical edge key, with ordered endpoints a < b
    struct EdgeKey
    {
        EntityKey a;
        EntityKey b;

        auto operator<=>(const EdgeKey &) const = default;

        struct Hash
        {
            std::size_t operator()(const EdgeKey &x) const
            {
                std::size_t h = 0;
                hash_combine_inplace(h, EntityKey::Hash{}(x.a));
                hash_combine_inplace(h, EntityKey::Hash{}(x.b));
                return h;
            }
        };
    };

    // Canonical physical face key built from sorted corner EntityKeys.  The key
    // identifies the geometric face; orientation is assigned later by comparing
    // each member's edge-boundary stencil with the selected owner face.
    struct FaceKey
    {
        EntityKey a;
        EntityKey b;
        EntityKey c;
        EntityKey d;

        auto operator<=>(const FaceKey &) const = default;

        struct Hash
        {
            std::size_t operator()(const FaceKey &x) const
            {
                std::size_t h = 0;
                hash_combine_inplace(h, EntityKey::Hash{}(x.a));
                hash_combine_inplace(h, EntityKey::Hash{}(x.b));
                hash_combine_inplace(h, EntityKey::Hash{}(x.c));
                hash_combine_inplace(h, EntityKey::Hash{}(x.d));
                return h;
            }
        };
    };

    // ============================================================
    // equivalence-class containers
    // ============================================================

    struct EquivMember
    {
        EntityKey entity{EntityDim::Cell, 0, 0, 0, 0, 0, EntityAxis::None};

        // member 相对 canonical orientation 的符号。
        int orient_sign = +1;

        bool is_owner = false;
    };

    inline StaggerLocation stagger_location(const EntityKey &entity)
    {
        switch (entity.dim)
        {
        case EntityDim::Node:
            return StaggerLocation::Node;
        case EntityDim::Cell:
            return StaggerLocation::Cell;
        case EntityDim::Edge:
            switch (entity.axis)
            {
            case EntityAxis::Xi:
                return StaggerLocation::EdgeXi;
            case EntityAxis::Eta:
                return StaggerLocation::EdgeEt;
            case EntityAxis::Zeta:
                return StaggerLocation::EdgeZe;
            default:
                break;
            }
            break;
        case EntityDim::Face:
            switch (entity.axis)
            {
            case EntityAxis::Xi:
                return StaggerLocation::FaceXi;
            case EntityAxis::Eta:
                return StaggerLocation::FaceEt;
            case EntityAxis::Zeta:
                return StaggerLocation::FaceZe;
            default:
                break;
            }
            break;
        }
        throw std::invalid_argument("TOPO::stagger_location: entity has an invalid dimension/axis combination.");
    }

    struct EquivClass
    {
        EntityDim dim = EntityDim::Node;

        int global_id = -1;

        EquivMember owner;
        std::vector<EquivMember> members;
    };

    struct NodeTopology
    {
        // local node representative -> quotient-node canonical representative
        std::unordered_map<EntityKey, EntityKey, EntityKey::Hash> local_to_rep;
        // canonical representative -> dimension-local quotient id
        std::unordered_map<EntityKey, int, EntityKey::Hash> rep_to_qid;
        // canonical representative -> gathered equivalence-class size
        std::unordered_map<EntityKey, int, EntityKey::Hash> rep_count;

        std::vector<EquivClass> classes;
    };

    struct EdgeTopology
    {
        // local edge representative -> canonical quotient edge key
        std::unordered_map<EntityKey, EdgeKey, EntityKey::Hash> local_to_qkey;
        // local edge representative -> sign relative to qkey canonical orientation
        std::unordered_map<EntityKey, int8_t, EntityKey::Hash> local_to_qsign;

        // Final topology view: qkey -> gathered representatives participating
        // in owner/alias management.  Members may belong to ranks other than
        // the current rank; construction-time local/candidate collections live
        // in EdgeBuildScratch, not here.
        std::unordered_map<EdgeKey, std::vector<EntityKey>, EdgeKey::Hash> qkey_to_members;
        std::unordered_map<EdgeKey, EntityKey, EdgeKey::Hash> qkey_to_owner;
        std::unordered_map<EntityKey, bool, EntityKey::Hash> local_is_owner;

        std::unordered_map<EntityKey, int, EntityKey::Hash> owner_to_gid;
        std::unordered_map<int, EntityKey> gid_to_owner;
        int n_local_owner = 0;
        int n_global_owner = 0;
        int owner_gid_begin = 0;
        int owner_gid_end = 0; // half-open: [begin, end)

        std::unordered_map<EdgeKey, int, EdgeKey::Hash> qkey_to_qid;
        std::vector<EquivClass> classes;
    };

    struct FaceTopology
    {
        // local face representative -> canonical quotient face key
        std::unordered_map<EntityKey, FaceKey, EntityKey::Hash> local_to_qkey;
        // local face representative -> sign relative to qkey canonical orientation
        std::unordered_map<EntityKey, int8_t, EntityKey::Hash> local_to_qsign;

        // Final topology view: qkey -> gathered representatives participating
        // in owner/alias management.  Members may belong to ranks other than
        // the current rank; construction-time local/candidate collections live
        // in FaceBuildScratch, not here.
        std::unordered_map<FaceKey, std::vector<EntityKey>, FaceKey::Hash> qkey_to_members;
        std::unordered_map<FaceKey, EntityKey, FaceKey::Hash> qkey_to_owner;
        std::unordered_map<EntityKey, bool, EntityKey::Hash> local_is_owner;

        std::unordered_map<EntityKey, int, EntityKey::Hash> owner_to_gid;
        std::unordered_map<int, EntityKey> gid_to_owner;
        int n_local_owner = 0;
        int n_global_owner = 0;
        int owner_gid_begin = 0;
        int owner_gid_end = 0; // half-open: [begin, end)

        std::unordered_map<FaceKey, int, FaceKey::Hash> qkey_to_qid;
        std::vector<EquivClass> classes;
    };

    struct CellTopology
    {
        // local cell representative -> dimension-local quotient id
        std::unordered_map<EntityKey, int, EntityKey::Hash> local_to_qid;
    };

    struct Topology
    {
        std::vector<InterfacePatch> inner_patches;
        std::vector<InterfacePatch> parallel_patches;
        std::vector<PhysicalPatch> physical_patches;
        std::vector<EdgePatch> inner_edge_patches;
        std::vector<EdgePatch> parallel_edge_patches;
        std::vector<EdgePatch> physical_edge_patches;
        std::vector<VertexPatch> inner_vertex_patches;
        std::vector<VertexPatch> parallel_vertex_patches;
        std::vector<VertexPatch> physical_vertex_patches;

        NodeTopology nodes;
        EdgeTopology edges;
        FaceTopology faces;
        CellTopology cells;

        bool has_node_equiv() const { return !nodes.classes.empty(); }
        bool has_edge_equiv() const
        {
            return !edges.classes.empty() || !edges.qkey_to_owner.empty();
        }
        bool has_face_equiv() const { return !faces.classes.empty(); }

        const std::vector<EquivClass> &classes(EntityDim dim) const;
        const std::vector<EquivClass> &edge_classes_view() const { return edges.classes; }
        const std::vector<EquivClass> &face_classes_view() const { return faces.classes; }

        EntityId id_of(const EntityKey &key) const;
        EntityKey owner_of(const EntityKey &key) const;
        int sign_to_owner(const EntityKey &key) const;
        bool is_owner(const EntityKey &key) const;
        EdgeKey edge_qkey(const EntityKey &edge) const;
        FaceKey face_qkey(const EntityKey &face) const;
        int edge_qsign(const EntityKey &edge) const;
        int face_qsign(const EntityKey &face) const;

        std::string dump_node(const EntityKey &node) const;
        std::string dump_edge(const EntityKey &edge) const;
        std::string dump_face(const EntityKey &face) const;
        std::string dump_edge_class(const EdgeKey &qkey) const;
        std::string dump_face_class(const FaceKey &qkey) const;

    };

    // ============================================================
    // small helpers
    // ============================================================

    // local edge endpoints
    // dir=1 (Xi)   : (i,j,k) -> (i+1,j,k)
    // dir=2 (Eta)  : (i,j,k) -> (i,j+1,k)
    // dir=3 (Zeta) : (i,j,k) -> (i,j,k+1)
    std::pair<EntityKey, EntityKey> endpoints(const EntityKey &e);

    // local face corners in oriented local order.
    // dir=1 (FaceXi) : (i,j,k), (i,j+1,k), (i,j,k+1), (i,j+1,k+1)
    // dir=2 (FaceEt) : (i,j,k), (i+1,j,k), (i,j,k+1), (i+1,j,k+1)
    // dir=3 (FaceZe) : (i,j,k), (i+1,j,k), (i,j+1,k), (i+1,j+1,k)
    std::array<EntityKey, 4> corners(const EntityKey &f);

    // build canonical EdgeKey from a local edge using nodes.local_to_rep
    // returns sign_to_canonical:
    //   +1 if local edge direction matches key.a -> key.b
    //   -1 otherwise
    EdgeKey make_edge_key(
        const EntityKey &e,
        const std::unordered_map<EntityKey, EntityKey, EntityKey::Hash> &local_to_rep,
        int8_t &sign_to_canonical);

    // build canonical FaceKey from a local face using nodes.local_to_rep
    // returns sign_to_canonical based on local corner ordering parity against
    // the sorted canonical corner ordering. Degenerate 2D faces use +1.
    FaceKey make_face_key(
        const EntityKey &f,
        const std::unordered_map<EntityKey, EntityKey, EntityKey::Hash> &local_to_rep,
        int8_t &sign_to_canonical);

    // Compares global edge boundary stencils induced by all members of each
    // FaceKey after normalization by their current faces.local_to_qsign value.  This is
    // a validation hook for the existing sorted-corner orientation rule; it
    // does not change that rule.
    bool validate_face_orientation_stencils(const Topology &topology,
                                            std::ostream &diagnostics);

} // namespace TOPO
