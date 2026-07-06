#pragma once
#include <compare>
#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iosfwd>
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

        // local node entity -> canonical node entity
        std::unordered_map<EntityKey, EntityKey, EntityKey::Hash> node2eq;
        // Canonical quotient ids used by the EntityKey facade.  These ids are
        // dimension-local and do not replace owner-sync gids below.
        std::unordered_map<EntityKey, int, EntityKey::Hash> node_eq_to_id;

        // local edge -> canonical edge key
        std::unordered_map<EntityKey, EdgeKey, EntityKey::Hash> edge2key;

        // local edge -> sign to canonical edge direction
        // +1 : local direction == key.a -> key.b
        // -1 : local direction == key.b -> key.a
        std::unordered_map<EntityKey, int8_t, EntityKey::Hash> edge2sign;

        // canonical edge key -> all local members on this rank
        std::unordered_map<EdgeKey, std::vector<EntityKey>, EdgeKey::Hash> edge_members;

        // canonical edge key -> chosen owner rep
        std::unordered_map<EdgeKey, EntityKey, EdgeKey::Hash> edge_owner;

        // local edge -> whether this rep is owner
        std::unordered_map<EntityKey, bool, EntityKey::Hash> edge_is_owner;

        std::unordered_map<EntityKey, int, EntityKey::Hash> edge_owner_gid;
        std::unordered_map<int, EntityKey> gid2edge_owner;
        int n_local_edge_owner = 0;
        int n_global_edge_owner = 0;
        int edge_owner_gid_begin = 0;
        int edge_owner_gid_end = 0; // half-open: [begin, end)
        std::unordered_map<EdgeKey, int, EdgeKey::Hash> edge_key_to_id;

        std::unordered_map<EntityKey, FaceKey, EntityKey::Hash> face2key;
        // Sign of a local face relative to the selected owner face orientation.
        // +1 means the member boundary stencil matches the owner; -1 means the
        // same stencil with all coefficients reversed.
        std::unordered_map<EntityKey, int8_t, EntityKey::Hash> face2sign;
        std::unordered_map<FaceKey, std::vector<EntityKey>, FaceKey::Hash> face_members;
        std::unordered_map<FaceKey, EntityKey, FaceKey::Hash> face_owner;
        std::unordered_map<EntityKey, bool, EntityKey::Hash> face_is_owner;

        std::unordered_map<EntityKey, int, EntityKey::Hash> face_owner_gid;
        std::unordered_map<int, EntityKey> gid2face_owner;
        int n_local_face_owner = 0;
        int n_global_face_owner = 0;
        int face_owner_gid_begin = 0;
        int face_owner_gid_end = 0; // half-open: [begin, end)
        std::unordered_map<FaceKey, int, FaceKey::Hash> face_key_to_id;

        // Cells are volume interiors and are not quotiented across block
        // interfaces, but they still need dimension-local EntityId values so
        // global d2 incidence can close the DEC complex.
        std::unordered_map<EntityKey, int, EntityKey::Hash> cell_to_id;

        std::vector<EquivClass> node_classes;
        std::vector<EquivClass> edge_classes;
        std::vector<EquivClass> face_classes;

        bool has_node_equiv() const { return !node_classes.empty(); }
        bool has_edge_equiv() const
        {
            return !edge_classes.empty() || !edge_owner.empty();
        }
        bool has_face_equiv() const { return !face_classes.empty(); }

        const std::vector<EquivClass> &classes(EntityDim dim) const;

        EntityId id_of(const EntityKey &key) const;
        EntityKey owner_of(const EntityKey &key) const;
        int sign_to_owner(const EntityKey &key) const;
        bool is_owner(const EntityKey &key) const;

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

    // build canonical EdgeKey from a local edge using node2eq
    // returns sign_to_canonical:
    //   +1 if local edge direction matches key.a -> key.b
    //   -1 otherwise
    EdgeKey make_edge_key(
        const EntityKey &e,
        const std::unordered_map<EntityKey, EntityKey, EntityKey::Hash> &node2eq,
        int8_t &sign_to_canonical);

    // build canonical FaceKey from a local face using node2eq
    // returns sign_to_canonical based on local corner ordering parity against
    // the sorted canonical corner ordering. Degenerate 2D faces use +1.
    FaceKey make_face_key(
        const EntityKey &f,
        const std::unordered_map<EntityKey, EntityKey, EntityKey::Hash> &node2eq,
        int8_t &sign_to_canonical);

    // Compares global edge boundary stencils induced by all members of each
    // FaceKey after normalization by their current face2sign value.  This is
    // a validation hook for the existing sorted-corner orientation rule; it
    // does not change that rule.
    bool validate_face_orientation_stencils(const Topology &topology,
                                            std::ostream &diagnostics);

} // namespace TOPO
