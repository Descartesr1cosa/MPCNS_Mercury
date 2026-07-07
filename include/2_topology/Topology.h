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
                auto combine = [&](std::size_t v)
                {
                    h ^= v + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
                };
                combine(EntityKey::Hash{}(x.a));
                combine(EntityKey::Hash{}(x.b));
                return h;
            }
        };
    };

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
                auto combine = [&](std::size_t v)
                {
                    h ^= v + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
                };
                combine(EntityKey::Hash{}(x.a));
                combine(EntityKey::Hash{}(x.b));
                combine(EntityKey::Hash{}(x.c));
                combine(EntityKey::Hash{}(x.d));
                return h;
            }
        };
    };

    struct EquivMember
    {
        EntityKey entity{EntityDim::Cell, 0, 0, 0, 0, 0, EntityAxis::None};

        // member 相对 canonical orientation 的符号。
        int orient_sign = +1;

        bool is_owner = false;
    };

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

    StaggerLocation stagger_location(const EntityKey &entity);
    std::pair<EntityKey, EntityKey> endpoints(const EntityKey &e);
    std::array<EntityKey, 4> corners(const EntityKey &f);

    EdgeKey make_edge_key(
        const EntityKey &e,
        const std::unordered_map<EntityKey, EntityKey, EntityKey::Hash> &local_to_rep,
        int8_t &sign_to_canonical);

    FaceKey make_face_key(
        const EntityKey &f,
        const std::unordered_map<EntityKey, EntityKey, EntityKey::Hash> &local_to_rep,
        int8_t &sign_to_canonical);

    bool validate_face_orientation_stencils(const Topology &topology,
                                            std::ostream &diagnostics);

} // namespace TOPO
