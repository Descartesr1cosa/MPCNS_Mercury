#include "2_topology/Topology.h"

#include "0_basic/Error.h"

#include <algorithm>
#include <sstream>
#include <stdexcept>

namespace TOPO
{
    StaggerLocation stagger_location(const EntityKey &entity)
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

    std::pair<EntityKey, EntityKey> endpoints(const EntityKey &e)
    {
        EntityKey n0 = make_node(e.rank, e.block, e.i, e.j, e.k);
        EntityKey n1 = n0;

        switch (e.axis)
        {
        case EntityAxis::Xi:
            ++n1.i;
            break;
        case EntityAxis::Eta:
            ++n1.j;
            break;
        case EntityAxis::Zeta:
            ++n1.k;
            break;
        default:
            throw std::runtime_error("TOPO::endpoints: invalid edge axis");
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
            break;
        case EntityAxis::Eta:
            ++n1.i;
            ++n2.k;
            ++n3.i;
            ++n3.k;
            break;
        case EntityAxis::Zeta:
            ++n1.i;
            ++n2.j;
            ++n3.i;
            ++n3.j;
            break;
        default:
            throw std::runtime_error("TOPO::corners: invalid face axis");
        }
        return {n0, n1, n2, n3};
    }

    EdgeKey make_edge_key(
        const EntityKey &e,
        const std::unordered_map<EntityKey, EntityKey, EntityKey::Hash> &local_to_rep,
        int8_t &sign_to_canonical)
    {
        const auto [ln0, ln1] = endpoints(e);

        auto it0 = local_to_rep.find(ln0);
        auto it1 = local_to_rep.find(ln1);
        if (it0 == local_to_rep.end() || it1 == local_to_rep.end())
            throw std::runtime_error("TOPO::make_edge_key: endpoint not found in nodes.local_to_rep.");

        const EntityKey &g0 = it0->second;
        const EntityKey &g1 = it1->second;
        if (g0 == g1)
            throw std::runtime_error("TOPO::make_edge_key: two endpoints collapse to the same EntityKey.");

        if (g0 < g1)
        {
            sign_to_canonical = +1;
            return EdgeKey{g0, g1};
        }

        sign_to_canonical = -1;
        return EdgeKey{g1, g0};
    }

    FaceKey make_face_key(
        const EntityKey &f,
        const std::unordered_map<EntityKey, EntityKey, EntityKey::Hash> &local_to_rep,
        int8_t &sign_to_canonical)
    {
        auto local_corners = corners(f);
        std::array<EntityKey, 4> corner_eq{};

        auto lookup = [&](const EntityKey &nid, EntityKey &out) -> bool
        {
            auto it = local_to_rep.find(nid);
            if (it == local_to_rep.end())
                return false;
            out = it->second;
            return true;
        };

        bool all_found = true;
        for (int n = 0; n < 4; ++n)
            all_found = lookup(local_corners[n], corner_eq[n]) && all_found;

        if (!all_found && (f.axis == EntityAxis::Xi || f.axis == EntityAxis::Eta))
        {
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
            throw std::runtime_error("TOPO::make_face_key: corner not found in nodes.local_to_rep.");

        auto sorted = corner_eq;
        std::sort(sorted.begin(), sorted.end());
        const bool degenerate = (sorted[0] == sorted[1]) || (sorted[1] == sorted[2]) || (sorted[2] == sorted[3]);

        sign_to_canonical = +1;
        if (!degenerate)
        {
            int inversions = 0;
            for (int i = 0; i < 4; ++i)
                for (int j = i + 1; j < 4; ++j)
                    if (corner_eq[j] < corner_eq[i])
                        ++inversions;
            sign_to_canonical = (inversions % 2 == 0) ? +1 : -1;
        }

        return FaceKey{sorted[0], sorted[1], sorted[2], sorted[3]};
    }

    EntityId Topology::id_of(const EntityKey &key) const
    {
        switch (key.dim)
        {
        case EntityDim::Node:
        {
            const auto eq_it = nodes.local_to_rep.find(key);
            if (eq_it == nodes.local_to_rep.end())
                throw std::runtime_error("Topology::id_of: node is not present in nodes.local_to_rep.");
            const auto id_it = nodes.rep_to_qid.find(eq_it->second);
            if (id_it == nodes.rep_to_qid.end())
                throw std::runtime_error("Topology::id_of: node quotient id is not built.");
            return EntityId{EntityDim::Node, id_it->second};
        }
        case EntityDim::Edge:
        {
            const auto key_it = edges.local_to_qkey.find(key);
            if (key_it == edges.local_to_qkey.end())
                throw std::runtime_error("Topology::id_of: edge is not present in edges.local_to_qkey.");
            const auto id_it = edges.qkey_to_qid.find(key_it->second);
            if (id_it == edges.qkey_to_qid.end())
                throw std::runtime_error("Topology::id_of: edge quotient id is not built.");
            return EntityId{EntityDim::Edge, id_it->second};
        }
        case EntityDim::Face:
        {
            const auto key_it = faces.local_to_qkey.find(key);
            if (key_it == faces.local_to_qkey.end())
                throw std::runtime_error("Topology::id_of: face is not present in faces.local_to_qkey.");
            const auto id_it = faces.qkey_to_qid.find(key_it->second);
            if (id_it == faces.qkey_to_qid.end())
                throw std::runtime_error("Topology::id_of: face quotient id is not built.");
            return EntityId{EntityDim::Face, id_it->second};
        }
        case EntityDim::Cell:
        {
            const auto id_it = cells.local_to_qid.find(key);
            if (id_it == cells.local_to_qid.end())
                throw std::runtime_error("Topology::id_of: cell quotient id is not built.");
            return EntityId{EntityDim::Cell, id_it->second};
        }
        }
        throw std::runtime_error("Topology::id_of: invalid entity dimension.");
    }

    EntityKey Topology::owner_of(const EntityKey &key) const
    {
        switch (key.dim)
        {
        case EntityDim::Node:
        {
            const auto it = nodes.local_to_rep.find(key);
            if (it == nodes.local_to_rep.end())
                throw std::runtime_error("Topology::owner_of: node is not present in nodes.local_to_rep.");
            return make_node(it->second.rank, it->second.block, it->second.i, it->second.j, it->second.k);
        }
        case EntityDim::Edge:
        {
            const auto key_it = edges.local_to_qkey.find(key);
            if (key_it == edges.local_to_qkey.end())
                throw std::runtime_error("Topology::owner_of: edge is not present in edges.local_to_qkey.");
            const auto owner_it = edges.qkey_to_owner.find(key_it->second);
            return owner_it == edges.qkey_to_owner.end() ? key : owner_it->second;
        }
        case EntityDim::Face:
        {
            const auto key_it = faces.local_to_qkey.find(key);
            if (key_it == faces.local_to_qkey.end())
                throw std::runtime_error("Topology::owner_of: face is not present in faces.local_to_qkey.");
            const auto owner_it = faces.qkey_to_owner.find(key_it->second);
            return owner_it == faces.qkey_to_owner.end() ? key : owner_it->second;
        }
        case EntityDim::Cell:
            if (cells.local_to_qid.find(key) == cells.local_to_qid.end())
                throw std::runtime_error("Topology::owner_of: cell quotient id is not built.");
            return key;
        }
        throw std::runtime_error("Topology::owner_of: invalid entity dimension.");
    }

    int Topology::sign_to_owner(const EntityKey &key) const
    {
        switch (key.dim)
        {
        case EntityDim::Node:
            (void)owner_of(key);
            return +1;
        case EntityDim::Edge:
        {
            const auto sign_it = edges.local_to_qsign.find(key);
            if (sign_it == edges.local_to_qsign.end())
                throw std::runtime_error("Topology::sign_to_owner: edge sign is not present.");
            const auto key_it = edges.local_to_qkey.find(key);
            if (key_it == edges.local_to_qkey.end())
                throw std::runtime_error("Topology::sign_to_owner: edge key is not present.");
            const auto owner_it = edges.qkey_to_owner.find(key_it->second);
            if (owner_it == edges.qkey_to_owner.end())
                return +1;
            const auto owner_sign_it = edges.local_to_qsign.find(owner_it->second);
            if (owner_sign_it == edges.local_to_qsign.end())
                throw std::runtime_error("Topology::sign_to_owner: owner edge sign is not present.");
            return static_cast<int>(sign_it->second) * static_cast<int>(owner_sign_it->second);
        }
        case EntityDim::Face:
        {
            const auto sign_it = faces.local_to_qsign.find(key);
            if (sign_it == faces.local_to_qsign.end())
                throw std::runtime_error("Topology::sign_to_owner: face sign is not present.");
            const auto key_it = faces.local_to_qkey.find(key);
            if (key_it == faces.local_to_qkey.end())
                throw std::runtime_error("Topology::sign_to_owner: face key is not present.");
            const auto owner_it = faces.qkey_to_owner.find(key_it->second);
            if (owner_it == faces.qkey_to_owner.end())
                return +1;
            const auto owner_sign_it = faces.local_to_qsign.find(owner_it->second);
            if (owner_sign_it == faces.local_to_qsign.end())
                throw std::runtime_error("Topology::sign_to_owner: owner face sign is not present.");
            return static_cast<int>(sign_it->second) * static_cast<int>(owner_sign_it->second);
        }
        case EntityDim::Cell:
            (void)owner_of(key);
            return +1;
        }
        throw std::runtime_error("Topology::sign_to_owner: invalid entity dimension.");
    }

    bool Topology::is_owner(const EntityKey &key) const
    {
        return owner_of(key) == key;
    }

    EdgeKey Topology::edge_qkey(const EntityKey &edge) const
    {
        if (edge.dim != EntityDim::Edge)
            throw std::invalid_argument("Topology::edge_qkey: entity is not an edge.");
        const auto it = edges.local_to_qkey.find(edge);
        if (it == edges.local_to_qkey.end())
            throw std::runtime_error("Topology::edge_qkey: edge is not present in edges.local_to_qkey.");
        return it->second;
    }

    FaceKey Topology::face_qkey(const EntityKey &face) const
    {
        if (face.dim != EntityDim::Face)
            throw std::invalid_argument("Topology::face_qkey: entity is not a face.");
        const auto it = faces.local_to_qkey.find(face);
        if (it == faces.local_to_qkey.end())
            throw std::runtime_error("Topology::face_qkey: face is not present in faces.local_to_qkey.");
        return it->second;
    }

    int Topology::edge_qsign(const EntityKey &edge) const
    {
        if (edge.dim != EntityDim::Edge)
            throw std::invalid_argument("Topology::edge_qsign: entity is not an edge.");
        const auto it = edges.local_to_qsign.find(edge);
        if (it == edges.local_to_qsign.end())
            throw std::runtime_error("Topology::edge_qsign: edge sign is not present.");
        return static_cast<int>(it->second);
    }

    int Topology::face_qsign(const EntityKey &face) const
    {
        if (face.dim != EntityDim::Face)
            throw std::invalid_argument("Topology::face_qsign: entity is not a face.");
        const auto it = faces.local_to_qsign.find(face);
        if (it == faces.local_to_qsign.end())
            throw std::runtime_error("Topology::face_qsign: face sign is not present.");
        return static_cast<int>(it->second);
    }

    const std::vector<EquivClass> &Topology::classes(EntityDim dim) const
    {
        switch (dim)
        {
        case EntityDim::Node:
            return nodes.classes;
        case EntityDim::Edge:
            return edges.classes;
        case EntityDim::Face:
            return faces.classes;
        case EntityDim::Cell:
            break;
        }

        ERROR::Abort("Topology::classes: equivalence classes are unavailable for this entity dimension");
        return nodes.classes;
    }
} // namespace TOPO
