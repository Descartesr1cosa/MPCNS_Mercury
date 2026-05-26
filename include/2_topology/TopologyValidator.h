#pragma once

#include <cstdlib>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "2_topology/GlobalIncidence.h"
#include "2_topology/LocalIncidence.h"
#include "2_topology/Topology.h"

namespace TOPO_VALIDATOR
{
    struct ValidationReport
    {
        std::vector<std::string> errors;

        bool ok() const { return errors.empty(); }

        void add(const std::string &message)
        {
            errors.push_back(message);
        }

        void merge(const ValidationReport &other)
        {
            errors.insert(errors.end(), other.errors.begin(), other.errors.end());
        }
    };

    inline const char *entity_dim_name(TOPO::EntityDim dim)
    {
        switch (dim)
        {
        case TOPO::EntityDim::Node:
            return "Node";
        case TOPO::EntityDim::Edge:
            return "Edge";
        case TOPO::EntityDim::Face:
            return "Face";
        case TOPO::EntityDim::Cell:
            return "Cell";
        }
        return "InvalidDim";
    }

    inline const char *entity_axis_name(TOPO::EntityAxis axis)
    {
        switch (axis)
        {
        case TOPO::EntityAxis::None:
            return "None";
        case TOPO::EntityAxis::Xi:
            return "Xi";
        case TOPO::EntityAxis::Eta:
            return "Eta";
        case TOPO::EntityAxis::Zeta:
            return "Zeta";
        }
        return "InvalidAxis";
    }

    inline std::string entity_where(const TOPO::EntityKey &key)
    {
        std::ostringstream oss;
        oss << "rank=" << key.rank << " block=" << key.block
            << " dim=" << entity_dim_name(key.dim)
            << " axis=" << entity_axis_name(key.axis)
            << " ijk=(" << key.i << "," << key.j << "," << key.k << ")";
        return oss.str();
    }

    inline std::string box_string(const Box3 &box)
    {
        std::ostringstream oss;
        oss << "[(" << box.lo.i << "," << box.lo.j << "," << box.lo.k << "),("
            << box.hi.i << "," << box.hi.j << "," << box.hi.k << "))";
        return oss.str();
    }

    inline bool is_valid_half_open_node_box(const Box3 &box)
    {
        return box.lo.i >= 0 && box.lo.j >= 0 && box.lo.k >= 0 &&
               box.lo.i < box.hi.i && box.lo.j < box.hi.j && box.lo.k < box.hi.k;
    }

    inline bool valid_direction(int direction)
    {
        const int axis = std::abs(direction);
        return axis >= 1 && axis <= 3;
    }

    inline bool valid_transform(const TOPO::IndexTransform &transform, int dimension)
    {
        if (dimension < 1 || dimension > 3)
            return false;

        bool used[3] = {false, false, false};
        for (int axis = 0; axis < 3; ++axis)
        {
            const int mapped = transform.perm[axis];
            if (mapped < 0 || mapped > 2 || used[mapped])
                return false;
            used[mapped] = true;

            const int sign = transform.sign[axis];
            if (axis < dimension)
            {
                if (sign != -1 && sign != +1)
                    return false;
            }
            else if (sign != -1 && sign != 0 && sign != +1)
            {
                return false;
            }
        }
        return true;
    }

    inline Box3 mapped_box(const Box3 &box, const TOPO::IndexTransform &transform)
    {
        const int lo[3] = {box.lo.i, box.lo.j, box.lo.k};
        const int hi[3] = {box.hi.i - 1, box.hi.j - 1, box.hi.k - 1};
        int mapped_lo[3] = {0, 0, 0};
        int mapped_hi[3] = {0, 0, 0};
        const int offset[3] = {transform.offset.i, transform.offset.j, transform.offset.k};

        for (int axis = 0; axis < 3; ++axis)
        {
            const int target = transform.perm[axis];
            const int first = transform.sign[axis] * lo[axis] + offset[axis];
            const int last = transform.sign[axis] * hi[axis] + offset[axis];
            mapped_lo[target] = first < last ? first : last;
            mapped_hi[target] = (first < last ? last : first) + 1;
        }
        return Box3{{mapped_lo[0], mapped_lo[1], mapped_lo[2]},
                    {mapped_hi[0], mapped_hi[1], mapped_hi[2]}};
    }

    inline bool same_box(const Box3 &lhs, const Box3 &rhs)
    {
        return lhs.lo.i == rhs.lo.i && lhs.lo.j == rhs.lo.j && lhs.lo.k == rhs.lo.k &&
               lhs.hi.i == rhs.hi.i && lhs.hi.j == rhs.hi.j && lhs.hi.k == rhs.hi.k;
    }

    inline void check_interface_patch(ValidationReport &report,
                                      const TOPO::InterfacePatch &patch,
                                      int dimension,
                                      const char *list_name,
                                      std::size_t index)
    {
        std::ostringstream where;
        where << list_name << "[" << index << "] rank=" << patch.this_rank
              << " block=" << patch.this_block << " direction=" << patch.direction
              << " nb_rank=" << patch.nb_rank << " nb_block=" << patch.nb_block
              << " nb_direction=" << patch.nb_direction;

        if (!is_valid_half_open_node_box(patch.this_box_node))
            report.add(where.str() + ": invalid this_box_node=" + box_string(patch.this_box_node));
        if (!is_valid_half_open_node_box(patch.nb_box_node))
            report.add(where.str() + ": invalid nb_box_node=" + box_string(patch.nb_box_node));
        if (!valid_direction(patch.direction) || !valid_direction(patch.nb_direction))
            report.add(where.str() + ": face directions must be nonzero +/-1, +/-2, or +/-3");
        if (!valid_transform(patch.trans, dimension))
        {
            report.add(where.str() + ": invalid IndexTransform");
        }
        else if (is_valid_half_open_node_box(patch.this_box_node) &&
                 is_valid_half_open_node_box(patch.nb_box_node))
        {
            const Box3 transformed = mapped_box(patch.this_box_node, patch.trans);
            if (!same_box(transformed, patch.nb_box_node))
            {
                report.add(where.str() + ": IndexTransform maps this_box_node to " +
                           box_string(transformed) + " instead of nb_box_node=" +
                           box_string(patch.nb_box_node));
            }
        }
    }

    inline void check_edge_patch(ValidationReport &report,
                                 const TOPO::EdgePatch &patch,
                                 int dimension,
                                 const char *list_name,
                                 std::size_t index)
    {
        std::ostringstream where;
        where << list_name << "[" << index << "] rank=" << patch.this_rank
              << " block=" << patch.this_block << " directions=("
              << patch.dir1 << "," << patch.dir2 << ")";
        if (!is_valid_half_open_node_box(patch.this_box_node))
            report.add(where.str() + ": invalid this_box_node=" + box_string(patch.this_box_node));
        if (!valid_direction(patch.dir1) || !valid_direction(patch.dir2) ||
            std::abs(patch.dir1) == std::abs(patch.dir2))
            report.add(where.str() + ": edge patch requires two valid distinct face axes");
        if (patch.kind != TOPO::PatchKind::Physical)
        {
            if (!is_valid_half_open_node_box(patch.nb_box_node))
                report.add(where.str() + ": invalid nb_box_node=" + box_string(patch.nb_box_node));
            if (!valid_transform(patch.trans, dimension))
                report.add(where.str() + ": invalid IndexTransform");
            else if (is_valid_half_open_node_box(patch.this_box_node) &&
                     is_valid_half_open_node_box(patch.nb_box_node) &&
                     !same_box(mapped_box(patch.this_box_node, patch.trans), patch.nb_box_node))
                report.add(where.str() + ": transformed edge node box differs from neighbor box");
        }
    }

    inline void check_physical_patch(ValidationReport &report,
                                     const TOPO::PhysicalPatch &patch,
                                     const char *list_name,
                                     std::size_t index)
    {
        std::ostringstream where;
        where << list_name << "[" << index << "] rank=" << patch.this_rank
              << " block=" << patch.this_block << " direction=" << patch.direction;
        if (!is_valid_half_open_node_box(patch.this_box_node))
            report.add(where.str() + ": invalid this_box_node=" + box_string(patch.this_box_node));
        if (!valid_direction(patch.direction))
            report.add(where.str() + ": physical face direction must be nonzero +/-1, +/-2, or +/-3");
    }

    inline void check_vertex_patch(ValidationReport &report,
                                   const TOPO::VertexPatch &patch,
                                   int dimension,
                                   const char *list_name,
                                   std::size_t index)
    {
        std::ostringstream where;
        where << list_name << "[" << index << "] rank=" << patch.this_rank
              << " block=" << patch.this_block << " directions=("
              << patch.dir1 << "," << patch.dir2 << "," << patch.dir3 << ")";
        if (!is_valid_half_open_node_box(patch.this_box_node))
            report.add(where.str() + ": invalid this_box_node=" + box_string(patch.this_box_node));
        if (!valid_direction(patch.dir1) || !valid_direction(patch.dir2) ||
            !valid_direction(patch.dir3) ||
            std::abs(patch.dir1) == std::abs(patch.dir2) ||
            std::abs(patch.dir1) == std::abs(patch.dir3) ||
            std::abs(patch.dir2) == std::abs(patch.dir3))
            report.add(where.str() + ": vertex patch requires three valid distinct face axes");
        if (patch.kind != TOPO::PatchKind::Physical)
        {
            if (!is_valid_half_open_node_box(patch.nb_box_node))
                report.add(where.str() + ": invalid nb_box_node=" + box_string(patch.nb_box_node));
            if (!valid_transform(patch.trans, dimension))
                report.add(where.str() + ": invalid IndexTransform");
            else if (is_valid_half_open_node_box(patch.this_box_node) &&
                     is_valid_half_open_node_box(patch.nb_box_node) &&
                     !same_box(mapped_box(patch.this_box_node, patch.trans), patch.nb_box_node))
                report.add(where.str() + ": transformed vertex node box differs from neighbor box");
        }
    }

    inline ValidationReport check_patch_connectivity(const TOPO::Topology &topology,
                                                      int dimension)
    {
        ValidationReport report;
        for (std::size_t n = 0; n < topology.inner_patches.size(); ++n)
            check_interface_patch(report, topology.inner_patches[n], dimension, "inner_patches", n);
        for (std::size_t n = 0; n < topology.parallel_patches.size(); ++n)
            check_interface_patch(report, topology.parallel_patches[n], dimension, "parallel_patches", n);
        for (std::size_t n = 0; n < topology.physical_patches.size(); ++n)
            check_physical_patch(report, topology.physical_patches[n], "physical_patches", n);
        for (std::size_t n = 0; n < topology.inner_edge_patches.size(); ++n)
            check_edge_patch(report, topology.inner_edge_patches[n], dimension, "inner_edge_patches", n);
        for (std::size_t n = 0; n < topology.parallel_edge_patches.size(); ++n)
            check_edge_patch(report, topology.parallel_edge_patches[n], dimension, "parallel_edge_patches", n);
        for (std::size_t n = 0; n < topology.physical_edge_patches.size(); ++n)
            check_edge_patch(report, topology.physical_edge_patches[n], dimension, "physical_edge_patches", n);
        for (std::size_t n = 0; n < topology.inner_vertex_patches.size(); ++n)
            check_vertex_patch(report, topology.inner_vertex_patches[n], dimension, "inner_vertex_patches", n);
        for (std::size_t n = 0; n < topology.parallel_vertex_patches.size(); ++n)
            check_vertex_patch(report, topology.parallel_vertex_patches[n], dimension, "parallel_vertex_patches", n);
        for (std::size_t n = 0; n < topology.physical_vertex_patches.size(); ++n)
            check_vertex_patch(report, topology.physical_vertex_patches[n], dimension, "physical_vertex_patches", n);
        return report;
    }

    inline ValidationReport check_entity(const TOPO::EntityKey &key,
                                         const TOPO::BlockTopoSize &size)
    {
        ValidationReport report;
        if (!TOPO::valid_axis_for_dim(key.dim, key.axis))
            report.add(entity_where(key) + ": axis is invalid for entity dimension");
        if (!TOPO::is_valid_entity_key_basic(key))
            report.add(entity_where(key) + ": invalid basic local entity key; ghost/negative indices are not topology entities");
        else if (!TOPO::is_valid_entity_key(key, size))
            report.add(entity_where(key) + ": index outside node-based block topology range; ghost indices are not topology entities");
        return report;
    }

    inline ValidationReport check_local_face_incidence(const TOPO::EntityKey &face)
    {
        ValidationReport report;
        try
        {
            if (!TOPO::check_boundary_boundary_face_zero(face))
                report.add(entity_where(face) + ": local boundary(boundary(face)) is nonzero");
        }
        catch (const std::exception &error)
        {
            report.add(entity_where(face) + ": cannot validate local face incidence: " + error.what());
        }
        return report;
    }

    inline ValidationReport check_local_cell_incidence(const TOPO::EntityKey &cell)
    {
        ValidationReport report;
        try
        {
            if (!TOPO::check_boundary_boundary_cell_zero(cell))
                report.add(entity_where(cell) + ": local boundary(boundary(cell)) is nonzero");
        }
        catch (const std::exception &error)
        {
            report.add(entity_where(cell) + ": cannot validate local cell incidence: " + error.what());
        }
        return report;
    }

    inline std::string edge_member_string(const TOPO::EntityKey &member)
    {
        std::ostringstream oss;
        oss << "rank=" << member.rank << " block=" << member.block
            << " dim=Edge axis=" << TOPO::axis_number(member.axis) << " ijk=("
            << member.i << "," << member.j << "," << member.k << ")";
        return oss.str();
    }

    inline std::string face_member_string(const TOPO::EntityKey &member)
    {
        std::ostringstream oss;
        oss << "rank=" << member.rank << " block=" << member.block
            << " dim=Face axis=" << TOPO::axis_number(member.axis) << " ijk=("
            << member.i << "," << member.j << "," << member.k << ")";
        return oss.str();
    }

    inline ValidationReport check_equivalence(const TOPO::Topology &equiv)
    {
        ValidationReport report;
        std::map<TOPO::EntityKey, std::set<TOPO::EntityKey>> node_members;
        for (const auto &[node, owner] : equiv.node2eq)
        {
            node_members[owner].insert(node);
            if (equiv.node_eq_to_id.find(owner) == equiv.node_eq_to_id.end())
            {
                report.add("Node equivalence owner rank=" + std::to_string(owner.rank) +
                           " block=" + std::to_string(owner.block) +
                           ": missing quotient EntityId");
            }
        }
        for (const auto &[owner, members] : node_members)
        {
            const TOPO::EntityKey canonical =
                TOPO::make_node(owner.rank, owner.block, owner.i, owner.j, owner.k);
            const auto local_owner = equiv.node2eq.find(canonical);
            if (local_owner != equiv.node2eq.end() && !(local_owner->second == owner))
                report.add("Node canonical owner maps to a different EntityKey at rank=" +
                           std::to_string(owner.rank) + " block=" + std::to_string(owner.block));
            (void)members;
        }

        for (const auto &[edge, sign] : equiv.edge2sign)
        {
            if (sign != +1 && sign != -1)
                report.add(edge_member_string(edge) + ": invalid edge sign=" + std::to_string(sign));
        }
        for (const auto &[face, sign] : equiv.face2sign)
        {
            if (sign != +1 && sign != -1)
                report.add(face_member_string(face) + ": invalid face sign=" + std::to_string(sign));
        }

        for (const auto &[key, members] : equiv.edge_members)
        {
            const auto owner_it = equiv.edge_owner.find(key);
            if (owner_it == equiv.edge_owner.end())
            {
                report.add("Edge equivalence class has no owner");
                continue;
            }
            std::set<TOPO::EntityKey> unique;
            int owner_count = 0;
            for (const TOPO::EntityKey &member : members)
            {
                if (!unique.insert(member).second)
                    report.add(edge_member_string(member) + ": duplicate edge class member");
                if (member == owner_it->second)
                    ++owner_count;
                const auto member_key = equiv.edge2key.find(member);
                if (member_key == equiv.edge2key.end() || !(member_key->second == key))
                    report.add(edge_member_string(member) + ": member does not map back to its EdgeKey");
            }
            if (owner_count != 1)
                report.add(edge_member_string(owner_it->second) +
                           ": edge class owner occurrence count=" + std::to_string(owner_count));
        }

        for (const auto &[key, members] : equiv.face_members)
        {
            const auto owner_it = equiv.face_owner.find(key);
            if (owner_it == equiv.face_owner.end())
            {
                report.add("Face equivalence class has no owner");
                continue;
            }
            std::set<TOPO::EntityKey> unique;
            int owner_count = 0;
            for (const TOPO::EntityKey &member : members)
            {
                if (!unique.insert(member).second)
                    report.add(face_member_string(member) + ": duplicate face class member");
                if (member == owner_it->second)
                    ++owner_count;
                const auto member_key = equiv.face2key.find(member);
                if (member_key == equiv.face2key.end() || !(member_key->second == key))
                    report.add(face_member_string(member) + ": member does not map back to its FaceKey");
            }
            if (owner_count != 1)
                report.add(face_member_string(owner_it->second) +
                           ": face class owner occurrence count=" + std::to_string(owner_count));
        }

        const auto check_declared_classes =
            [&report](const std::vector<TOPO::EquivClass> &classes, const char *kind)
        {
            for (std::size_t n = 0; n < classes.size(); ++n)
            {
                const TOPO::EquivClass &cls = classes[n];
                using MemberKey = std::tuple<int, int, int, int, int, int>;
                const auto member_key = [](const TOPO::EquivMember &member)
                {
                    return MemberKey{member.rank, member.block,
                                     static_cast<int>(member.location),
                                     member.i, member.j, member.k};
                };

                std::set<MemberKey> unique;
                int flagged_owners = 0;
                bool owner_is_member = false;
                for (const TOPO::EquivMember &member : cls.members)
                {
                    if (!unique.insert(member_key(member)).second)
                    {
                        report.add(std::string(kind) + " class[" + std::to_string(n) +
                                   "]: duplicate member rank=" + std::to_string(member.rank) +
                                   " block=" + std::to_string(member.block) +
                                   " ijk=(" + std::to_string(member.i) + "," +
                                   std::to_string(member.j) + "," + std::to_string(member.k) + ")");
                    }
                    if (member.is_owner)
                        ++flagged_owners;
                    if (member_key(member) == member_key(cls.owner))
                        owner_is_member = true;
                }
                if (flagged_owners != 1)
                    report.add(std::string(kind) + " class[" + std::to_string(n) +
                               "]: owner flag count=" + std::to_string(flagged_owners));
                if (!owner_is_member)
                    report.add(std::string(kind) + " class[" + std::to_string(n) +
                               "]: declared owner is not present in members");
            }
        };
        check_declared_classes(equiv.node_classes, "Node");
        check_declared_classes(equiv.edge_classes_general, "Edge");
        check_declared_classes(equiv.face_classes, "Face");
        return report;
    }

    inline ValidationReport check_face_orientation(const TOPO::Topology &equiv)
    {
        ValidationReport report;
        std::ostringstream diagnostics;
        if (!TOPO::validate_face_orientation_stencils(equiv, diagnostics))
            report.add(diagnostics.str());
        return report;
    }

    inline ValidationReport check_global_incidence(const TOPO::GlobalIncidence &incidence,
                                                   bool require_cell_incidence = false)
    {
        ValidationReport report;
        if (!incidence.check_d1_d0_zero())
            report.add("Global incidence validation failed: D1 * D0 is nonzero");
        if (require_cell_incidence)
        {
            try
            {
                if (!incidence.check_d2_d1_zero())
                    report.add("Global incidence validation failed: D2 * D1 is nonzero");
            }
            catch (const std::exception &error)
            {
                report.add(std::string("Global incidence D2 * D1 unavailable: ") + error.what());
            }
        }
        return report;
    }
} // namespace TOPO_VALIDATOR
