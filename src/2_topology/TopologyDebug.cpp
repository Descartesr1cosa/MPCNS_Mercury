
#include "2_topology/TopologyDebug.h"
#include "2_topology/LocalIncidence.h"

#include "0_basic/BoxOps.h"
#include "0_basic/Direction.h"
#include "0_basic/Error.h"
#include "0_basic/LayoutTraits.h"
#include "1_grid/BlockTraits.h"
#include "1_grid/1_MPCNS_Grid.h"

#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>

namespace
{
    Box3 node_domain_of_block(const Block &blk)
    {
        const Int3 ncells = GRID_TRAITS::cell_counts(blk);
        const Int3 nnodes = LAYOUT::node_size_from_cells(ncells);

        return Box3{{0, 0, 0}, {nnodes.i, nnodes.j, nnodes.k}};
    }

    std::string patch_where(const char *kind, int block, const Box3 &box)
    {
        std::ostringstream oss;
        oss << kind << " block=" << block
            << " box=" << BOX::to_string(box);
        return oss.str();
    }

    void validate_block_id(int ib, int nblock, const char *where)
    {
        if (ib < 0 || ib >= nblock)
        {
            std::ostringstream oss;
            oss << where << ": invalid block id " << ib
                << ", nblock=" << nblock;
            ERROR::Abort(oss.str());
        }
    }

    void validate_node_box_inside_block(const Block &blk,
                                        const Box3 &box,
                                        const char *where)
    {
        const Box3 domain = node_domain_of_block(blk);
        if (!BOX::contains(domain, box))
        {
            std::ostringstream oss;
            oss << where
                << ": node box outside block node domain. "
                << "box=" << BOX::to_string(box)
                << " domain=" << BOX::to_string(domain);
            ERROR::Abort(oss.str());
        }

        BOX::assert_nonempty(box, where);
    }

    void validate_transform(const TOPO::IndexTransform &tr,
                            int dimension,
                            const char *where)
    {
        bool used[3] = {false, false, false};

        for (int d = 0; d < 3; ++d)
        {
            if (tr.perm[d] < 0 || tr.perm[d] > 2)
            {
                std::ostringstream oss;
                oss << where << ": invalid perm[" << d << "]=" << tr.perm[d];
                ERROR::Abort(oss.str());
            }

            if (used[tr.perm[d]])
            {
                std::ostringstream oss;
                oss << where << ": duplicated perm value " << tr.perm[d];
                ERROR::Abort(oss.str());
            }

            used[tr.perm[d]] = true;
        }

        for (int d = 0; d < 3; ++d)
        {
            const int s = tr.sign[d];

            if (d < dimension)
            {
                if (s != -1 && s != +1)
                {
                    std::ostringstream oss;
                    oss << where << ": invalid active sign[" << d << "]=" << s;
                    ERROR::Abort(oss.str());
                }
            }
            else
            {
                // 2D 情况下第三维可能为 0，也允许 ±1。
                if (s != -1 && s != 0 && s != +1)
                {
                    std::ostringstream oss;
                    oss << where << ": invalid inactive sign[" << d << "]=" << s;
                    ERROR::Abort(oss.str());
                }
            }
        }
    }

    void validate_interface_patch(const TOPO::InterfacePatch &p,
                                  Grid &grid,
                                  int dimension,
                                  const char *list_name,
                                  int index)
    {
        std::ostringstream where;
        where << list_name << "[" << index << "]";

        validate_block_id(p.this_block, grid.nblock, where.str().c_str());

        const Block &blk = grid.grids(p.this_block);

        validate_node_box_inside_block(
            blk,
            p.this_box_node,
            patch_where(where.str().c_str(), p.this_block, p.this_box_node).c_str());

        if (!DIR::is_valid(p.direction))
        {
            std::ostringstream oss;
            oss << where.str() << ": invalid direction=" << p.direction;
            ERROR::Abort(oss.str());
        }

        if (!DIR::is_valid(p.nb_direction))
        {
            std::ostringstream oss;
            oss << where.str() << ": invalid nb_direction=" << p.nb_direction;
            ERROR::Abort(oss.str());
        }

        validate_transform(p.trans, dimension, where.str().c_str());
    }

    void validate_physical_patch(const TOPO::PhysicalPatch &p,
                                 Grid &grid,
                                 const char *list_name,
                                 int index)
    {
        std::ostringstream where;
        where << list_name << "[" << index << "]";

        validate_block_id(p.this_block, grid.nblock, where.str().c_str());

        const Block &blk = grid.grids(p.this_block);

        validate_node_box_inside_block(
            blk,
            p.this_box_node,
            patch_where(where.str().c_str(), p.this_block, p.this_box_node).c_str());

        if (!DIR::is_valid(p.direction))
        {
            std::ostringstream oss;
            oss << where.str() << ": invalid direction=" << p.direction;
            ERROR::Abort(oss.str());
        }
    }

    void validate_edge_patch(const TOPO::EdgePatch &p,
                             Grid &grid,
                             const char *list_name,
                             int index)
    {
        std::ostringstream where;
        where << list_name << "[" << index << "]";

        validate_block_id(p.this_block, grid.nblock, where.str().c_str());

        const Block &blk = grid.grids(p.this_block);

        validate_node_box_inside_block(
            blk,
            p.this_box_node,
            patch_where(where.str().c_str(), p.this_block, p.this_box_node).c_str());

        if (!DIR::distinct_axes(p.dir1, p.dir2))
        {
            std::ostringstream oss;
            oss << where.str() << ": edge directions share the same axis: "
                << p.dir1 << ", " << p.dir2;
            ERROR::Abort(oss.str());
        }
    }

    void validate_vertex_patch(const TOPO::VertexPatch &p,
                               Grid &grid,
                               const char *list_name,
                               int index)
    {
        std::ostringstream where;
        where << list_name << "[" << index << "]";

        validate_block_id(p.this_block, grid.nblock, where.str().c_str());

        const Block &blk = grid.grids(p.this_block);

        validate_node_box_inside_block(
            blk,
            p.this_box_node,
            patch_where(where.str().c_str(), p.this_block, p.this_box_node).c_str());

        if (!DIR::distinct_axes(p.dir1, p.dir2, p.dir3))
        {
            std::ostringstream oss;
            oss << where.str() << ": vertex directions are not three distinct axes: "
                << p.dir1 << ", " << p.dir2 << ", " << p.dir3;
            ERROR::Abort(oss.str());
        }
    }
}

namespace TOPO_DEBUG
{
    void dump_topology_summary(const TOPO::Topology &topo, int my_rank)
    {
        if (my_rank != 0)
            return;

        std::cout << "\n========== Topology summary ==========\n";
        std::cout << "faces:\n";
        std::cout << "  inner    = " << topo.inner_patches.size() << "\n";
        std::cout << "  parallel = " << topo.parallel_patches.size() << "\n";
        std::cout << "  physical = " << topo.physical_patches.size() << "\n";

        std::cout << "edges:\n";
        std::cout << "  inner    = " << topo.inner_edge_patches.size() << "\n";
        std::cout << "  parallel = " << topo.parallel_edge_patches.size() << "\n";
        std::cout << "  physical = " << topo.physical_edge_patches.size() << "\n";

        std::cout << "vertices:\n";
        std::cout << "  inner    = " << topo.inner_vertex_patches.size() << "\n";
        std::cout << "  parallel = " << topo.parallel_vertex_patches.size() << "\n";
        std::cout << "  physical = " << topo.physical_vertex_patches.size() << "\n";
        std::cout << "======================================\n\n";
    }

    void validate_topology_or_abort(TOPO::Topology &topo,
                                    Grid &grid,
                                    int my_rank,
                                    int dimension)
    {
        (void)my_rank;

        for (int n = 0; n < static_cast<int>(topo.inner_patches.size()); ++n)
            validate_interface_patch(topo.inner_patches[n], grid, dimension, "inner_patches", n);

        for (int n = 0; n < static_cast<int>(topo.parallel_patches.size()); ++n)
            validate_interface_patch(topo.parallel_patches[n], grid, dimension, "parallel_patches", n);

        for (int n = 0; n < static_cast<int>(topo.physical_patches.size()); ++n)
            validate_physical_patch(topo.physical_patches[n], grid, "physical_patches", n);

        for (int n = 0; n < static_cast<int>(topo.inner_edge_patches.size()); ++n)
            validate_edge_patch(topo.inner_edge_patches[n], grid, "inner_edge_patches", n);

        for (int n = 0; n < static_cast<int>(topo.parallel_edge_patches.size()); ++n)
            validate_edge_patch(topo.parallel_edge_patches[n], grid, "parallel_edge_patches", n);

        for (int n = 0; n < static_cast<int>(topo.physical_edge_patches.size()); ++n)
            validate_edge_patch(topo.physical_edge_patches[n], grid, "physical_edge_patches", n);

        for (int n = 0; n < static_cast<int>(topo.inner_vertex_patches.size()); ++n)
            validate_vertex_patch(topo.inner_vertex_patches[n], grid, "inner_vertex_patches", n);

        for (int n = 0; n < static_cast<int>(topo.parallel_vertex_patches.size()); ++n)
            validate_vertex_patch(topo.parallel_vertex_patches[n], grid, "parallel_vertex_patches", n);

        for (int n = 0; n < static_cast<int>(topo.physical_vertex_patches.size()); ++n)
            validate_vertex_patch(topo.physical_vertex_patches[n], grid, "physical_vertex_patches", n);
    }
}

namespace TOPO
{
namespace
{
    const char *dim_name(EntityDim dim)
    {
        switch (dim)
        {
        case EntityDim::Node:
        return "node";
        case EntityDim::Edge:
        return "edge";
        case EntityDim::Face:
        return "face";
        case EntityDim::Cell:
        return "cell";
        }
        return "invalid";
    }

    
    std::string entity_string(const EntityKey &entity)
    {
        std::ostringstream oss;
        oss << dim_name(entity.dim) << "("
        << "rank=" << entity.rank
        << ",block=" << entity.block
        << ",i=" << entity.i
        << ",j=" << entity.j
        << ",k=" << entity.k
        << ",axis=" << axis_number(entity.axis)
        << ")";
        return oss.str();
    }

    
    std::string edge_key_string(const EdgeKey &key)
    {
        std::ostringstream oss;
        oss << "edge_qkey{a=" << entity_string(key.a)
        << ",b=" << entity_string(key.b) << "}";
        return oss.str();
    }

    
    std::string face_key_string(const FaceKey &key)
    {
        std::ostringstream oss;
        oss << "face_qkey{a=" << entity_string(key.a)
        << ",b=" << entity_string(key.b)
        << ",c=" << entity_string(key.c)
        << ",d=" << entity_string(key.d) << "}";
        return oss.str();
    }

    
    void append_members(std::ostream &os, const std::vector<EntityKey> &members)
    {
        os << "[";
        for (std::size_t n = 0; n < members.size(); ++n)
        {
        if (n != 0)
            os << ",";
        os << entity_string(members[n]);
        }
        os << "]";
    }

    
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

std::string Topology::dump_node(const EntityKey &node) const
    {
        if (node.dim != EntityDim::Node)
            throw std::invalid_argument("Topology::dump_node: entity is not a node.");

        const auto rep_it = nodes.local_to_rep.find(node);
        if (rep_it == nodes.local_to_rep.end())
            throw std::runtime_error("Topology::dump_node: node is not present in nodes.local_to_rep.");

        const EntityKey &rep = rep_it->second;
        const auto qid = id_of(node);
        const auto count_it = nodes.rep_count.find(rep);

        std::ostringstream oss;
        oss << "node local=" << entity_string(node)
            << " rep=" << entity_string(rep)
            << " qid=" << qid.id
            << " owner=" << entity_string(owner_of(node))
            << " sign_to_owner=" << sign_to_owner(node)
            << " rep_count=" << (count_it == nodes.rep_count.end() ? 1 : count_it->second);
        return oss.str();
    }

    
std::string Topology::dump_edge(const EntityKey &edge) const
    {
        const EdgeKey qkey = edge_qkey(edge);
        const auto qid = id_of(edge);
        const EntityKey owner = owner_of(edge);
        const auto gid_it = edges.owner_to_gid.find(owner);
        const auto members_it = edges.qkey_to_members.find(qkey);

        std::ostringstream oss;
        oss << "edge local=" << entity_string(edge)
            << " qkey=" << edge_key_string(qkey)
            << " qsign=" << edge_qsign(edge)
            << " qid=" << qid.id
            << " owner=" << entity_string(owner)
            << " sign_to_owner=" << sign_to_owner(edge)
            << " owner_gid=" << (gid_it == edges.owner_to_gid.end() ? -1 : gid_it->second)
            << " members=";
        if (members_it == edges.qkey_to_members.end())
            oss << "[]";
        else
            append_members(oss, members_it->second);
        return oss.str();
    }

    
std::string Topology::dump_face(const EntityKey &face) const
    {
        const FaceKey qkey = face_qkey(face);
        const auto qid = id_of(face);
        const EntityKey owner = owner_of(face);
        const auto gid_it = faces.owner_to_gid.find(owner);
        const auto members_it = faces.qkey_to_members.find(qkey);

        std::ostringstream oss;
        oss << "face local=" << entity_string(face)
            << " qkey=" << face_key_string(qkey)
            << " qsign=" << face_qsign(face)
            << " qid=" << qid.id
            << " owner=" << entity_string(owner)
            << " sign_to_owner=" << sign_to_owner(face)
            << " owner_gid=" << (gid_it == faces.owner_to_gid.end() ? -1 : gid_it->second)
            << " members=";
        if (members_it == faces.qkey_to_members.end())
            oss << "[]";
        else
            append_members(oss, members_it->second);
        return oss.str();
    }

    
std::string Topology::dump_edge_class(const EdgeKey &qkey) const
    {
        const auto owner_it = edges.qkey_to_owner.find(qkey);
        const auto members_it = edges.qkey_to_members.find(qkey);
        const auto qid_it = edges.qkey_to_qid.find(qkey);
        const int owner_gid = (owner_it == edges.qkey_to_owner.end())
                                  ? -1
                                  : [&]() {
                                        const auto gid_it = edges.owner_to_gid.find(owner_it->second);
                                        return gid_it == edges.owner_to_gid.end() ? -1 : gid_it->second;
                                    }();

        std::ostringstream oss;
        oss << "edge_class qkey=" << edge_key_string(qkey)
            << " qid=" << (qid_it == edges.qkey_to_qid.end() ? -1 : qid_it->second)
            << " owner=" << (owner_it == edges.qkey_to_owner.end() ? "<none>" : entity_string(owner_it->second))
            << " owner_gid=" << owner_gid
            << " members=";
        if (members_it == edges.qkey_to_members.end())
            oss << "[]";
        else
            append_members(oss, members_it->second);
        return oss.str();
    }

    
std::string Topology::dump_face_class(const FaceKey &qkey) const
    {
        const auto owner_it = faces.qkey_to_owner.find(qkey);
        const auto members_it = faces.qkey_to_members.find(qkey);
        const auto qid_it = faces.qkey_to_qid.find(qkey);
        const int owner_gid = (owner_it == faces.qkey_to_owner.end())
                                  ? -1
                                  : [&]() {
                                        const auto gid_it = faces.owner_to_gid.find(owner_it->second);
                                        return gid_it == faces.owner_to_gid.end() ? -1 : gid_it->second;
                                    }();

        std::ostringstream oss;
        oss << "face_class qkey=" << face_key_string(qkey)
            << " qid=" << (qid_it == faces.qkey_to_qid.end() ? -1 : qid_it->second)
            << " owner=" << (owner_it == faces.qkey_to_owner.end() ? "<none>" : entity_string(owner_it->second))
            << " owner_gid=" << owner_gid
            << " members=";
        if (members_it == faces.qkey_to_members.end())
            oss << "[]";
        else
            append_members(oss, members_it->second);
        return oss.str();
    }

    
bool validate_face_orientation_stencils(const Topology &equiv,
                                            std::ostream &diagnostics)
    {
        bool valid = true;
        for (const auto &[key, members] : equiv.faces.qkey_to_members)
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

    
} // namespace TOPO
