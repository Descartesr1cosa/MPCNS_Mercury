#include "2_topology/TopologyView.h"

#include "2_topology/TopologyOps.h"
#include "0_basic/Error.h"

namespace
{
    TOPO_VIEW::FacePatchView make_face_view(const TOPO::InterfacePatch &p)
    {
        TOPO_VIEW::FacePatchView v;

        v.kind = p.kind;
        v.has_neighbor = true;
        v.is_coupling = p.is_coupling;

        v.this_rank = p.this_rank;
        v.nb_rank = p.nb_rank;

        v.this_block = p.this_block;
        v.nb_block = p.nb_block;

        v.this_block_name = p.this_block_name;
        v.nb_block_name = p.nb_block_name;

        v.this_box_node = p.this_box_node;
        v.nb_box_node = p.nb_box_node;

        v.direction = p.direction;
        v.nb_direction = p.nb_direction;

        v.trans = p.trans;

        v.send_flag = p.send_flag;
        v.recv_flag = p.recv_flag;

        v.bc_id = 0;
        v.bc_name.clear();

        v.interface_patch = &p;
        v.physical_patch = nullptr;

        return v;
    }

    TOPO_VIEW::FacePatchView make_face_view(const TOPO::PhysicalPatch &p)
    {
        TOPO_VIEW::FacePatchView v;

        v.kind = TOPO::PatchKind::Physical;
        v.has_neighbor = false;
        v.is_coupling = false;

        v.this_rank = p.this_rank;
        v.nb_rank = -1;

        v.this_block = p.this_block;
        v.nb_block = -1;

        v.this_block_name = p.this_block_name;
        v.nb_block_name.clear();

        v.this_box_node = p.this_box_node;
        v.nb_box_node = Box3{};

        v.direction = p.direction;
        v.nb_direction = 0;

        v.trans = TOPO::identity_transform();

        v.send_flag = 0;
        v.recv_flag = 0;

        v.bc_id = p.bc_id;
        v.bc_name = p.bc_name;

        v.interface_patch = nullptr;
        v.physical_patch = &p;

        return v;
    }

    void append_interface_faces(std::vector<TOPO_VIEW::FacePatchView> &out,
                                const std::vector<TOPO::InterfacePatch> &patches)
    {
        for (const auto &p : patches)
            out.push_back(make_face_view(p));
    }

    void append_physical_faces(std::vector<TOPO_VIEW::FacePatchView> &out,
                               const std::vector<TOPO::PhysicalPatch> &patches)
    {
        for (const auto &p : patches)
            out.push_back(make_face_view(p));
    }
}

namespace TOPO_VIEW
{
    std::vector<FacePatchView> inner_faces(const TOPO::Topology &topo)
    {
        std::vector<FacePatchView> out;
        out.reserve(topo.inner_patches.size());
        append_interface_faces(out, topo.inner_patches);
        return out;
    }

    std::vector<FacePatchView> parallel_faces(const TOPO::Topology &topo)
    {
        std::vector<FacePatchView> out;
        out.reserve(topo.parallel_patches.size());
        append_interface_faces(out, topo.parallel_patches);
        return out;
    }

    std::vector<FacePatchView> physical_faces(const TOPO::Topology &topo)
    {
        std::vector<FacePatchView> out;
        out.reserve(topo.physical_patches.size());
        append_physical_faces(out, topo.physical_patches);
        return out;
    }

    std::vector<FacePatchView> interface_faces(const TOPO::Topology &topo)
    {
        std::vector<FacePatchView> out;
        out.reserve(topo.inner_patches.size() + topo.parallel_patches.size());
        append_interface_faces(out, topo.inner_patches);
        append_interface_faces(out, topo.parallel_patches);
        return out;
    }

    std::vector<FacePatchView> all_faces(const TOPO::Topology &topo)
    {
        std::vector<FacePatchView> out;
        out.reserve(topo.inner_patches.size() + topo.parallel_patches.size() + topo.physical_patches.size());
        append_interface_faces(out, topo.inner_patches);
        append_interface_faces(out, topo.parallel_patches);
        append_physical_faces(out, topo.physical_patches);
        return out;
    }

    std::vector<FacePatchView> coupling_interfaces(const TOPO::Topology &topo)
    {
        std::vector<FacePatchView> out;
        out.reserve(topo.inner_patches.size() + topo.parallel_patches.size());

        for (const auto &p : topo.inner_patches)
        {
            if (p.is_coupling)
                out.push_back(make_face_view(p));
        }

        for (const auto &p : topo.parallel_patches)
        {
            if (p.is_coupling)
                out.push_back(make_face_view(p));
        }

        return out;
    }

    std::vector<FacePatchView> noncoupling_interfaces(const TOPO::Topology &topo)
    {
        std::vector<FacePatchView> out;
        out.reserve(topo.inner_patches.size() + topo.parallel_patches.size());

        for (const auto &p : topo.inner_patches)
        {
            if (!p.is_coupling)
                out.push_back(make_face_view(p));
        }

        for (const auto &p : topo.parallel_patches)
        {
            if (!p.is_coupling)
                out.push_back(make_face_view(p));
        }

        return out;
    }

    std::vector<FacePatchView> faces_on_block(const TOPO::Topology &topo, int this_block)
    {
        std::vector<FacePatchView> out;
        out.reserve(topo.inner_patches.size() + topo.parallel_patches.size() + topo.physical_patches.size());

        for (const auto &p : topo.inner_patches)
        {
            if (p.this_block == this_block)
                out.push_back(make_face_view(p));
        }

        for (const auto &p : topo.parallel_patches)
        {
            if (p.this_block == this_block)
                out.push_back(make_face_view(p));
        }

        for (const auto &p : topo.physical_patches)
        {
            if (p.this_block == this_block)
                out.push_back(make_face_view(p));
        }

        return out;
    }

    const std::vector<TOPO::EdgePatch> &edge_patches(const TOPO::Topology &topo, TOPO::PatchKind kind)
    {
        switch (kind)
        {
        case TOPO::PatchKind::Inner:
            return topo.inner_edge_patches;
        case TOPO::PatchKind::Parallel:
            return topo.parallel_edge_patches;
        case TOPO::PatchKind::Physical:
            return topo.physical_edge_patches;
        }

        ERROR::Abort("TOPO_VIEW::edge_patches: invalid PatchKind");
        return topo.physical_edge_patches;
    }

    const std::vector<TOPO::VertexPatch> &vertex_patches(const TOPO::Topology &topo, TOPO::PatchKind kind)
    {
        switch (kind)
        {
        case TOPO::PatchKind::Inner:
            return topo.inner_vertex_patches;
        case TOPO::PatchKind::Parallel:
            return topo.parallel_vertex_patches;
        case TOPO::PatchKind::Physical:
            return topo.physical_vertex_patches;
        }

        ERROR::Abort("TOPO_VIEW::vertex_patches: invalid PatchKind");
        return topo.physical_vertex_patches;
    }

    std::vector<const TOPO::EdgePatch *> edge_patches_on_block(const TOPO::Topology &topo, TOPO::PatchKind kind, int this_block)
    {
        const auto &patches = edge_patches(topo, kind);
        std::vector<const TOPO::EdgePatch *> out;
        out.reserve(patches.size());

        for (const auto &p : patches)
        {
            if (p.this_block == this_block)
                out.push_back(&p);
        }

        return out;
    }

    std::vector<const TOPO::VertexPatch *> vertex_patches_on_block(const TOPO::Topology &topo, TOPO::PatchKind kind, int this_block)
    {
        const auto &patches = vertex_patches(topo, kind);
        std::vector<const TOPO::VertexPatch *> out;
        out.reserve(patches.size());

        for (const auto &p : patches)
        {
            if (p.this_block == this_block)
                out.push_back(&p);
        }

        return out;
    }

    bool is_interface(const FacePatchView &p)
    {
        return p.has_neighbor;
    }

    bool is_physical(const FacePatchView &p)
    {
        return p.kind == TOPO::PatchKind::Physical;
    }

    bool is_inner(const FacePatchView &p)
    {
        return p.kind == TOPO::PatchKind::Inner;
    }

    bool is_parallel(const FacePatchView &p)
    {
        return p.kind == TOPO::PatchKind::Parallel;
    }
}
