#pragma once

#include <string>
#include <vector>

#include "2_topology/TopologyTypes.h"

namespace TOPO_VIEW
{
    // Query facade over topology-side patch adjacency.  These views expose
    // codim-1 faces and the underlying codim-2/codim-3 patch lists for halo
    // pattern construction; they do not describe ghost storage or transfer
    // buffers.  A future halo layer should introduce StorageAddress or
    // HaloAddress when expanding these regions for a particular field.
    struct FacePatchView
    {
        TOPO::PatchKind kind;
        bool has_neighbor = false;
        bool is_coupling = false;

        int this_rank = 0;
        int nb_rank = -1;

        int this_block = 0;
        int nb_block = -1;

        std::string this_block_name;
        std::string nb_block_name;

        Box3 this_box_node{};
        Box3 nb_box_node{};

        int direction = 0;
        int nb_direction = 0;

        TOPO::IndexTransform trans{};

        int32_t send_flag = 0;
        int32_t recv_flag = 0;

        int bc_id = 0;
        std::string bc_name;

        const TOPO::InterfacePatch *interface_patch = nullptr;
        const TOPO::PhysicalPatch *physical_patch = nullptr;
    };

    std::vector<FacePatchView> inner_faces(const TOPO::Topology &topo);
    std::vector<FacePatchView> parallel_faces(const TOPO::Topology &topo);
    std::vector<FacePatchView> physical_faces(const TOPO::Topology &topo);
    std::vector<FacePatchView> interface_faces(const TOPO::Topology &topo);
    std::vector<FacePatchView> all_faces(const TOPO::Topology &topo);
    std::vector<FacePatchView> coupling_interfaces(const TOPO::Topology &topo);
    std::vector<FacePatchView> noncoupling_interfaces(const TOPO::Topology &topo);
    std::vector<FacePatchView> faces_on_block(const TOPO::Topology &topo, int this_block);

    // Codim-2/codim-3 adjacency accessors for direct edge/corner halo plans.
    // They intentionally return topology patches, not ghost entities.
    const std::vector<TOPO::EdgePatch> &edge_patches(const TOPO::Topology &topo, TOPO::PatchKind kind);
    const std::vector<TOPO::VertexPatch> &vertex_patches(const TOPO::Topology &topo, TOPO::PatchKind kind);
    std::vector<const TOPO::EdgePatch *> edge_patches_on_block(const TOPO::Topology &topo, TOPO::PatchKind kind, int this_block);
    std::vector<const TOPO::VertexPatch *> vertex_patches_on_block(const TOPO::Topology &topo, TOPO::PatchKind kind, int this_block);

    bool is_interface(const FacePatchView &p);
    bool is_physical(const FacePatchView &p);
    bool is_inner(const FacePatchView &p);
    bool is_parallel(const FacePatchView &p);
}
