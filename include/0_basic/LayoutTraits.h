#pragma once

#include "0_basic/TYPES.h"
#include "0_basic/StaggerLocation.h"

namespace LAYOUT
{
    // ------------------------------------------------------------
    // 基础约定
    // ------------------------------------------------------------
    //
    // ncells = {mx, my, mz}
    //
    // node-space half-open domain:
    //   [0, mx+1) x [0, my+1) x [0, mz+1)
    //
    // cell-space:
    //   [0, mx) x [0, my) x [0, mz)
    //

    Int3 node_size_from_cells(const Int3 &ncells);

    // delta = node_size - dof_size
    //
    // Cell   : {1,1,1}
    // Node   : {0,0,0}
    // FaceXi : {0,1,1}
    // FaceEt : {1,0,1}
    // FaceZe : {1,1,0}
    // EdgeXi : {1,0,0}
    // EdgeEt : {0,1,0}
    // EdgeZe : {0,0,1}
    Int3 dof_delta(StaggerLocation loc);

    // owned DOF box，不含 ghost
    Box3 owned_box_from_cells(const Int3 &ncells,
                              StaggerLocation loc);

    // allocated DOF box，含 ghost
    Box3 allocated_box_from_cells(const Int3 &ncells,
                                  StaggerLocation loc,
                                  int nghost);

    // node-space patch -> field DOF-space patch
    Box3 node_box_to_dof_box(StaggerLocation loc,
                             const Box3 &node_box);

    // ------------------------------------------------------------
    // Direction code:
    //   -1 XMinus, +1 XPlus
    //   -2 YMinus, +2 YPlus
    //   -3 ZMinus, +3 ZPlus
    // ------------------------------------------------------------
    bool is_valid_dir(int dir_code);
    int axis_from_dir(int dir_code); // 0/1/2
    int sign_from_dir(int dir_code); // -1/+1

    // ------------------------------------------------------------
    // Halo box rules
    // ------------------------------------------------------------

    Box3 neighbor_inner_strip_from_node_box(StaggerLocation loc,
                                            const Box3 &node_box,
                                            int dir_code,
                                            int nghost);

    Box3 ghost_strip_from_node_box(StaggerLocation loc,
                                   const Box3 &node_box,
                                   int dir_code,
                                   int nghost);

    Box3 corner2_ghost_from_node_box(StaggerLocation loc,
                                     const Box3 &node_box,
                                     int dir1,
                                     int dir2,
                                     int nghost);

    Box3 corner2_innerghost_from_node_box(StaggerLocation loc,
                                          const Box3 &node_box,
                                          int inner_dir,
                                          int ghost_dir,
                                          int nghost);

    Box3 corner3_ghost_from_node_box(StaggerLocation loc,
                                     const Box3 &node_box,
                                     int dir1,
                                     int dir2,
                                     int dir3,
                                     int nghost);

    Box3 corner3_innerghost_from_node_box(StaggerLocation loc,
                                          const Box3 &node_box,
                                          int inner_dir,
                                          int ghost_dir1,
                                          int ghost_dir2,
                                          int nghost);

    // ------------------------------------------------------------
    // Boundary / Coupling rules
    // ------------------------------------------------------------

    Box3 boundary_inner_slab_one_layer_from_cells(const Int3 &ncells,
                                                  StaggerLocation loc,
                                                  const Box3 &face_node_box,
                                                  int dir_code);

    Box3 ghost_slab_from_inner(const Box3 &inner_slab,
                               int dir_code,
                               int nghost);

    Box3 coupling_face_ghost_slab_from_cells(const Int3 &ncells,
                                             StaggerLocation loc,
                                             const Box3 &face_node_box,
                                             int dir_code,
                                             int nghost);

    Box3 coupling_edge_ghost_slab_from_cells(const Int3 &ncells,
                                             StaggerLocation loc,
                                             const Box3 &edge_node_box,
                                             int dir1,
                                             int dir2,
                                             int nghost);

    Box3 coupling_vertex_ghost_slab_from_cells(const Int3 &ncells,
                                               StaggerLocation loc,
                                               int dir1,
                                               int dir2,
                                               int dir3,
                                               int nghost);

    // ------------------------------------------------------------
    // Debug / validation
    // ------------------------------------------------------------

    bool box_empty(const Box3 &b);
    Int3 box_size(const Box3 &b);

    void assert_valid_box(const Box3 &b,
                          const char *where);

    void assert_box_inside(const Box3 &inner,
                           const Box3 &outer,
                           const char *where);
}
