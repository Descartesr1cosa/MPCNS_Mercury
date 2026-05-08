// include/4_halo/detail/halo_build_boxmakers.h
#pragma once

#include "4_halo/Halo_Type.h"

namespace HALO_BOX
{

    Box3 make_1DCorner_inner_box(StaggerLocation loc, const Box3 &face_node, Direction dir, int nghost);
    Box3 make_1DCorner_ghost_box(StaggerLocation loc, const Box3 &face_node, Direction dir, int nghost);

    // edge_node 是该 edge 在 node 空间里的交集区域（一般是一条线段 [lo,hi)）
    // 接受可以直接两个方向均为ghost，代表角区
    Box3 make_2DCorner_ghost_box(StaggerLocation loc, const Box3 &edge_node,
                                 Direction dir1, Direction dir2, int nghost);
    // 发送则一定是一个为inner一个为ghost，约定，dir1为inner
    Box3 make_2DCorner_innerghost_box(StaggerLocation loc, const Box3 &edge_node,
                                      Direction dir1, Direction dir2, int nghost);

    // 在 Cell 空间里，为三维角区（vertex：三个方向 ghost）构造 box
    // vertex_node 是顶点附近的 node 区域（通常 [i,i+1)×[j,j+1)×[k,k+1)）
    Box3 make_3DCorner_ghost_box(StaggerLocation loc, const Box3 &vertex_node,
                                 Direction dir1, Direction dir2, Direction dir3, int nghost);
    Box3 make_3DCorner_innerghost_box(StaggerLocation loc, const Box3 &vertex_node,
                                      Direction dir1, Direction dir2, Direction dir3, int nghost);
} // namespace HALO_BOX