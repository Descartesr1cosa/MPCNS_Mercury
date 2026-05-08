#include "4_halo/detail/HaloBuildBoxMakers.h"
#include "0_basic/LayoutTraits.h"

namespace
{
    int dir_code(Direction dir)
    {
        switch (dir)
        {
        case Direction::XMinus:
            return -1;
        case Direction::XPlus:
            return +1;
        case Direction::YMinus:
            return -2;
        case Direction::YPlus:
            return +2;
        case Direction::ZMinus:
            return -3;
        case Direction::ZPlus:
            return +3;
        }

        return 0;
    }
}

namespace HALO_BOX
{
    Box3 make_1DCorner_inner_box(StaggerLocation loc, const Box3 &face_node, Direction dir, int nghost)
    {
        return LAYOUT::neighbor_inner_strip_from_node_box(loc, face_node, dir_code(dir), nghost);
    }

    Box3 make_1DCorner_ghost_box(StaggerLocation loc, const Box3 &face_node, Direction dir, int nghost)
    {
        return LAYOUT::ghost_strip_from_node_box(loc, face_node, dir_code(dir), nghost);
    }

    Box3 make_2DCorner_ghost_box(StaggerLocation loc, const Box3 &edge_node,
                                 Direction dir1, Direction dir2, int nghost)
    {
        return LAYOUT::corner2_ghost_from_node_box(
            loc, edge_node, dir_code(dir1), dir_code(dir2), nghost);
    }

    // 发送则一定是一个为inner一个为ghost，约定，dir1为inner
    Box3 make_2DCorner_innerghost_box(StaggerLocation loc, const Box3 &edge_node,
                                      Direction dir1, Direction dir2, int nghost)
    {
        return LAYOUT::corner2_innerghost_from_node_box(
            loc, edge_node, dir_code(dir1), dir_code(dir2), nghost);
    }

    Box3 make_3DCorner_ghost_box(StaggerLocation loc, const Box3 &vertex_node,
                                 Direction dir1, Direction dir2, Direction dir3, int nghost)
    {
        return LAYOUT::corner3_ghost_from_node_box(
            loc, vertex_node, dir_code(dir1), dir_code(dir2), dir_code(dir3), nghost);
    }

    // 发送则一定是一个为inner两个为ghost，约定，dir1为inner
    Box3 make_3DCorner_innerghost_box(StaggerLocation loc, const Box3 &vertex_node,
                                      Direction dir1, Direction dir2, Direction dir3, int nghost)
    {
        return LAYOUT::corner3_innerghost_from_node_box(
            loc, vertex_node, dir_code(dir1), dir_code(dir2), dir_code(dir3), nghost);
    }
} // namespace HALO_BOX
