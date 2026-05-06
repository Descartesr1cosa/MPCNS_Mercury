#include "4_halo/detail/halo_build_boxmakers.h"
#include <stdexcept>
#include "0_basic/Error.h"

namespace HALO_BOX
{
    namespace detail
    {
        static inline int delta_on_axis(StaggerLocation loc, int ax) // ax:0=i,1=j,2=k
        {
            // delta=1: cell-like; delta=0: node-like
            switch (loc)
            {
            case StaggerLocation::Cell:
                return 1;
            case StaggerLocation::Node:
                return 0;
            case StaggerLocation::FaceXi:
                return (ax == 0 ? 0 : 1);
            case StaggerLocation::FaceEt:
                return (ax == 1 ? 0 : 1);
            case StaggerLocation::FaceZe:
                return (ax == 2 ? 0 : 1);
            case StaggerLocation::EdgeXi:
                return (ax == 0 ? 1 : 0);
            case StaggerLocation::EdgeEt:
                return (ax == 1 ? 1 : 0);
            case StaggerLocation::EdgeZe:
                return (ax == 2 ? 1 : 0);
            default:
                return 1;
            }
        }

        static inline int axis_from_dir(Direction dir) // 0/1/2
        {
            switch (dir)
            {
            case Direction::XMinus:
            case Direction::XPlus:
                return 0;
            case Direction::YMinus:
            case Direction::YPlus:
                return 1;
            case Direction::ZMinus:
            case Direction::ZPlus:
                return 2;
            }
            return 0;
        }

        static void switch_range_inner(Box3 &box, StaggerLocation loc, Direction dir, int nghost)
        {
            const int ax = axis_from_dir(dir);
            const int delta = delta_on_axis(loc, ax);
            const int exclude = (delta == 0) ? 1 : 0; // node-like 法向：排除接口层；cell-like：不排除

            auto apply_minus = [&](int &lo, int &hi)
            {
                lo = lo + exclude;
                hi = lo + nghost;
            };
            auto apply_plus = [&](int &lo, int &hi)
            {
                hi = hi - exclude;
                lo = hi - nghost;
            };

            switch (dir)
            {
            case Direction::XMinus:
                apply_minus(box.lo.i, box.hi.i);
                break;
            case Direction::XPlus:
                apply_plus(box.lo.i, box.hi.i);
                break;
            case Direction::YMinus:
                apply_minus(box.lo.j, box.hi.j);
                break;
            case Direction::YPlus:
                apply_plus(box.lo.j, box.hi.j);
                break;
            case Direction::ZMinus:
                apply_minus(box.lo.k, box.hi.k);
                break;
            case Direction::ZPlus:
                apply_plus(box.lo.k, box.hi.k);
                break;
            }
        }
        // static void switch_range_inner(Box3 &box, Direction dir, int nghost)
        // {
        //     switch (dir)
        //     {
        //     case Direction::XMinus:
        //         // 左 inner：i ∈ [0 .. g-1] = [0,g)
        //         box.lo.i = 0;
        //         box.hi.i = nghost;
        //         break;

        //     case Direction::XPlus:
        //         // 右 inner strip：i ∈ [Ni-g .. Ni) = [hi_int-g, hi_int)
        //         box.lo.i = box.hi.i - nghost;
        //         box.hi.i = box.hi.i;
        //         break;

        //     case Direction::YMinus:
        //         box.lo.j = 0;
        //         box.hi.j = nghost;
        //         break;

        //     case Direction::YPlus:
        //         box.lo.j = box.hi.j - nghost;
        //         box.hi.j = box.hi.j;
        //         break;

        //     case Direction::ZMinus:
        //         box.lo.k = 0;
        //         box.hi.k = nghost;
        //         break;

        //     case Direction::ZPlus:
        //         box.lo.k = box.hi.k - nghost;
        //         box.hi.k = box.hi.k;
        //         break;
        //     }
        // }

        static void switch_range_ghost(Box3 &box, Direction dir, int nghost)
        {
            // 这里输入的Box3 &face为特定StaggerLocation的坐标
            // Cell [0,mx-1] = [0,mx)
            switch (dir)
            {
            case Direction::XMinus:
                // 左 ghost：i ∈ [-g .. -1] = [-g,0)
                box.lo.i = -nghost;
                box.hi.i = 0;
                break;

            case Direction::XPlus:
                // 右 ghost：i ∈ [Ni .. Ni+g) = [hi_int, hi_int+g)
                box.lo.i = box.hi.i;
                box.hi.i = box.hi.i + nghost;
                break;

            case Direction::YMinus:
                box.lo.j = -nghost;
                box.hi.j = 0;
                break;

            case Direction::YPlus:
                box.lo.j = box.hi.j;
                box.hi.j = box.hi.j + nghost;
                break;

            case Direction::ZMinus:
                box.lo.k = -nghost;
                box.hi.k = 0;
                break;

            case Direction::ZPlus:
                box.lo.k = box.hi.k;
                box.hi.k = box.hi.k + nghost;
                break;
            }
        }

        // 把 face_node 映射成该 loc 对应的“内部 DOF box”
        // 这一步只做 hi 的 -1/0 规则
        static Box3 make_base_dof_box_from_node(StaggerLocation loc, const Box3 &face_node)
        {
            Box3 box;
            box.lo = face_node.lo;

            // face_node 是 node 空间 [lo,hi)，此处构造该 loc 的“内部 DOF”范围 [lo,hi)
            // 规则完全等价于你现有实现：
            //   Cell  : hi = (i-1,j-1,k-1)
            //   FaceXi: hi = (i  ,j-1,k-1)
            //   FaceEt: hi = (i-1,j  ,k-1)
            //   FaceZe: hi = (i-1,j-1,k  )
            //   EdgeXi: hi = (i-1,j  ,k  )
            //   EdgeEt: hi = (i  ,j-1,k  )
            //   EdgeZe: hi = (i  ,j  ,k-1)
            switch (loc)
            {
            case StaggerLocation::Cell:
                box.hi = {face_node.hi.i - 1, face_node.hi.j - 1, face_node.hi.k - 1};
                break;
            case StaggerLocation::FaceXi:
                box.hi = {face_node.hi.i, face_node.hi.j - 1, face_node.hi.k - 1};
                break;
            case StaggerLocation::FaceEt:
                box.hi = {face_node.hi.i - 1, face_node.hi.j, face_node.hi.k - 1};
                break;
            case StaggerLocation::FaceZe:
                box.hi = {face_node.hi.i - 1, face_node.hi.j - 1, face_node.hi.k};
                break;
            case StaggerLocation::EdgeXi:
                box.hi = {face_node.hi.i - 1, face_node.hi.j, face_node.hi.k};
                break;
            case StaggerLocation::EdgeEt:
                box.hi = {face_node.hi.i, face_node.hi.j - 1, face_node.hi.k};
                break;
            case StaggerLocation::EdgeZe:
                box.hi = {face_node.hi.i, face_node.hi.j, face_node.hi.k - 1};
                break;
            default:
                ERROR::Abort("make_base_dof_box_1D: unsupported StaggerLocation");
            }

            return box;
        }

        // 用于安全检查：避免两个方向在同一轴上
        static int dir_axis(Direction d)
        {
            switch (d)
            {
            case Direction::XMinus:
            case Direction::XPlus:
                return 1;
            case Direction::YMinus:
            case Direction::YPlus:
                return 2;
            case Direction::ZMinus:
            case Direction::ZPlus:
                return 3;
            }
            return 0;
        }

    } // namespace detail

    Box3 make_1DCorner_inner_box(StaggerLocation loc, const Box3 &face_node, Direction dir, int nghost)
    {
        Box3 box = detail::make_base_dof_box_from_node(loc, face_node);
        detail::switch_range_inner(box, loc, dir, nghost);
        return box;
    }

    Box3 make_1DCorner_ghost_box(StaggerLocation loc, const Box3 &face_node, Direction dir, int nghost)
    {
        Box3 box = detail::make_base_dof_box_from_node(loc, face_node);
        detail::switch_range_ghost(box, dir, nghost);
        return box;
    }

    Box3 make_2DCorner_ghost_box(StaggerLocation loc, const Box3 &edge_node,
                                 Direction dir1, Direction dir2, int nghost)
    {

        if (detail::dir_axis(dir1) == detail::dir_axis(dir2))
            ERROR::Abort("make_2DCorner_ghost_box: dir1/dir2 are on the same axis");

        // edge_node 是“沿棱的一条 node strip”：[lo, hi)
        Box3 box = detail::make_base_dof_box_from_node(loc, edge_node);

        // 在两个出界方向上依次改 ghost 范围
        detail::switch_range_ghost(box, dir1, nghost);
        detail::switch_range_ghost(box, dir2, nghost);
        return box;
    }

    // 发送则一定是一个为inner一个为ghost，约定，dir1为inner
    Box3 make_2DCorner_innerghost_box(StaggerLocation loc, const Box3 &edge_node,
                                      Direction dir1, Direction dir2, int nghost)
    {
        if (detail::dir_axis(dir1) == detail::dir_axis(dir2))
            ERROR::Abort("make_2DCorner_innerghost_box: dir1/dir2 are on the same axis");

        Box3 box = detail::make_base_dof_box_from_node(loc, edge_node);

        // 先在第一个方向取 inner strip，再在第二个方向取 ghost strip
        detail::switch_range_inner(box, loc, dir1, nghost);
        detail::switch_range_ghost(box, dir2, nghost);
        return box;
    }

    //=========================================================================
    // For 3D Corner

    // vertex_node 是顶点附近的 node 区域（通常 [i,i+1)×[j,j+1)×[k,k+1)）
    Box3 make_3DCorner_ghost_box(StaggerLocation loc, const Box3 &vertex_node,
                                 Direction dir1, Direction dir2, Direction dir3, int nghost)
    {
        // vertex_node 是一个 [i,i+1)×[j,j+1)×[k,k+1) 的小 node 立方
        // 先映射到内部 那个顶点，再沿三个方向扩展 ghost
        Box3 box = detail::make_base_dof_box_from_node(loc, vertex_node);
        // 三个方向依次做 ghost 扩展，就得到三维角区:
        // 例如 XMinus+YMinus+ZMinus -> i∈[-g,0), j∈[-g,0), k∈[-g,0)
        detail::switch_range_ghost(box, dir1, nghost);
        detail::switch_range_ghost(box, dir2, nghost);
        detail::switch_range_ghost(box, dir3, nghost);
        return box;
    }

    // 发送则一定是一个为inner两个为ghost，约定，dir1为inner
    Box3 make_3DCorner_innerghost_box(StaggerLocation loc, const Box3 &vertex_node,
                                      Direction dir1, Direction dir2, Direction dir3, int nghost)
    {
        Box3 box = detail::make_base_dof_box_from_node(loc, vertex_node);

        // 三个方向依次裁剪为 inner strip：
        // 例如 XMinus+YMinus+ZMinus -> i∈[0,g), j∈[0,g), k∈[0,g)
        detail::switch_range_inner(box, loc, dir1, nghost); // 约定 dir1 是 inner
        detail::switch_range_ghost(box, dir2, nghost);
        detail::switch_range_ghost(box, dir3, nghost);
        return box;
    }
    //=========================================================================

}