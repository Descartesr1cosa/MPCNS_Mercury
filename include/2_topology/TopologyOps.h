#pragma once

#include "2_topology/TopologyTypes.h"

namespace TOPO
{
    // ------------------------------------------------------------
    // Int3 component helpers
    // ------------------------------------------------------------

    int get_comp(const Int3 &x, int axis0);
    void set_comp(Int3 &x, int axis0, int value);

    // ------------------------------------------------------------
    // Box / node-space helpers
    // ------------------------------------------------------------

    // sub/sup 是 Grid_Boundary::Pre_process 之后的 0-based 闭区间端点，
    // 这里统一转为 node-space 半开区间 [lo, hi)。
    Box3 make_node_box_from_subsup(const int sub[3], const int sup[3]);

    // 保留旧接口的实现入口，方便逐步迁移。
    void fill_node_box_from_subsup(const int sub[3], const int sup[3], Box3 &box);

    // ------------------------------------------------------------
    // IndexTransform helpers
    // ------------------------------------------------------------

    IndexTransform identity_transform();

    // nb[ perm[a] ] = sign[a] * local[a] + offset[a]
    Int3 map_node_point(const Int3 &p_local, const IndexTransform &tr);

    // 对 node-space half-open box 做 transform。
    // 输入输出均为 [lo, hi)。
    Box3 transform_node_box(const Box3 &box, const IndexTransform &tr);

    // 根据 Grid_Boundary 中的 Transform / tar_Transform / sub / sup 构造 IndexTransform。
    // this_transform / nb_transform 是 Pre_process 后的 Transform 数组，取值 0/1/2。
    IndexTransform make_index_transform_from_boundary_arrays(
        const int this_sub[3],
        const int this_sup[3],
        const int this_transform[3],
        const int nb_sub[3],
        const int nb_sup[3],
        const int nb_transform[3],
        int this_direction,
        int nb_direction,
        int dimension,
        const char *context);

    // ------------------------------------------------------------
    // Patch priority
    // ------------------------------------------------------------

    int patch_kind_priority(PatchKind k);

    // 目前沿用你已有策略：
    //   非 coupling 优先于 coupling；
    //   同类下 Inner > Parallel > Physical。
    int patch_priority(PatchKind k, bool is_coupling);

    // ------------------------------------------------------------
    // Direction helpers for edge / vertex patches
    // ------------------------------------------------------------

    // 给定两个 face direction，返回 edge 自身沿哪个轴，返回 1/2/3。
    // 例如 dir1 = X, dir2 = Y，则 edge 沿 Z，返回 3。
    int edge_axis1_from_face_dirs(int dir1, int dir2);

    // 把两个方向整理为：
    //   dir1 == owner_dir
    //   dir2 == the other dir
    void order_edge_dirs_by_owner(int owner_dir, int &dir1, int &dir2);

    // 把 vertex 的三个方向整理为：
    //   dir1 == owner_dir
    //   dir2/dir3 保持其余两个方向
    void order_vertex_dirs_by_owner(int owner_dir,
                                    int &dir1,
                                    int &dir2,
                                    int &dir3);
}