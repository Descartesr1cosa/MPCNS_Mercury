// include/4_halo/detail/halo_build_tools.h
#pragma once

#include <cstddef> // for nullptr_t

#include "4_halo/Halo_Type.h"
#include "3_field/1_Field_Block.h"
#include "3_field/Coupling_Type.h"

namespace HALO_TOOLS
{
    //=========================================================================
    // 返回 {mx,my,mz}
    Int3 block_node_size(const Block &blk);

    // 从 node face box 推断方向（1DCorner）
    // where: 用于把报错信息标上调用位置，便于调试（可传 nullptr）
    Direction detect_face_direction(const Box3 &face, const Int3 &blk_mxyz, const char *where = nullptr);
    //=========================================================================

    //=========================================================================
    // For 2D Corner

    // 方向编码：±1/±2/±3 -> Direction
    Direction int_to_direction(int direction);

    // IndexTransform 求逆：用于 sender 侧构造（你 2D/3D parallel send 都会用）
    TOPO::IndexTransform inverse_transform(const TOPO::IndexTransform &tr);

    // 把“本地方向编码”映射到“邻居方向编码”
    // dir_local: ±1/±2/±3（表示 local 的 X/Y/Z 及符号）
    // is_norm: 是否为法向（法向要多一次符号翻转，你现在已有此规则）
    int map_dir_to_neighbor(int dir_local, const TOPO::IndexTransform &tr, bool is_norm);

    // 2D: 从 edge_node 推断两个方向编码（±1/±2/±3）
    // dir1/dir2: 输出，顺序不保证；
    // blk_mxyz： block_node_size = {mx,my,mz}
    // where：用于异常信息定位
    void detect_edge_direction(const Box3 &edge_node,
                               const Int3 &blk_mxyz,
                               int &dir1,
                               int &dir2,
                               const char *where);
    //=========================================================================

    //=========================================================================
    // For 3D Corner

    // 3D: 从 vertex_node box 推断三个方向编码（±1/±2/±3）
    // 假设：该 vertex_node 在三个轴上厚度为 1，并且分别贴在 0 或 Ni_node 上
    void detect_vertex_direction(const Box3 &vertex_node,
                                 const Int3 &blk_mxyz,
                                 int &dir1,
                                 int &dir2,
                                 int &dir3,
                                 const char *where);

    //=========================================================================

    //=========================================================================
    // T: (i,j,k) --> (io,jo,ko)
    void apply_transform(const TOPO::IndexTransform &T,
                         int i, int j, int k,
                         int &io, int &jo, int &ko);

    int face_2form_orientation_sign(const TOPO::IndexTransform &T,
                                    int source_axis);

    // For Coupling 1D 2D 3D Corner
    // pack：按“邻居坐标顺序”打包 fb --> out
    void pack_to_neighbor_order(FieldBlock &fb,
                                const Box3 &sb,
                                int ncomp,
                                const TOPO::IndexTransform &T, // this -> nb
                                std::vector<double> &out);

    // 按 recv_box 的 i/j/k 顺序解包到 CouplingBufferBlock in --> bufblk
    void unpack_to_coupling_buffer(CouplingBufferBlock &bufblk,
                                   const Box3 &rb,
                                   int ncomp,
                                   const std::vector<double> &in);

    bool box_equal(const Box3 &a, const Box3 &b);
} // namespace HALO_TOOLS
