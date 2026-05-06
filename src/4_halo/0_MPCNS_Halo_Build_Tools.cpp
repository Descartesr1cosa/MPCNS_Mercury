#include "4_halo/detail/halo_build_tools.h"
#include "0_basic/Error.h"

namespace HALO_TOOLS
{
    //=========================================================================
    Int3 block_node_size(const Block &blk)
    {
        return {blk.mx, blk.my, blk.mz};
    }

    Direction detect_face_direction(const Box3 &face, const Int3 &blk_mxyz, const char *where)
    {
        // 下面的判断逻辑基于假设：
        // 接口在某一方向上厚度为 1 层节点，而且刚好贴在 0 或 Ni_node 上。

        if (face.lo.i == 0 && face.hi.i == 1)
            return Direction::XMinus;
        if (face.lo.i == blk_mxyz.i && face.hi.i == blk_mxyz.i + 1)
            return Direction::XPlus;

        if (face.lo.j == 0 && face.hi.j == 1)
            return Direction::YMinus;
        if (face.lo.j == blk_mxyz.j && face.hi.j == blk_mxyz.j + 1)
            return Direction::YPlus;

        if (face.lo.k == 0 && face.hi.k == 1)
            return Direction::ZMinus;
        if (face.lo.k == blk_mxyz.k && face.hi.k == blk_mxyz.k + 1)
            return Direction::ZPlus;

        std::string oss;
        oss = "detect_face_direction: cannot determine direction from node_box";
        if (where)
            oss = oss + " in " + where;
        ERROR::Abort(oss);
    }
    //=========================================================================

    //=========================================================================
    // For 2D Corner

    // 方向编码：±1/±2/±3 -> Direction
    Direction int_to_direction(int direction)
    {
        switch (direction)
        {
        case +1:
            return Direction::XPlus;
        case -1:
            return Direction::XMinus;
        case +2:
            return Direction::YPlus;
        case -2:
            return Direction::YMinus;
        case +3:
            return Direction::ZPlus;
        case -3:
            return Direction::ZMinus;
        default:
            ERROR::Abort("int_to_direction: invalid direction");
        }
    }

    // IndexTransform 求逆：用于 sender 侧构造（你 2D/3D parallel send 都会用）
    TOPO::IndexTransform inverse_transform(const TOPO::IndexTransform &tr)
    {
        TOPO::IndexTransform inv;
        // nb[perm[a]] = sign[a]*loc[a] + offset[a]
        // => loc[a] = sign[a]*(nb[perm[a]] - offset[a])
        // 可以推到 inv.perm / inv.sign / inv.offset
        int inv_perm[3];
        for (int a = 0; a < 3; ++a)
            inv_perm[tr.perm[a]] = a;

        for (int b = 0; b < 3; ++b)
        {
            int a = inv_perm[b]; // nb 的第 b 轴，对应原来的 local 第 a 轴

            inv.perm[b] = a;          // nb[b] -> local[a]
            inv.sign[b] = tr.sign[a]; // 系数同 sign[a]
            // offset[a] 是原来公式 nb[...] = sign[a]*local[a] + offset[a] 里的 offset[a]
            // 逆公式里变成 -sign[a]*offset[a]
            int off_a = (a == 0   ? tr.offset.i
                         : a == 1 ? tr.offset.j
                                  : tr.offset.k);
            int off_b = -tr.sign[a] * off_a;

            if (b == 0)
                inv.offset.i = off_b;
            if (b == 1)
                inv.offset.j = off_b;
            if (b == 2)
                inv.offset.k = off_b;
        }
        return inv;
    }

    // 把“本地方向编码”映射到“邻居方向编码”
    // dir_local: ±1/±2/±3（表示 local 的 X/Y/Z 及符号）
    // is_norm: 是否为法向（法向要多一次符号翻转，你现在已有此规则）
    int map_dir_to_neighbor(int dir_local, const TOPO::IndexTransform &tr, bool is_norm)
    {
        int axis_local = std::abs(dir_local) - 1;
        int axis_nb = tr.perm[axis_local];
        int s_local = (dir_local > 0) ? tr.sign[axis_local] : -tr.sign[axis_local]; // sign为翻转，dir_local正负为大小号面
        int s_multi = (is_norm) ? -s_local : s_local;                               // 法向面需要额外翻转
        return s_multi * (axis_nb + 1);
    }

    static inline void push_dir(int code, int &dir1, int &dir2, int &count)
    {
        if (count == 0)
            dir1 = code;
        else if (count == 1)
            dir2 = code;
        ++count;
    }

    // 2D: 从 edge_node 推断两个方向编码（±1/±2/±3）
    // dir1/dir2: 输出，顺序不保证；
    // blk_mxyz： block_node_size = {mx,my,mz}
    // where：用于异常信息定位
    void detect_edge_direction(const Box3 &edge_node,
                               const Int3 &blk_mxyz,
                               int &dir1,
                               int &dir2,
                               const char *where)
    {
        dir1 = 0;
        dir2 = 0;
        int count = 0;

        // X-
        if (edge_node.lo.i == 0 && edge_node.hi.i == 1)
            push_dir(-1, dir1, dir2, count);
        // X+
        if (edge_node.lo.i == blk_mxyz.i && edge_node.hi.i == blk_mxyz.i + 1)
            push_dir(+1, dir1, dir2, count);

        // Y-
        if (edge_node.lo.j == 0 && edge_node.hi.j == 1)
            push_dir(-2, dir1, dir2, count);
        // Y+
        if (edge_node.lo.j == blk_mxyz.j && edge_node.hi.j == blk_mxyz.j + 1)
            push_dir(+2, dir1, dir2, count);

        // Z-
        if (edge_node.lo.k == 0 && edge_node.hi.k == 1)
            push_dir(-3, dir1, dir2, count);
        // Z+
        if (edge_node.lo.k == blk_mxyz.k && edge_node.hi.k == blk_mxyz.k + 1)
            push_dir(+3, dir1, dir2, count);

        // Edge 必须恰好贴两条边（两个方向）
        if (count != 2)
        {
            std::string oss;
            oss += ("detect_edge_dirs: expected 2 boundary directions, got " + std::to_string(count));
            if (where)
                oss = oss + " in " + where;

            ERROR::Abort(oss);
        }

        // 防止同一轴同时匹配 +/-（理论上不该发生，但加一道保险）
        if (std::abs(dir1) == std::abs(dir2))
        {
            std::string oss;
            oss += ("detect_edge_dirs: invalid dirs (same axis) dir1=" + std::to_string(dir1) + " dir2=" + std::to_string(dir2));
            if (where)
                oss = oss + " in " + where;
            ERROR::Abort(oss);
        }
    }
    //=========================================================================

    //=========================================================================
    // For 3D Corner

    static inline void push_dir3(int code, int &d1, int &d2, int &d3, int &cnt)
    {
        if (cnt == 0)
            d1 = code;
        else if (cnt == 1)
            d2 = code;
        else if (cnt == 2)
            d3 = code;
        ++cnt;
    }

    void detect_vertex_direction(const Box3 &vertex_node,
                                 const Int3 &blk_mxyz,
                                 int &dir1,
                                 int &dir2,
                                 int &dir3,
                                 const char *where)
    {
        dir1 = 0;
        dir2 = 0;
        dir3 = 0;
        int cnt = 0;

        // X-
        if (vertex_node.lo.i == 0 && vertex_node.hi.i == 1)
            push_dir3(-1, dir1, dir2, dir3, cnt);
        // X+
        if (vertex_node.lo.i == blk_mxyz.i && vertex_node.hi.i == blk_mxyz.i + 1)
            push_dir3(+1, dir1, dir2, dir3, cnt);

        // Y-
        if (vertex_node.lo.j == 0 && vertex_node.hi.j == 1)
            push_dir3(-2, dir1, dir2, dir3, cnt);
        // Y+
        if (vertex_node.lo.j == blk_mxyz.j && vertex_node.hi.j == blk_mxyz.j + 1)
            push_dir3(+2, dir1, dir2, dir3, cnt);

        // Z-
        if (vertex_node.lo.k == 0 && vertex_node.hi.k == 1)
            push_dir3(-3, dir1, dir2, dir3, cnt);
        // Z+
        if (vertex_node.lo.k == blk_mxyz.k && vertex_node.hi.k == blk_mxyz.k + 1)
            push_dir3(+3, dir1, dir2, dir3, cnt);

        // Vertex 必须恰好贴三条边（三个方向）
        if (cnt != 3)
        {
            std::string oss;
            oss += ("detect_vertex_direction: expected 3 boundary directions, got " + std::to_string(cnt));
            if (where)
                oss = oss + " in " + where;

            ERROR::Abort(oss);
        }

        // 防御：不允许重复轴
        int a1 = std::abs(dir1), a2 = std::abs(dir2), a3 = std::abs(dir3);
        if (a1 == a2 || a1 == a3 || a2 == a3)
        {
            std::string oss;
            oss += ("detect_vertex_direction: invalid dirs (repeated axis)  dir1=" + std::to_string(dir1) + " dir2=" + std::to_string(dir2) + " dir3=" + std::to_string(dir3));
            if (where)
                oss = oss + " in " + where;
            ERROR::Abort(oss);
        }
    }
    //=========================================================================

    //=========================================================================
    // T: (i,j,k) --> (io,jo,ko)
    void apply_transform(const TOPO::IndexTransform &T,
                         int i, int j, int k,
                         int &io, int &jo, int &ko)
    {
        const int loc[3] = {i, j, k};
        int tar[3] = {0, 0, 0};
        const int off[3] = {T.offset.i, T.offset.j, T.offset.k};

        for (int d = 0; d < 3; ++d)
            tar[T.perm[d]] = T.sign[d] * loc[d] + off[d];

        io = tar[0];
        jo = tar[1];
        ko = tar[2];
    }

    void pack_to_neighbor_order(FieldBlock &fb,
                                const Box3 &sb,
                                int ncomp,
                                const TOPO::IndexTransform &T, // this -> nb
                                std::vector<double> &out)
    {
        const int32_t n_total =
            (sb.hi.i - sb.lo.i) * (sb.hi.j - sb.lo.j) * (sb.hi.k - sb.lo.k) * ncomp;
        out.resize(n_total);

        // int loc_lo[3] = {sb.lo.i, sb.lo.j, sb.lo.k};
        // int loc_hi[3] = {sb.hi.i - 1, sb.hi.j - 1, sb.hi.k - 1}; // closed interval
        int len_loc[3] = {sb.hi.i - sb.lo.i, sb.hi.j - sb.lo.j, sb.hi.k - sb.lo.k};

        int offset[3] = {T.offset.i, T.offset.j, T.offset.k};
        // int tar1[3], tar2[3], tar_ref[3];

        // for (int d = 0; d < 3; ++d)
        //     tar1[T.perm[d]] = T.sign[d] * loc_lo[d] + offset[d];
        // for (int d = 0; d < 3; ++d)
        //     tar2[T.perm[d]] = T.sign[d] * loc_hi[d] + offset[d];

        // tar_ref[0] = (tar1[0] <= tar2[0]) ? tar1[0] : tar2[0];
        // tar_ref[1] = (tar1[1] <= tar2[1]) ? tar1[1] : tar2[1];
        // tar_ref[2] = (tar1[2] <= tar2[2]) ? tar1[2] : tar2[2];
        int tar_ref[3] = {0, 0, 0}; // nb box 的“最小角”锚点：用 sb.lo 映射得到
        {
            int loc0[3] = {sb.lo.i, sb.lo.j, sb.lo.k};
            for (int a = 0; a < 3; ++a)
                tar_ref[T.perm[a]] = T.sign[a] * loc0[a] + offset[a];
        }

        int len_nb[3] = {0, 0, 0};
        for (int d = 0; d < 3; ++d)
            len_nb[T.perm[d]] = len_loc[d];

        int ijk[3], tar_ijk[3];
        for (ijk[0] = sb.lo.i; ijk[0] < sb.hi.i; ++ijk[0])
            for (ijk[1] = sb.lo.j; ijk[1] < sb.hi.j; ++ijk[1])
                for (ijk[2] = sb.lo.k; ijk[2] < sb.hi.k; ++ijk[2])
                {
                    for (int ii = 0; ii < 3; ++ii)
                        tar_ijk[T.perm[ii]] = T.sign[ii] * ijk[ii] + offset[ii];

                    int ri = tar_ijk[0] - tar_ref[0];
                    int rj = tar_ijk[1] - tar_ref[1];
                    int rk = tar_ijk[2] - tar_ref[2];

                    int32_t base = ((ri * len_nb[1] + rj) * len_nb[2] + rk) * ncomp;

                    for (int m = 0; m < ncomp; ++m)
                        out[base + m] = fb(ijk[0], ijk[1], ijk[2], m);
                }
    }

    void unpack_to_coupling_buffer(CouplingBufferBlock &bufblk,
                                   const Box3 &rb,
                                   int ncomp,
                                   const std::vector<double> &in)
    {
        const int ni = rb.hi.i - rb.lo.i;
        const int nj = rb.hi.j - rb.lo.j;
        const int nk = rb.hi.k - rb.lo.k;

        const int32_t n_total = ni * nj * nk * ncomp;
        if ((int32_t)in.size() < n_total)
        {
            std::cout << "Fatal Error!!! unpack_to_coupling_buffer size mismatch\n";
            std::exit(-1);
        }

        for (int ii = 0; ii < ni; ++ii)
            for (int jj = 0; jj < nj; ++jj)
                for (int kk = 0; kk < nk; ++kk)
                {
                    int i = rb.lo.i + ii;
                    int j = rb.lo.j + jj;
                    int k = rb.lo.k + kk;

                    int32_t base = (((ii * nj + jj) * nk) + kk) * ncomp;
                    for (int m = 0; m < ncomp; ++m)
                        bufblk(i, j, k, m) = in[base + m];
                }
    }

    bool box_equal(const Box3 &a, const Box3 &b)
    {
        return (a.lo.i == b.lo.i && a.lo.j == b.lo.j && a.lo.k == b.lo.k &&
                a.hi.i == b.hi.i && a.hi.j == b.hi.j && a.hi.k == b.hi.k);
    }
}