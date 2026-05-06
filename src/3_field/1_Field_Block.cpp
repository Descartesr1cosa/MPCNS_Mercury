#include "3_field/1_Field_Block.h"

void FieldBlock::allocate(const Block &blk, const FieldDescriptor &desc_in)
{
    block = &blk;
    desc = desc_in;

    // 1. 先算逻辑 index 范围 [lo, hi)
    compute_extent(blk);

    // 2. 把逻辑范围转换成底层 Vector 的 dim1,dim2,dim3
    const int nx = hi.i - lo.i; // 总点数 = 内部 + ghost + stagger额外
    const int ny = hi.j - lo.j;
    const int nz = hi.k - lo.k;

    const int ngh = desc.nghost;

    // 对 FieldVector 来说：
    //   SetSize(dim1,dim2,dim3, ghost, ncomp);
    // 有效索引范围是：
    //   n1 ∈ [-ghost, dim1-ghost-1] 等价于 [lo, hi)
    data.SetSize(nx, ny, nz, ngh, desc.ncomp);

    allocated_ = true;
}

// 不分配数据，但绑定 block/desc/extent，供上层统一查询 lo/hi
void FieldBlock::bind_inactive(const Block &blk, const FieldDescriptor &desc_in)
{
    block = &blk;
    desc = desc_in;
    compute_extent(blk);

    data = Vector{}; // 或 data.Clear()/SetSize(0,0,0,0)
    allocated_ = false;
}

// 根据 desc.location 计算逻辑 index 范围 [lo,hi)
// 约定：
//   Block 有 Ni,Nj,Nk 个 cell；
//   Cell-centered:    i = [0 .. Ni-1] 内区间
//   Node:             i = [0 .. Ni]   (Ni+1)
//   FaceXi:           i-face = [0 .. Ni]       (Ni+1), j,k 和 cell 一样
//   FaceEt:           j-face = [0 .. Nj]       (Nj+1)
//   FaceZe:           k-face = [0 .. Nk]       (Nk+1)
//   EdgeXi:   沿 xi 的棱： i=[0..Ni-1], j,k 是 node：Nj+1, Nk+1
//   EdgeEt:   沿 eta 的棱：j=[0..Nj-1], i,k 是 node：Ni+1, Nk+1
//   EdgeZe:   沿 zeta 的棱：k=[0..Nk-1], i,j 是 node：Ni+1, Nj+1
//
// 然后再在三方向都各加 desc.nghost 层 ghost。
void FieldBlock::compute_extent(const Block &blk)
{
    auto block_node_size = [](const Block &blk) -> Int3
    {
        return {blk.mx, blk.my, blk.mz};
    };
    Int3 nc = block_node_size(blk); // {Ni, Nj, Nk}
    const int Ni = nc.i;
    const int Nj = nc.j;
    const int Nk = nc.k;

    const int g = desc.nghost;

    switch (desc.location)
    {
    case StaggerLocation::Cell:
        // cell-centered：内部 [0..Ni-1]，加 ghost 后：
        // 逻辑 index i ∈ [-g .. Ni+g-1]，写成 [lo, hi) = [-g, Ni+g)
        lo = Int3{-g, -g, -g};
        hi = Int3{Ni + g, Nj + g, Nk + g};
        break;

    case StaggerLocation::Node:
        // node: 内部 [0..Ni]，共 Ni+1 个点
        // 逻辑 index i ∈ [-g .. Ni+g]，即 [lo, hi) = [-g, Ni+1+g)
        lo = Int3{-g, -g, -g};
        hi = Int3{Ni + 1 + g, Nj + 1 + g, Nk + 1 + g};
        break;

    case StaggerLocation::FaceXi:
        // Xi-face：i 有 Ni+1 个，j,k 和 cell 一样
        // i: [-g .. Ni+g], j: [-g .. Nj+g-1], k: [-g .. Nk+g-1]
        lo = Int3{-g, -g, -g};
        hi = Int3{Ni + 1 + g, Nj + g, Nk + g};
        break;

    case StaggerLocation::FaceEt:
        // Eta-face：j 有 Nj+1 个
        // i: [-g .. Ni+g-1], j: [-g .. Nj+g], k: [-g .. Nk+g-1]
        lo = Int3{-g, -g, -g};
        hi = Int3{Ni + g, Nj + 1 + g, Nk + g};
        break;

    case StaggerLocation::FaceZe:
        // Zeta-face：k 有 Nk+1 个
        // i: [-g .. Ni+g-1], j: [-g .. Nj+g-1], k: [-g .. Nk+g]
        lo = Int3{-g, -g, -g};
        hi = Int3{Ni + g, Nj + g, Nk + 1 + g};
        break;

    case StaggerLocation::EdgeXi:
        // 沿 xi 的棱：
        //   i cell-based: [0..Ni-1]
        //   j,k node-based: [0..Nj], [0..Nk]
        // 加 ghost 后：
        lo = Int3{-g, -g, -g};
        hi = Int3{Ni + g, Nj + 1 + g, Nk + 1 + g};
        break;

    case StaggerLocation::EdgeEt:
        // 沿 eta 的棱：
        //   j cell-based: [0..Nj-1]
        //   i,k node-based: [0..Ni], [0..Nk]
        lo = Int3{-g, -g, -g};
        hi = Int3{Ni + 1 + g, Nj + g, Nk + 1 + g};
        break;

    case StaggerLocation::EdgeZe:
        // 沿 zeta 的棱：
        //   k cell-based: [0..Nk-1]
        //   i,j node-based: [0..Ni], [0..Nj]
        lo = Int3{-g, -g, -g};
        hi = Int3{Ni + 1 + g, Nj + 1 + g, Nk + g};
        break;

    default:
        // 防御式编程：如果以后加枚举忘了处理
        lo = Int3{0, 0, 0};
        hi = Int3{0, 0, 0};
        break;
    }
}