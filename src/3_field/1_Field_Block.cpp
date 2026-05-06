#include "3_field/1_Field_Block.h"
#include "0_basic/LayoutTraits.h"
#include "1_grid/BlockTraits.h"

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
    const Box3 box = LAYOUT::allocated_box_from_cells(
        GRID_TRAITS::cell_counts(blk),
        desc.location,
        desc.nghost);

    lo = box.lo;
    hi = box.hi;
}
