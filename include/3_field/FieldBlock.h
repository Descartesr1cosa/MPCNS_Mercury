#pragma once
#include "0_basic/TYPES.h"
#include "1_grid/1_MPCNS_Grid.h"
#include "3_field/FieldDescriptor.h"
#include "3_field/FieldArray.h"

// 知道自己 belong to 哪个 Block；
// 知道自己是 cell / face / edge；
// 知道自己在逻辑 index 上的范围（含 ghost）；
// 提供 operator()(i,j,k,comp) 访问。

class FieldBlock
{
public:
    FieldBlock() {}
    ~FieldBlock() = default;

    // 真正分配数据
    void allocate(const Block &blk, const FieldDescriptor &desc_in);

    // 不分配数据，但绑定 block/desc/extent，供上层统一查询 lo/hi
    void bind_inactive(const Block &blk, const FieldDescriptor &desc_in);

    bool is_allocated() const { return allocated_; }

    // 逻辑 index 范围（包含 ghost）：[lo, hi)，半开区间
    const Int3 &get_lo() const { return lo; }
    const Int3 &get_hi() const { return hi; }

    const FieldDescriptor &descriptor() const { return desc; }
    const Block &get_block() const { return *block; }

    // 计算域的范围
    Int3 inner_lo() const { return Int3{0, 0, 0}; }
    Int3 inner_hi() const { return Int3{hi.i - desc.nghost, hi.j - desc.nghost, hi.k - desc.nghost}; }

    // 访问接口：传入“逻辑 index”，内部自动加偏移
    inline double &operator()(int i, int j, int k, int m)
    {
        // int ii = i - lo.i;
        // int jj = j - lo.j;
        // int kk = k - lo.k;
        return data(i, j, k, m);
    }

    inline double operator()(int i, int j, int k, int m) const
    {
        return const_cast<Vector &>(data)(i, j, k, m);
    }

    // 底层数据（如果某些老代码想直接操作数组也可以拿到）
    Vector &raw_data() { return data; }

private:
    const Block *block = nullptr; // 指向 Grid 中的 Block
    FieldDescriptor desc;         // 拷贝一份 descriptor

    Int3 lo; // 逻辑 index 下界
    Int3 hi; // 逻辑 index 上界（不含）

    // 你可以用自己的多维数组类型，比如类似 Phy_Tensor 那样：
    Vector data; // 假定有 SetSize(Ni, Nj, Nk, ncomp) / operator()(i,j,k,m)

    bool allocated_ = false;

    void compute_extent(const Block &blk); // 根据 desc.location 计算 lo/hi
};
