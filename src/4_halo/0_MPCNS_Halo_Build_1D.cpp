#include "4_halo/1_MPCNS_Halo.h"
#include "4_halo/detail/halo_build_tools.h"
#include "4_halo/detail/halo_build_boxmakers.h"

void Halo::build_inner_1DCorner_pattern(StaggerLocation loc, int nghost)
{
    PatternKey key{loc, nghost};
    // 已经 build 过了，直接返回，避免重复构造
    if (inner_patterns_.find(key) != inner_patterns_.end())
        return;

    HaloPattern pat;
    pat.location = loc;
    pat.nghost = nghost;
    pat.regions.clear();

    // 遍历所有 inner_patches（同一 rank 的块-块接口）
    for (const TOPO::InterfacePatch &patch : topo_->inner_patches)
    {
        if (patch.is_coupling)
            continue; // 跳过耦合接口

        const int this_b = patch.this_block;
        const int nb_b = patch.nb_block;

        // 取出两个 block
        const Block &blk_this = fld_->grd->grids(this_b); // blocks_ 是 Field 里的
        const Block &blk_nb = fld_->grd->grids(nb_b);

        // 取出面的范围，以node为坐标
        Box3 my_face_range = patch.this_box_node;
        Box3 tar_face_range = patch.nb_box_node;

        // 各自的 cell 范围
        Int3 my_blk = HALO_TOOLS::block_node_size(blk_this);
        Int3 tar_blk = HALO_TOOLS::block_node_size(blk_nb);

        // 通过 this_box_node 判断接口在本块的方向 XMinus XPlus...
        Direction dir_my = HALO_TOOLS::detect_face_direction(my_face_range, my_blk, "build_inner_1DCorner_pattern");
        Direction dir_tar = HALO_TOOLS::detect_face_direction(tar_face_range, tar_blk, "build_inner_1DCorner_pattern");

        // --------- region 1: this_block 的 inner 到 nb_block 的 ghost  ---------
        HaloRegion r1;
        r1.this_block = this_b;   // 发送方
        r1.neighbor_block = nb_b; // 接收方
        r1.this_rank = patch.this_rank;
        r1.neighbor_rank = patch.nb_rank;
        r1.trans = patch.trans; // 从 this_block 逻辑坐标 -> neighbor_block 逻辑坐标

        // ---------处理范围 ---------
        r1.send_box = HALO_BOX::make_1DCorner_inner_box(loc, my_face_range, dir_my, nghost);
        r1.recv_box = HALO_BOX::make_1DCorner_ghost_box(loc, tar_face_range, dir_tar, nghost);

        // ---------加入regions ---------
        pat.regions.push_back(r1);
    }
    // 存储供exchange_inner使用
    inner_patterns_[key] = std::move(pat);
}

void Halo::build_parallel_1DCorner_pattern(StaggerLocation loc, int nghost)
{
    PatternKey key{loc, nghost};
    if (parallel_patterns_.find(key) != parallel_patterns_.end())
        return; // 已经建过，直接返回

    HaloPattern pat;
    pat.location = loc;
    pat.nghost = nghost;
    pat.regions.clear();

    // 遍历所有 parallel_patches（跨 rank 的接口）
    for (const TOPO::InterfacePatch &patch : topo_->parallel_patches)
    {
        if (patch.is_coupling)
            continue; // 跳过耦合接口

        const int this_b = patch.this_block;

        // 只知道本 rank 上的 block，邻居 block 不在本 rank 上
        const Block &blk_this = fld_->grd->grids(this_b);

        // 本块 interface 的 node 区域
        Box3 my_face_range = patch.this_box_node;
        Int3 my_blk = HALO_TOOLS::block_node_size(blk_this);

        // 判定在本块上的方向
        Direction dir_my = HALO_TOOLS::detect_face_direction(my_face_range, my_blk, "build_parallel_1DCorner_pattern");

        HaloRegion r;
        r.this_block = this_b;             // 本块（既是发送方也是接收方，从 MPI 角度看）
        r.neighbor_block = patch.nb_block; // 对面的 block index，Parallel 情况下一般=-1，仅做记录
        r.this_rank = patch.this_rank;
        r.neighbor_rank = patch.nb_rank; // 真实的邻居 rank
        r.trans = patch.trans;           // 从 this_block 逻辑坐标 -> neighbor_block 逻辑坐标
        r.send_flag = patch.send_flag;
        r.recv_flag = patch.recv_flag;

        // 对 parallel 来说，send_box / recv_box 都在本块的 FieldBlock 空间里：
        //   recv_box: 本块 ghost 区域（MPI 接收后，往这里 unpack）
        //   send_box: 本块 inner strip（打包后，MPI 发送给 neighbor_rank）
        r.recv_box = HALO_BOX::make_1DCorner_ghost_box(loc, my_face_range, dir_my, nghost);
        r.send_box = HALO_BOX::make_1DCorner_inner_box(loc, my_face_range, dir_my, nghost);

        pat.regions.push_back(r);
    }

    // 存起来，供 exchange_parallel 使用
    parallel_patterns_[key] = std::move(pat);
}
