#include "4_halo/1_MPCNS_Halo.h"
#include "0_basic/MPI_WRAPPER.h"
#include "4_halo/detail/halo_build_boxmakers.h"
#include "4_halo/detail/halo_build_tools.h"

void Halo::build_parallel_3DCorner_pattern(StaggerLocation loc, int nghost)
{
    int myid;
    PARALLEL::mpi_rank(&myid);

    PatternKey key{loc, nghost};
    if (parallel_vertex_patterns_recv.count(key))
        return;

    HaloPattern pat_recv;
    pat_recv.location = loc;
    pat_recv.nghost = nghost;

    // 每个 neighbor_rank 一个 meta 列表
    std::map<int, std::vector<VertexMeta>> meta_to_send;

    // 每个 neighbor_rank 自增 tag
    std::map<int, int> next_tag;

    for (TOPO::VertexPatch &vp : topo_->parallel_vertex_patches)
    {
        if (vp.is_coupling)
            continue;
        //=========================
        // 1. 先搞清发送侧的方向
        //=========================
        // 本块上的 edge 方向（两个出界轴）
        int this_dir1 = vp.dir1;
        int this_dir2 = vp.dir2;
        int this_dir3 = vp.dir3;

        // 映射到邻居块坐标的方向
        int tar_dir1_i = HALO_TOOLS::map_dir_to_neighbor(this_dir1, vp.trans, true);
        int tar_dir2_i = HALO_TOOLS::map_dir_to_neighbor(this_dir2, vp.trans, false);
        int tar_dir3_i = HALO_TOOLS::map_dir_to_neighbor(this_dir3, vp.trans, false);

        // 保证 tar_dir1 是「接触面的法向」，和 ep.dir1 变换一致
        if (vp.trans.perm[std::abs(vp.dir1) - 1] + 1 != std::abs(tar_dir1_i))
        {
            if (vp.trans.perm[std::abs(vp.dir1) - 1] + 1 == std::abs(tar_dir2_i))
                std::swap(tar_dir1_i, tar_dir2_i);
            else if (vp.trans.perm[std::abs(vp.dir1) - 1] + 1 == std::abs(tar_dir3_i))
                std::swap(tar_dir1_i, tar_dir3_i);
            else
            {
                std::cout << "Error for parallel Vertex Processing, mapping broken!\n";
                std::exit(-1);
            }
        }

        Direction this_d1 = HALO_TOOLS::int_to_direction(vp.dir1);
        Direction this_d2 = HALO_TOOLS::int_to_direction(vp.dir2);
        Direction this_d3 = HALO_TOOLS::int_to_direction(vp.dir3);

        Direction send_d1 = HALO_TOOLS::int_to_direction(tar_dir1_i);
        Direction send_d2 = HALO_TOOLS::int_to_direction(tar_dir2_i);
        Direction send_d3 = HALO_TOOLS::int_to_direction(tar_dir3_i);

        //=========================
        // 2. 构造 recv 侧 HaloRegion
        //=========================
        HaloRegion r;
        r.this_block = vp.this_block;   // recv
        r.neighbor_block = vp.nb_block; // send
        r.this_rank = vp.this_rank;
        r.neighbor_rank = vp.nb_rank;
        r.trans = vp.trans; // recv -> send

        // 接收块上的 ghost edge 区域
        r.recv_box = HALO_BOX::make_3DCorner_ghost_box(loc, vp.this_box_node, this_d1, this_d2, this_d3, nghost);

        // recv 侧其实不会用 send_box，这里可以留空
        r.send_box = Box3{};

        // 给这一条 edge 分配一个本地 tag
        int nb_rank = vp.nb_rank;
        int tag = next_tag[nb_rank]++; // 默认 0 起，不重复就行

        r.send_flag = tag;
        r.recv_flag = tag;

        pat_recv.regions.push_back(r);

        //=========================
        // 3. 打一个 EdgeMeta，准备发给 send_rank
        //=========================
        VertexMeta meta;
        meta.key = key;
        meta.recv_rank = vp.this_rank;
        meta.send_rank = vp.nb_rank;
        meta.recv_block = vp.this_block;
        meta.send_block = vp.nb_block;
        meta.vertex_node_on_send = vp.nb_box_node;
        meta.dir1_send = tar_dir1_i;
        meta.dir2_send = tar_dir2_i;
        meta.dir3_send = tar_dir3_i;
        meta.trans_recv_to_send = vp.trans; // recv->send
        meta.tag = tag;

        meta_to_send[vp.nb_rank].push_back(meta);
    }

    parallel_vertex_patterns_recv[key] = std::move(pat_recv);

    //=========================
    // 4. MPI 交换 EdgeMeta
    //=========================
    std::vector<VertexMeta> recv_metas;
    mpi_exchange_vertex_meta(meta_to_send, recv_metas);

    //=========================
    // 5. 根据收到的 meta 构建 parallel_edge_patterns_send
    //=========================
    HaloPattern pat_send;
    pat_send.location = loc;
    pat_send.nghost = nghost;

    for (const VertexMeta &m : recv_metas)
    {

        if (m.send_rank != myid)
        {
            std::cout << "Fatal Error, received data is not for my process id: " << myid << "\t but data is for " << m.send_rank << std::endl;
            exit(-1);
        }

        HaloRegion r;
        r.this_block = m.send_block; // 在发送侧，本块是 send_block
        r.neighbor_block = m.recv_block;
        r.this_rank = m.send_rank;
        r.neighbor_rank = m.recv_rank;

        // 对发送侧来说，需要 trans: this(send) -> nb(recv)
        r.trans = HALO_TOOLS::inverse_transform(m.trans_recv_to_send);

        Direction send_d1 = HALO_TOOLS::int_to_direction(m.dir1_send);
        Direction send_d2 = HALO_TOOLS::int_to_direction(m.dir2_send);
        Direction send_d3 = HALO_TOOLS::int_to_direction(m.dir3_send);

        // send_box：发送侧 inner strip + ghost strip（按约定，dir1 为 inner）
        r.send_box = HALO_BOX::make_3DCorner_innerghost_box(loc,
                                                            m.vertex_node_on_send, send_d1, send_d2, send_d3, nghost);

        // 发送侧通常不需要用 recv_box，这里可以空着
        r.recv_box = Box3{};

        // tag 必须和 recv 侧一致
        r.send_flag = m.tag;
        r.recv_flag = m.tag;

        pat_send.regions.push_back(r);
    }

    parallel_vertex_patterns_send[key] = std::move(pat_send);
}

void Halo::build_inner_3DCorner_pattern(StaggerLocation loc, int nghost)
{
    PatternKey key{loc, nghost};
    if (inner_vertex_patterns_.count(key))
        return;

    HaloPattern pat;
    pat.location = loc;
    pat.nghost = nghost;

    for (TOPO::VertexPatch &vp : topo_->inner_vertex_patches)
    {
        if (vp.is_coupling)
            continue;
        //=====================================================================
        // 获取邻居block中edge的dir1 dir2
        const int nb_b = vp.nb_block;
        // 取出 block
        const Block &blk_nb = fld_->grd->grids(nb_b);
        // 取出Edge的范围，以node为坐标
        Box3 tar_face_range = vp.nb_box_node;
        // tar的 cell 范围
        Int3 tar_blk = HALO_TOOLS::block_node_size(blk_nb);
        // 通过 this_box_node 判断接口在本块的方向 XMinus XPlus...
        int tar_dir1, tar_dir2, tar_dir3;
        HALO_TOOLS::detect_vertex_direction(tar_face_range, tar_blk, tar_dir1, tar_dir2, tar_dir3, "build_inner_3DCorner_pattern");
        // 保证dir1为neighbor block与this block接触面的法向
        if (vp.trans.perm[abs(vp.dir1) - 1] + 1 != abs(tar_dir1))
        {
            if (vp.trans.perm[abs(vp.dir1) - 1] + 1 == abs(tar_dir2))
            {
                int temp = tar_dir1;
                tar_dir1 = tar_dir2;
                tar_dir2 = temp;
            }
            else if (vp.trans.perm[abs(vp.dir1) - 1] + 1 == abs(tar_dir3))
            {
                int temp = tar_dir1;
                tar_dir1 = tar_dir3;
                tar_dir3 = temp;
            }
            else
            {
                std::cout << "Error for inner Edge Processing, the corresponding relation is broken! !\n";
                exit(-1);
            }
        }

        //=====================================================================
        // 构造HaloRegion
        HaloRegion r;
        r.this_block = vp.this_block;   // 顶点隶属块，接收 ghost
        r.neighbor_block = vp.nb_block; // 提供 inner
        r.this_rank = vp.this_rank;
        r.neighbor_rank = vp.nb_rank;

        r.trans = vp.trans; // this -> nb

        Direction this_d1 = HALO_TOOLS::int_to_direction(vp.dir1);
        Direction this_d2 = HALO_TOOLS::int_to_direction(vp.dir2);
        Direction this_d3 = HALO_TOOLS::int_to_direction(vp.dir3);

        Direction tar_d1 = HALO_TOOLS::int_to_direction(tar_dir1);
        Direction tar_d2 = HALO_TOOLS::int_to_direction(tar_dir2);
        Direction tar_d3 = HALO_TOOLS::int_to_direction(tar_dir3);

        // recv_box：this_block 顶点 ghost 角点
        r.recv_box = HALO_BOX::make_3DCorner_ghost_box(loc, vp.this_box_node, this_d1, this_d2, this_d3, nghost);

        // send_box：neighbor_block 对应 inner 角点, tar_d1为接触面法向，取inner，其他ghost
        r.send_box = HALO_BOX::make_3DCorner_innerghost_box(loc, vp.nb_box_node, tar_d1, tar_d2, tar_d3, nghost);

        pat.regions.push_back(r);
    }

    inner_vertex_patterns_[key] = std::move(pat);
}
