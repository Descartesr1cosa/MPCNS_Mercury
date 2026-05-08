#include "4_halo/Halo.h"
#include "2_topology/TopologyView.h"
#include "0_basic/MPI_WRAPPER.h"
#include "4_halo/detail/HaloBuildTools.h"
#include "4_halo/detail/HaloBuildBoxMakers.h"

void Halo::build_parallel_2DCorner_pattern(StaggerLocation loc, int nghost)
{
    int myid;
    PARALLEL::mpi_rank(&myid);

    PatternKey key{loc, nghost};
    if (parallel_edge_patterns_recv.count(key))
        return;

    HaloPattern pat_recv;
    pat_recv.location = loc;
    pat_recv.nghost = nghost;

    // 每个 neighbor_rank 一个 meta 列表
    std::map<int, std::vector<EdgeMeta>> meta_to_send;

    // 每个 neighbor_rank 自增 tag
    std::map<int, int> next_tag;

    const auto &parallel_edges = TOPO_VIEW::edge_patches(*topo_, TOPO::PatchKind::Parallel);
    for (const auto &ep : parallel_edges)
    {
        if (ep.is_coupling)
            continue;
        //=========================
        // 1. 先搞清发送侧的方向
        //=========================
        // 本块上的 edge 方向（两个出界轴）
        int this_dir1 = ep.dir1;
        int this_dir2 = ep.dir2;

        // 映射到邻居块坐标的方向
        int tar_dir1_i = HALO_TOOLS::map_dir_to_neighbor(this_dir1, ep.trans, true);
        int tar_dir2_i = HALO_TOOLS::map_dir_to_neighbor(this_dir2, ep.trans, false);

        // 保证 tar_dir1 是「接触面的法向」，和 ep.dir1 变换一致
        if (ep.trans.perm[std::abs(ep.dir1) - 1] + 1 != std::abs(tar_dir1_i))
        {
            if (ep.trans.perm[std::abs(ep.dir1) - 1] + 1 == std::abs(tar_dir2_i))
                std::swap(tar_dir1_i, tar_dir2_i);
            else
            {
                std::cout << "Error for parallel Edge Processing, mapping broken!\n";
                std::exit(-1);
            }
        }

        Direction this_d1 = HALO_TOOLS::int_to_direction(ep.dir1);
        Direction this_d2 = HALO_TOOLS::int_to_direction(ep.dir2);
        Direction send_d1 = HALO_TOOLS::int_to_direction(tar_dir1_i);
        Direction send_d2 = HALO_TOOLS::int_to_direction(tar_dir2_i);

        //=========================
        // 2. 构造 recv 侧 HaloRegion
        //=========================
        HaloRegion r;
        r.this_block = ep.this_block;   // recv
        r.neighbor_block = ep.nb_block; // send
        r.this_rank = ep.this_rank;
        r.neighbor_rank = ep.nb_rank;
        r.trans = ep.trans; // recv -> send

        // 接收块上的 ghost edge 区域
        r.recv_box = HALO_BOX::make_2DCorner_ghost_box(loc, ep.this_box_node, this_d1, this_d2, nghost);

        // recv 侧其实不会用 send_box，这里可以留空
        r.send_box = Box3{};

        // 给这一条 edge 分配一个本地 tag
        int nb_rank = ep.nb_rank;
        int tag = next_tag[nb_rank]++; // 默认 0 起，不重复就行

        r.send_flag = tag;
        r.recv_flag = tag;

        pat_recv.regions.push_back(r);

        //=========================
        // 3. 打一个 EdgeMeta，准备发给 send_rank
        //=========================
        EdgeMeta meta;
        meta.key = key;
        meta.recv_rank = ep.this_rank;
        meta.send_rank = ep.nb_rank;
        meta.recv_block = ep.this_block;
        meta.send_block = ep.nb_block;
        meta.edge_node_on_send = ep.nb_box_node;
        meta.dir1_send = tar_dir1_i;
        meta.dir2_send = tar_dir2_i;
        meta.trans_recv_to_send = ep.trans; // recv->send
        meta.tag = tag;

        meta_to_send[ep.nb_rank].push_back(meta);
    }

    parallel_edge_patterns_recv[key] = std::move(pat_recv);

    //=========================
    // 4. MPI 交换 EdgeMeta
    //=========================
    std::vector<EdgeMeta> recv_metas;
    mpi_exchange_edge_meta(meta_to_send, recv_metas);

    //=========================
    // 5. 根据收到的 meta 构建 parallel_edge_patterns_send
    //=========================
    HaloPattern pat_send;
    pat_send.location = loc;
    pat_send.nghost = nghost;

    for (const EdgeMeta &m : recv_metas)
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

        // send_box：发送侧 inner strip + ghost strip（按约定，dir1 为 inner）
        r.send_box = HALO_BOX::make_2DCorner_innerghost_box(loc,
                                                            m.edge_node_on_send, send_d1, send_d2, nghost);

        // 发送侧通常不需要用 recv_box，这里可以空着
        r.recv_box = Box3{};

        // tag 必须和 recv 侧一致
        r.send_flag = m.tag;
        r.recv_flag = m.tag;

        pat_send.regions.push_back(r);
    }

    parallel_edge_patterns_send[key] = std::move(pat_send);
}

void Halo::build_inner_2DCorner_pattern(StaggerLocation loc, int nghost)
{
    PatternKey key{loc, nghost};
    if (inner_edge_patterns_.find(key) != inner_edge_patterns_.end())
        return; // 已经建过

    HaloPattern pat;
    pat.location = loc;
    pat.nghost = nghost;
    pat.regions.clear();

    const auto &inner_edges = TOPO_VIEW::edge_patches(*topo_, TOPO::PatchKind::Inner);
    for (const auto &ep : inner_edges)
    {
        if (ep.is_coupling)
            continue;
        //=====================================================================
        // 获取邻居block中edge的dir1 dir2
        const int nb_b = ep.nb_block;
        // 取出 block
        const Block &blk_nb = fld_->grd->grids(nb_b);
        // 取出Edge的范围，以node为坐标
        Box3 tar_face_range = ep.nb_box_node;
        // tar的 cell 范围
        Int3 tar_blk = HALO_TOOLS::block_node_size(blk_nb);
        // 通过 this_box_node 判断接口在本块的方向 XMinus XPlus...
        int tar_dir1, tar_dir2;
        HALO_TOOLS::detect_edge_direction(tar_face_range, tar_blk, tar_dir1, tar_dir2, "build_inner_edge_pattern");
        // 保证dir1为neighbor block与this block接触面的法向
        if (ep.trans.perm[abs(ep.dir1) - 1] + 1 != abs(tar_dir1))
        {
            if (ep.trans.perm[abs(ep.dir1) - 1] + 1 == abs(tar_dir2))
            {
                // 安全检测
                int temp = tar_dir1;
                tar_dir1 = tar_dir2;
                tar_dir2 = temp;
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

        // edge 隶属块 = 接收方（this_block）
        r.this_block = ep.this_block;   // 这里填 ghost
        r.neighbor_block = ep.nb_block; // 这里提供 inner
        r.this_rank = ep.this_rank;
        r.neighbor_rank = ep.nb_rank;

        // trans: 按 HaloTypes 的语义，是 this -> nb 的映射
        r.trans = ep.trans;

        Direction this_d1 = HALO_TOOLS::int_to_direction(ep.dir1);
        Direction this_d2 = HALO_TOOLS::int_to_direction(ep.dir2);

        Direction tar_d1 = HALO_TOOLS::int_to_direction(tar_dir1);
        Direction tar_d2 = HALO_TOOLS::int_to_direction(tar_dir2);

        // recv_box 在 this_block（edge 所属块）上
        // 注意，recv的两个方向都是虚网格故可直接调用
        r.recv_box = HALO_BOX::make_2DCorner_ghost_box(loc, ep.this_box_node, this_d1, this_d2, nghost);

        // send_box 在 neighbor_block 上，从邻居的 inner strip 取数据
        // 注意tar_d1为两块交界面的法向，处理角区时应该取inner
        r.send_box = HALO_BOX::make_2DCorner_innerghost_box(loc, ep.nb_box_node, tar_d1, tar_d2, nghost);

        pat.regions.push_back(r);
    }

    inner_edge_patterns_[key] = std::move(pat);
}
