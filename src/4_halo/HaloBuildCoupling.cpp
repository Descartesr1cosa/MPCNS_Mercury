#include "4_halo/Halo.h"
#include "2_topology/TopologyView.h"
#include "0_basic/MPI_WRAPPER.h"
#include "4_halo/detail/HaloBuildTools.h"
#include "4_halo/detail/HaloBuildBoxMakers.h"

//==============================================================================
// Coupling (Parallel 2D Corner): build edge-corner send/recv patterns for a
// directed physical pair (src -> dst).
//
// - recv-pattern is stored on dst ranks (this=dst), and only contains recv_box.
// - send-pattern is built on src ranks from exchanged EdgeMeta, and only
//   contains send_box.
// - Tag policy: identical to normal parallel-edge halo: per-neighbor-rank local
//   sequential tags starting from 0. This is safe because coupling exchanges
//   are required to be serial (no outstanding requests across different calls).
//==============================================================================
void Halo::build_coupling_parallel_2DCorner_pattern(const std::string &src,
                                                    const std::string &dst,
                                                    StaggerLocation loc,
                                                    int nghost)
{
    int myid;
    PARALLEL::mpi_rank(&myid);

    CouplingPatternKey ckey{src, dst, loc, nghost};
    if (coupling_parallel_edge_patterns_recv.count(ckey))
        return;

    HaloPattern pat_recv;
    pat_recv.location = loc;
    pat_recv.nghost = nghost;

    // Each neighbor_rank has its own meta list and tag counter.
    std::map<int, std::vector<EdgeMeta>> meta_to_send;
    std::map<int, int> next_tag;

    // ----------------------------------------------------------------------
    // 1) Build recv side (dst) regions and EdgeMeta to send to src ranks.
    // ----------------------------------------------------------------------
    const auto &parallel_edges = TOPO_VIEW::edge_patches(*topo_, TOPO::PatchKind::Parallel);
    for (const auto &ep : parallel_edges)
    {
        if (!ep.is_coupling)
            continue;

        // Only handle direction src -> dst (nb=src, this=dst)
        if (ep.this_block_name != dst || ep.nb_block_name != src)
            continue;

        // --- map directions to neighbor (send) coordinate system ---
        int this_dir1 = ep.dir1;
        int this_dir2 = ep.dir2;

        int tar_dir1_i = HALO_TOOLS::map_dir_to_neighbor(this_dir1, ep.trans, true);
        int tar_dir2_i = HALO_TOOLS::map_dir_to_neighbor(this_dir2, ep.trans, false);

        // Ensure tar_dir1 is the mapped normal of ep.dir1
        if (ep.trans.perm[std::abs(ep.dir1) - 1] + 1 != std::abs(tar_dir1_i))
        {
            if (ep.trans.perm[std::abs(ep.dir1) - 1] + 1 == std::abs(tar_dir2_i))
                std::swap(tar_dir1_i, tar_dir2_i);
            else
            {
                std::cout << "Error for coupling parallel Edge Processing, mapping broken!\n";
                std::exit(-1);
            }
        }

        Direction this_d1 = HALO_TOOLS::int_to_direction(ep.dir1);
        Direction this_d2 = HALO_TOOLS::int_to_direction(ep.dir2);

        // --- recv region on dst side ---
        HaloRegion r;
        r.this_block = ep.this_block;   // recv (dst)
        r.neighbor_block = ep.nb_block; // send (src)
        r.this_rank = ep.this_rank;
        r.neighbor_rank = ep.nb_rank;
        r.trans = ep.trans; // dst(this) -> src(nb)

        r.recv_box = HALO_BOX::make_2DCorner_ghost_box(loc, ep.this_box_node, this_d1, this_d2, nghost);
        r.send_box = Box3{};

        const int nb_rank = ep.nb_rank;
        const int tag = next_tag[nb_rank]++; // 0,1,2,... per neighbor rank
        r.send_flag = tag;
        r.recv_flag = tag;

        pat_recv.regions.push_back(r);

        // --- meta for src side to build send_box ---
        EdgeMeta meta;
        meta.key = PatternKey{loc, nghost};
        meta.recv_rank = ep.this_rank;
        meta.send_rank = ep.nb_rank;
        meta.recv_block = ep.this_block;
        meta.send_block = ep.nb_block;
        meta.edge_node_on_send = ep.nb_box_node;
        meta.dir1_send = tar_dir1_i;
        meta.dir2_send = tar_dir2_i;
        meta.trans_recv_to_send = ep.trans;
        meta.tag = tag;

        meta_to_send[ep.nb_rank].push_back(meta);
    }

    coupling_parallel_edge_patterns_recv[ckey] = std::move(pat_recv);

    // ----------------------------------------------------------------------
    // 2) Exchange EdgeMeta and build send side (src) regions.
    // ----------------------------------------------------------------------
    std::vector<EdgeMeta> recv_metas;
    mpi_exchange_edge_meta(meta_to_send, recv_metas);

    HaloPattern pat_send;
    pat_send.location = loc;
    pat_send.nghost = nghost;

    for (const EdgeMeta &m : recv_metas)
    {
        if (m.send_rank != myid)
        {
            std::cout << "Fatal Error, received coupling EdgeMeta not for my rank: "
                      << myid << "\t but data is for " << m.send_rank << std::endl;
            std::exit(-1);
        }

        HaloRegion r;
        r.this_block = m.send_block; // src side block
        r.neighbor_block = m.recv_block;
        r.this_rank = m.send_rank;
        r.neighbor_rank = m.recv_rank;

        // send side needs transform: src(this) -> dst(nb)
        r.trans = HALO_TOOLS::inverse_transform(m.trans_recv_to_send);

        Direction send_d1 = HALO_TOOLS::int_to_direction(m.dir1_send);
        Direction send_d2 = HALO_TOOLS::int_to_direction(m.dir2_send);

        r.send_box = HALO_BOX::make_2DCorner_innerghost_box(loc,
                                                            m.edge_node_on_send,
                                                            send_d1, send_d2,
                                                            nghost);
        r.recv_box = Box3{};

        r.send_flag = m.tag;
        r.recv_flag = m.tag;

        pat_send.regions.push_back(r);
    }

    coupling_parallel_edge_patterns_send[ckey] = std::move(pat_send);
}

//==============================================================================
// Coupling (Parallel 3D Corner): build vertex-corner send/recv patterns for a
// directed physical pair (src -> dst).
//
// - recv-pattern lives on dst ranks and describes where received data lands
//   (ghost vertex box on dst blocks).
// - send-pattern lives on src ranks and describes what to pack (inner+ghost
//   vertex box on src blocks).
//
// Tag matching is ensured via VertexMeta exchange; no global tag base is used.
//==============================================================================
void Halo::build_coupling_parallel_3DCorner_pattern(const std::string &src,
                                                    const std::string &dst,
                                                    StaggerLocation loc,
                                                    int nghost)
{
    int myid;
    PARALLEL::mpi_rank(&myid);

    CouplingPatternKey ckey{src, dst, loc, nghost};
    if (coupling_parallel_vertex_patterns_recv.count(ckey))
        return;

    HaloPattern pat_recv;
    pat_recv.location = loc;
    pat_recv.nghost = nghost;

    std::map<int, std::vector<VertexMeta>> meta_to_send;
    std::map<int, int> next_tag;

    const auto &parallel_vertices = TOPO_VIEW::vertex_patches(*topo_, TOPO::PatchKind::Parallel);
    for (const auto &vp : parallel_vertices)
    {
        if (!vp.is_coupling)
            continue;

        // only build recv-side regions for this directed pair: src -> dst
        if (vp.this_block_name != dst || vp.nb_block_name != src)
            continue;

        //=========================
        // 1) map directions to sender (src) coordinate system
        //=========================
        int this_dir1 = vp.dir1;
        int this_dir2 = vp.dir2;
        int this_dir3 = vp.dir3;

        int tar_dir1_i = HALO_TOOLS::map_dir_to_neighbor(this_dir1, vp.trans, true);
        int tar_dir2_i = HALO_TOOLS::map_dir_to_neighbor(this_dir2, vp.trans, false);
        int tar_dir3_i = HALO_TOOLS::map_dir_to_neighbor(this_dir3, vp.trans, false);

        // enforce: tar_dir1 matches mapping of vp.dir1 (contact normal)
        if (vp.trans.perm[std::abs(vp.dir1) - 1] + 1 != std::abs(tar_dir1_i))
        {
            if (vp.trans.perm[std::abs(vp.dir1) - 1] + 1 == std::abs(tar_dir2_i))
                std::swap(tar_dir1_i, tar_dir2_i);
            else if (vp.trans.perm[std::abs(vp.dir1) - 1] + 1 == std::abs(tar_dir3_i))
                std::swap(tar_dir1_i, tar_dir3_i);
            else
            {
                std::cout << "Error for coupling parallel Vertex Processing, mapping broken!\n";
                std::exit(-1);
            }
        }

        Direction this_d1 = HALO_TOOLS::int_to_direction(vp.dir1);
        Direction this_d2 = HALO_TOOLS::int_to_direction(vp.dir2);
        Direction this_d3 = HALO_TOOLS::int_to_direction(vp.dir3);

        //=========================
        // 2) build recv-side HaloRegion (dst side)
        //=========================
        HaloRegion r;
        r.this_block = vp.this_block;   // dst local block
        r.neighbor_block = vp.nb_block; // src remote block
        r.this_rank = vp.this_rank;
        r.neighbor_rank = vp.nb_rank;
        r.trans = vp.trans; // dst(this) -> src(nb)

        r.recv_box = HALO_BOX::make_3DCorner_ghost_box(loc,
                                                       vp.this_box_node,
                                                       this_d1, this_d2, this_d3,
                                                       nghost);
        r.send_box = Box3{};

        int nb_rank = vp.nb_rank;
        int tag = next_tag[nb_rank]++; // unique per neighbor rank within this build
        r.send_flag = tag;
        r.recv_flag = tag;

        pat_recv.regions.push_back(r);

        //=========================
        // 3) send VertexMeta to src rank so it can build its send pattern
        //=========================
        VertexMeta meta;
        meta.key = PatternKey{loc, nghost};
        meta.recv_rank = vp.this_rank; // dst
        meta.send_rank = vp.nb_rank;   // src
        meta.recv_block = vp.this_block;
        meta.send_block = vp.nb_block;
        meta.vertex_node_on_send = vp.nb_box_node;
        meta.dir1_send = tar_dir1_i;
        meta.dir2_send = tar_dir2_i;
        meta.dir3_send = tar_dir3_i;
        meta.trans_recv_to_send = vp.trans; // dst->src
        meta.tag = tag;

        meta_to_send[vp.nb_rank].push_back(meta);
    }

    coupling_parallel_vertex_patterns_recv[ckey] = std::move(pat_recv);

    //=========================
    // 4) MPI exchange metas
    //=========================
    std::vector<VertexMeta> recv_metas;
    mpi_exchange_vertex_meta(meta_to_send, recv_metas);

    //=========================
    // 5) build send-side pattern on src ranks
    //=========================
    HaloPattern pat_send;
    pat_send.location = loc;
    pat_send.nghost = nghost;

    for (const VertexMeta &m : recv_metas)
    {
        if (m.send_rank != myid)
        {
            std::cout << "Fatal Error, received coupling VertexMeta not for my rank: "
                      << myid << "\t but data is for " << m.send_rank << std::endl;
            std::exit(-1);
        }

        HaloRegion r;
        r.this_block = m.send_block; // src side block
        r.neighbor_block = m.recv_block;
        r.this_rank = m.send_rank;
        r.neighbor_rank = m.recv_rank;

        // send side needs transform: src(this) -> dst(nb)
        r.trans = HALO_TOOLS::inverse_transform(m.trans_recv_to_send);

        Direction send_d1 = HALO_TOOLS::int_to_direction(m.dir1_send);
        Direction send_d2 = HALO_TOOLS::int_to_direction(m.dir2_send);
        Direction send_d3 = HALO_TOOLS::int_to_direction(m.dir3_send);

        r.send_box = HALO_BOX::make_3DCorner_innerghost_box(loc,
                                                            m.vertex_node_on_send,
                                                            send_d1, send_d2, send_d3,
                                                            nghost);
        r.recv_box = Box3{};

        r.send_flag = m.tag;
        r.recv_flag = m.tag;

        pat_send.regions.push_back(r);
    }

    coupling_parallel_vertex_patterns_send[ckey] = std::move(pat_send);
}
