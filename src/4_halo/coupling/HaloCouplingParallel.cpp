#include "4_halo/Halo.h"
#include "2_topology/TopologyView.h"
#include "0_basic/MPI_WRAPPER.h"
#include "4_halo/detail/HaloBuildTools.h"
#include "4_halo/detail/HaloBuildBoxMakers.h"

namespace
{
void pack_to_neighbor_order_scaled(FieldBlock &fb,
                                   const Box3 &sb,
                                   int ncomp,
                                   const TOPO::IndexTransform &T,
                                   double factor,
                                   std::vector<double> &out)
{
    HALO_TOOLS::pack_to_neighbor_order(fb, sb, ncomp, T, out);
    if (factor == 1.0)
        return;

    for (double &v : out)
        v *= factor;
}
} // namespace

void Halo::coupling_parallel_face(const std::string &src, const std::string &dst)
{
    if (!fld_->has_coupling_pair(src, dst))
        return;

    CouplingBuffersForPair &cb = fld_->coupling_buffers(src, dst);
    const CouplingPairDesc &desc = fld_->coupling_pair(src, dst);

    if (cb.parallel_face.empty())
        return;

    // 强制串行：按 channel 逐个 stage
    for (size_t cid = 0; cid < desc.channels.size(); ++cid)
    {
        const CouplingChannelSpec &ch = desc.channels[cid];
        const int ncomp = ch.ncomp;

        // ---------------- collect send/recv items ----------------
        struct SendItem
        {
            const TOPO::InterfacePatch *p;
            int fid;
            int this_block;
            Box3 sb;
            double factor = 1.0;
            int32_t len = 0;
        };
        struct RecvItem
        {
            const TOPO::InterfacePatch *p;
            CouplingBufferBlock *buf;
            Box3 rb;
            int32_t len = 0;
        };

        std::vector<SendItem> sends;
        std::vector<RecvItem> recvs;

        const auto parallel_faces = TOPO_VIEW::parallel_faces(*topo_);
        for (size_t ip = 0; ip < parallel_faces.size(); ++ip)
        {
            const auto &p = parallel_faces[ip];
            if (!p.is_coupling)
                continue;

            // 发送侧：src(this) -> dst(nb)
            if (p.this_block_name == src && p.nb_block_name == dst)
            {
                std::string tag = ch.tag;
                double factor = 1.0;
                StaggerLocation src_loc = ch.location;
                if (coupling_channel_needs_form_transfer_(ch))
                {
                    const int dst_axis = coupling_form_axis_from_location_(ch.location);
                    const int src_axis = coupling_src_axis_from_src_to_dst_transform_(p.interface_patch->trans, dst_axis);
                    tag = find_triplet_field_name_(ch.tag, ch.value_kind, src_axis);
                    src_loc = coupling_form_location_from_axis_(ch.value_kind, src_axis);
                    factor = static_cast<double>(
                        coupling_form_orientation_sign_src_to_dst_(ch, p.interface_patch->trans, src_axis));
                }
                const int fid = fld_->field_id(tag);

                const Block &blk = fld_->grd->grids(p.this_block);
                const Int3 blk_node = HALO_TOOLS::block_node_size(blk);

                Direction dir = HALO_TOOLS::detect_face_direction(
                    p.this_box_node, blk_node, "coupling_parallel_face");

                Box3 sb = HALO_BOX::make_1DCorner_inner_box(
                    src_loc, p.this_box_node, dir, ch.nghost);

                const int Ni = sb.hi.i - sb.lo.i;
                const int Nj = sb.hi.j - sb.lo.j;
                const int Nk = sb.hi.k - sb.lo.k;
                if (Ni <= 0 || Nj <= 0 || Nk <= 0)
                    continue;

                SendItem it;
                it.p = p.interface_patch;
                it.fid = fid;
                it.this_block = p.this_block;
                it.sb = sb;
                it.factor = factor;
                sends.push_back(std::move(it));
            }

            // 接收侧：dst(this) 接收 src(nb)
            if (p.this_block_name == dst && p.nb_block_name == src)
            {
                if (cid >= cb.parallel_face.size() || ip >= cb.parallel_face[cid].size())
                {
                    std::cout << "Fatal: coupling_parallel_face buffer index out of range\n";
                    std::exit(-1);
                }

                CouplingBufferBlock &buf = cb.parallel_face[cid][ip];
                if (!buf.allocated)
                    continue; // 若你希望更严格，这里可以直接 fatal

                RecvItem it;
                it.p = p.interface_patch;
                it.buf = &buf;
                it.rb = buf.box;
                recvs.push_back(std::move(it));
            }
        }

        // ---------------- ensure reusable buffers ----------------
        if (send_buf.size() < sends.size())
            send_buf.resize(sends.size());
        if (req_send.size() < sends.size())
            req_send.resize(sends.size());
        if (stat_send.size() < sends.size())
            stat_send.resize(sends.size());
        if (length.size() < sends.size())
            length.resize(sends.size());

        if (recv_buf.size() < recvs.size())
            recv_buf.resize(recvs.size());
        if (req_recv.size() < recvs.size())
            req_recv.resize(recvs.size());
        if (stat_recv.size() < recvs.size())
            stat_recv.resize(recvs.size());
        if (length_corner_recv.size() < recvs.size())
            length_corner_recv.resize(recvs.size());

        // ---------------- pack sends ----------------
        for (size_t is = 0; is < sends.size(); ++is)
        {
            SendItem &it = sends[is];
            FieldBlock &fb = fld_->field(it.fid, it.this_block);

            // this(src) -> nb(dst)
            pack_to_neighbor_order_scaled(fb, it.sb, ncomp, it.p->trans, it.factor, send_buf[is]);

            it.len = (int32_t)send_buf[is].size();
            length[is] = it.len;
        }

        // ---------------- prepare recv buffers ----------------
        for (size_t ir = 0; ir < recvs.size(); ++ir)
        {
            RecvItem &it = recvs[ir];
            const Box3 &rb = it.rb;

            const int32_t n_total =
                (rb.hi.i - rb.lo.i) * (rb.hi.j - rb.lo.j) * (rb.hi.k - rb.lo.k) * ncomp;

            it.len = n_total;
            length_corner_recv[ir] = n_total;

            if ((int32_t)recv_buf[ir].size() < n_total)
                recv_buf[ir].resize(n_total);
        }

        // ---------------- MPI (recv first, then send) ----------------
        for (size_t ir = 0; ir < recvs.size(); ++ir)
        {
            const TOPO::InterfacePatch *p = recvs[ir].p;
            PARALLEL::mpi_data_recv(p->nb_rank, p->recv_flag,
                                    recv_buf[ir].data(), length_corner_recv[ir],
                                    &(req_recv[ir]));
        }

        for (size_t is = 0; is < sends.size(); ++is)
        {
            const TOPO::InterfacePatch *p = sends[is].p;
            PARALLEL::mpi_data_send(p->nb_rank, p->send_flag,
                                    send_buf[is].data(), length[is],
                                    &(req_send[is]));
        }

        if (!sends.empty())
        {
            int nsend = (int)sends.size();
            PARALLEL::mpi_wait(nsend, req_send.data(), stat_send.data());
        }
        if (!recvs.empty())
        {
            int nrecv = (int)recvs.size();
            PARALLEL::mpi_wait(nrecv, req_recv.data(), stat_recv.data());
        }

        // 全局严格流水线，保留 barrier
        PARALLEL::mpi_barrier();

        // ---------------- unpack to coupling buffers ----------------
        for (size_t ir = 0; ir < recvs.size(); ++ir)
        {
            RecvItem &it = recvs[ir];
            HALO_TOOLS::unpack_to_coupling_buffer(*(it.buf), it.rb, ncomp, recv_buf[ir]);
        }

        PARALLEL::mpi_barrier();
    }
}

void Halo::coupling_parallel_face(const std::string &src, const std::string &dst, std::vector<int32_t> &field_cids)
{
    if (!fld_->has_coupling_pair(src, dst))
        return;

    CouplingBuffersForPair &cb = fld_->coupling_buffers(src, dst);
    const CouplingPairDesc &desc = fld_->coupling_pair(src, dst);

    if (cb.parallel_face.empty())
        return;

    // Find corresponding channel to field_name cid
    for (auto cid_for_string : field_cids)
    {
        const CouplingChannelSpec &ch = desc.channels[cid_for_string];
        const int ncomp = ch.ncomp;

        // ---------------- collect send/recv items ----------------
        struct SendItem
        {
            const TOPO::InterfacePatch *p;
            int fid;
            int this_block;
            Box3 sb;
            double factor = 1.0;
            int32_t len = 0;
        };
        struct RecvItem
        {
            const TOPO::InterfacePatch *p;
            CouplingBufferBlock *buf;
            Box3 rb;
            int32_t len = 0;
        };

        std::vector<SendItem> sends;
        std::vector<RecvItem> recvs;

        const auto parallel_faces = TOPO_VIEW::parallel_faces(*topo_);
        for (size_t ip = 0; ip < parallel_faces.size(); ++ip)
        {
            const auto &p = parallel_faces[ip];
            if (!p.is_coupling)
                continue;

            // 发送侧：src(this) -> dst(nb)
            if (p.this_block_name == src && p.nb_block_name == dst)
            {
                std::string tag = ch.tag;
                double factor = 1.0;
                StaggerLocation src_loc = ch.location;
                if (coupling_channel_needs_form_transfer_(ch))
                {
                    const int dst_axis = coupling_form_axis_from_location_(ch.location);
                    const int src_axis = coupling_src_axis_from_src_to_dst_transform_(p.interface_patch->trans, dst_axis);
                    tag = find_triplet_field_name_(ch.tag, ch.value_kind, src_axis);
                    src_loc = coupling_form_location_from_axis_(ch.value_kind, src_axis);
                    factor = static_cast<double>(
                        coupling_form_orientation_sign_src_to_dst_(ch, p.interface_patch->trans, src_axis));
                }
                const int fid = fld_->field_id(tag);

                const Block &blk = fld_->grd->grids(p.this_block);
                const Int3 blk_node = HALO_TOOLS::block_node_size(blk);

                Direction dir = HALO_TOOLS::detect_face_direction(
                    p.this_box_node, blk_node, "coupling_parallel_face");

                Box3 sb = HALO_BOX::make_1DCorner_inner_box(
                    src_loc, p.this_box_node, dir, ch.nghost);

                const int Ni = sb.hi.i - sb.lo.i;
                const int Nj = sb.hi.j - sb.lo.j;
                const int Nk = sb.hi.k - sb.lo.k;
                if (Ni <= 0 || Nj <= 0 || Nk <= 0)
                    continue;

                SendItem it;
                it.p = p.interface_patch;
                it.fid = fid;
                it.this_block = p.this_block;
                it.sb = sb;
                it.factor = factor;
                sends.push_back(std::move(it));
            }

            // 接收侧：dst(this) 接收 src(nb)
            if (p.this_block_name == dst && p.nb_block_name == src)
            {
                if (cid_for_string >= cb.parallel_face.size() || ip >= cb.parallel_face[cid_for_string].size())
                {
                    std::cout << "Fatal: coupling_parallel_face buffer index out of range\n";
                    std::exit(-1);
                }

                CouplingBufferBlock &buf = cb.parallel_face[cid_for_string][ip];
                if (!buf.allocated)
                    continue; // 若你希望更严格，这里可以直接 fatal

                RecvItem it;
                it.p = p.interface_patch;
                it.buf = &buf;
                it.rb = buf.box;
                recvs.push_back(std::move(it));
            }
        }

        // ---------------- ensure reusable buffers ----------------
        if (send_buf.size() < sends.size())
            send_buf.resize(sends.size());
        if (req_send.size() < sends.size())
            req_send.resize(sends.size());
        if (stat_send.size() < sends.size())
            stat_send.resize(sends.size());
        if (length.size() < sends.size())
            length.resize(sends.size());

        if (recv_buf.size() < recvs.size())
            recv_buf.resize(recvs.size());
        if (req_recv.size() < recvs.size())
            req_recv.resize(recvs.size());
        if (stat_recv.size() < recvs.size())
            stat_recv.resize(recvs.size());
        if (length_corner_recv.size() < recvs.size())
            length_corner_recv.resize(recvs.size());

        // ---------------- pack sends ----------------
        for (size_t is = 0; is < sends.size(); ++is)
        {
            SendItem &it = sends[is];
            FieldBlock &fb = fld_->field(it.fid, it.this_block);

            auto in_box = [](const Box3 &b, int i, int j, int k)
            {
                return (i >= b.lo.i && i < b.hi.i &&
                        j >= b.lo.j && j < b.hi.j &&
                        k >= b.lo.k && k < b.hi.k);
            };

            static bool printed = false;
            if (!printed && ch.tag == "Eface_eta" && in_box(it.sb, 49, 49, 0))
            {
                printed = true;
                int io, jo, ko;
                HALO_TOOLS::apply_transform(it.p->trans, 49, 50, 0, io, jo, ko);

                std::cout << "[PROBE] src=" << src << " dst=" << dst
                          << " this_block=" << it.p->this_block
                          << " nb_block=" << it.p->nb_block << "\n";
                std::cout << "[PROBE] sb: i[" << it.sb.lo.i << "," << it.sb.hi.i << ") "
                          << "j[" << it.sb.lo.j << "," << it.sb.hi.j << ") "
                          << "k[" << it.sb.lo.k << "," << it.sb.hi.k << ")\n";
                std::cout << "[PROBE] this_box_node: i[" << it.p->this_box_node.lo.i << "," << it.p->this_box_node.hi.i << ") "
                          << "j[" << it.p->this_box_node.lo.j << "," << it.p->this_box_node.hi.j << ") "
                          << "k[" << it.p->this_box_node.lo.k << "," << it.p->this_box_node.hi.k << ")\n";
                std::cout << "[PROBE] nb_box_node: i[" << it.p->nb_box_node.lo.i << "," << it.p->nb_box_node.hi.i << ") "
                          << "j[" << it.p->nb_box_node.lo.j << "," << it.p->nb_box_node.hi.j << ") "
                          << "k[" << it.p->nb_box_node.lo.k << "," << it.p->nb_box_node.hi.k << ")\n";

                std::cout << "[PROBE] T perm=(" << it.p->trans.perm[0] << "," << it.p->trans.perm[1] << "," << it.p->trans.perm[2] << ") "
                          << "sign=(" << it.p->trans.sign[0] << "," << it.p->trans.sign[1] << "," << it.p->trans.sign[2] << ") "
                          << "off=(" << it.p->trans.offset.i << "," << it.p->trans.offset.j << "," << it.p->trans.offset.k << ")\n";

                std::cout << "[PROBE] map (49,50,0) -> (" << io << "," << jo << "," << ko << ")\n";
            }

            // this(src) -> nb(dst)
            pack_to_neighbor_order_scaled(fb, it.sb, ncomp, it.p->trans, it.factor, send_buf[is]);

            it.len = (int32_t)send_buf[is].size();
            length[is] = it.len;
        }

        // ---------------- prepare recv buffers ----------------
        for (size_t ir = 0; ir < recvs.size(); ++ir)
        {
            RecvItem &it = recvs[ir];
            const Box3 &rb = it.rb;

            const int32_t n_total =
                (rb.hi.i - rb.lo.i) * (rb.hi.j - rb.lo.j) * (rb.hi.k - rb.lo.k) * ncomp;

            it.len = n_total;
            length_corner_recv[ir] = n_total;

            if ((int32_t)recv_buf[ir].size() < n_total)
                recv_buf[ir].resize(n_total);
        }

        // ---------------- MPI (recv first, then send) ----------------
        for (size_t ir = 0; ir < recvs.size(); ++ir)
        {
            const TOPO::InterfacePatch *p = recvs[ir].p;
            PARALLEL::mpi_data_recv(p->nb_rank, p->recv_flag,
                                    recv_buf[ir].data(), length_corner_recv[ir],
                                    &(req_recv[ir]));
        }

        for (size_t is = 0; is < sends.size(); ++is)
        {
            const TOPO::InterfacePatch *p = sends[is].p;
            PARALLEL::mpi_data_send(p->nb_rank, p->send_flag,
                                    send_buf[is].data(), length[is],
                                    &(req_send[is]));
        }

        if (!sends.empty())
        {
            int nsend = (int)sends.size();
            PARALLEL::mpi_wait(nsend, req_send.data(), stat_send.data());
        }
        if (!recvs.empty())
        {
            int nrecv = (int)recvs.size();
            PARALLEL::mpi_wait(nrecv, req_recv.data(), stat_recv.data());
        }

        // 全局严格流水线，保留 barrier
        PARALLEL::mpi_barrier();

        // ---------------- unpack to coupling buffers ----------------
        for (size_t ir = 0; ir < recvs.size(); ++ir)
        {
            RecvItem &it = recvs[ir];
            HALO_TOOLS::unpack_to_coupling_buffer(*(it.buf), it.rb, ncomp, recv_buf[ir]);
        }

        PARALLEL::mpi_barrier();
    }
}

void Halo::coupling_parallel_edge(const std::string &src, const std::string &dst)
{
    // 2D+ 才有 edge corner
    const int dim = fld_->grd->dimension;
    if (dim < 2)
        return;

    if (!fld_->has_coupling_pair(src, dst))
        return;

    CouplingBuffersForPair &cb = fld_->coupling_buffers(src, dst);
    const CouplingPairDesc &desc = fld_->coupling_pair(src, dst);

    if (cb.parallel_edge.empty())
        return;

    // 强制串行：按 channel 逐个 stage
    for (size_t cid = 0; cid < desc.channels.size(); ++cid)
    {
        const CouplingChannelSpec &ch = desc.channels[cid];
        const int ncomp = ch.ncomp;

        // ---- 若 topo 根本没有该 pair 的 parallel coupling edge patch，就跳过（避免无谓 fatal）----
        bool has_any_patch = false;
        const auto &parallel_edges = TOPO_VIEW::edge_patches(*topo_, TOPO::PatchKind::Parallel);
        for (const auto &ep : parallel_edges)
        {
            if (!ep.is_coupling)
                continue;
            if ((ep.this_block_name == dst && ep.nb_block_name == src) ||
                (ep.this_block_name == src && ep.nb_block_name == dst))
            {
                has_any_patch = true;
                break;
            }
        }
        if (!has_any_patch)
            continue;

        // ---- pattern 取出（必须已在 build 阶段建好并缓存）----
        CouplingPatternKey pkey{src, dst, ch.location, ch.nghost};

        auto it_send = coupling_parallel_edge_patterns_send.find(pkey);
        auto it_recv = coupling_parallel_edge_patterns_recv.find(pkey);

        if (it_send == coupling_parallel_edge_patterns_send.end() ||
            it_recv == coupling_parallel_edge_patterns_recv.end())
        {
            std::cout << "Fatal: coupling_parallel_edge pattern not built for pair "
                      << src << " -> " << dst << " (loc=" << (int)ch.location
                      << ", nghost=" << ch.nghost << ")\n";
            std::exit(-1);
        }

        const HaloPattern &pat_send = it_send->second;
        const HaloPattern &pat_recv = it_recv->second;

        // ---------------- collect send/recv items ----------------
        struct SendItem
        {
            const HaloRegion *r;
            int fid;
            int this_block;
            Box3 sb;
            double factor = 1.0;
            int32_t len = 0;
        };
        struct RecvItem
        {
            const HaloRegion *r;
            CouplingBufferBlock *buf;
            Box3 rb;
            int32_t len = 0;
        };

        std::vector<SendItem> sends;
        std::vector<RecvItem> recvs;

        std::vector<const HaloPattern *> send_patterns;
        if (coupling_channel_needs_form_transfer_(ch))
        {
            for (int axis = 0; axis < 3; ++axis)
            {
                CouplingPatternKey skey{src, dst,
                                        coupling_form_location_from_axis_(ch.value_kind, axis),
                                        ch.nghost};
                auto sit = coupling_parallel_edge_patterns_send.find(skey);
                if (sit != coupling_parallel_edge_patterns_send.end())
                    send_patterns.push_back(&sit->second);
            }
        }
        else
        {
            send_patterns.push_back(&pat_send);
        }

        // send side regions（src -> dst）
        for (const HaloPattern *src_pat : send_patterns)
        for (const HaloRegion &r : src_pat->regions)
        {
            int src_axis = -1;
            std::string tag = ch.tag;
            double factor = 1.0;
            if (coupling_channel_needs_form_transfer_(ch))
            {
                const int dst_axis = coupling_form_axis_from_location_(ch.location);
                src_axis = coupling_src_axis_from_src_to_dst_transform_(r.trans, dst_axis);
                if (coupling_form_location_from_axis_(ch.value_kind, src_axis) != src_pat->location)
                    continue;
                tag = find_triplet_field_name_(ch.tag, ch.value_kind, src_axis);
                factor = static_cast<double>(
                    coupling_form_orientation_sign_src_to_dst_(ch, r.trans, src_axis));
            }

            // send_box 为空就跳过
            const Box3 &sb = r.send_box;
            const int Ni = sb.hi.i - sb.lo.i;
            const int Nj = sb.hi.j - sb.lo.j;
            const int Nk = sb.hi.k - sb.lo.k;
            if (Ni <= 0 || Nj <= 0 || Nk <= 0)
                continue;

            const int fid = fld_->field_id(tag);

            // 发送侧取 src field 数据；若该块未分配该 field，可选择跳过或 fatal
            if (!field_active_(fid, r.this_block))
                continue;

            SendItem it;
            it.r = &r;
            it.fid = fid;
            it.this_block = r.this_block;
            it.sb = sb;
            it.factor = factor;
            sends.push_back(std::move(it));
        }

        // recv side regions（dst 接收 src）
        for (const HaloRegion &r : pat_recv.regions)
        {
            const Box3 &rb = r.recv_box;
            const int Ni = rb.hi.i - rb.lo.i;
            const int Nj = rb.hi.j - rb.lo.j;
            const int Nk = rb.hi.k - rb.lo.k;
            if (Ni <= 0 || Nj <= 0 || Nk <= 0)
                continue;

            // 在 cb.parallel_edge[cid][*] 里找对应 buffer（用 this_block + box 精确匹配）
            CouplingBufferBlock *found = nullptr;
            if (cid < cb.parallel_edge.size())
            {
                for (auto &blk : cb.parallel_edge[cid])
                {
                    if (!blk.allocated)
                        continue;
                    if (blk.this_block != r.this_block)
                        continue;
                    if (blk.location != ch.location)
                        continue;
                    if (blk.nghost != ch.nghost)
                        continue;
                    if (!HALO_TOOLS::box_equal(blk.box, rb))
                        continue;

                    found = &blk;
                    break;
                }
            }
            if (found == nullptr)
            {
                std::cout << "Fatal: coupling_parallel_edge cannot find matching CouplingBufferBlock "
                          << "for pair " << src << " -> " << dst
                          << " (cid=" << cid << ", loc=" << (int)ch.location
                          << ", nghost=" << ch.nghost << ")\n";
                std::exit(-1);
            }

            RecvItem it;
            it.r = &r;
            it.buf = found;
            it.rb = rb;
            recvs.push_back(std::move(it));
        }

        // ---------------- ensure reusable buffers ----------------
        if (send_buf.size() < sends.size())
            send_buf.resize(sends.size());
        if (req_send.size() < sends.size())
            req_send.resize(sends.size());
        if (stat_send.size() < sends.size())
            stat_send.resize(sends.size());
        if (length.size() < sends.size())
            length.resize(sends.size());

        if (recv_buf.size() < recvs.size())
            recv_buf.resize(recvs.size());
        if (req_recv.size() < recvs.size())
            req_recv.resize(recvs.size());
        if (stat_recv.size() < recvs.size())
            stat_recv.resize(recvs.size());
        if (length_corner_recv.size() < recvs.size())
            length_corner_recv.resize(recvs.size());

        // ---------------- pack sends (neighbor order) ----------------
        for (size_t is = 0; is < sends.size(); ++is)
        {
            SendItem &it = sends[is];
            FieldBlock &fb = fld_->field(it.fid, it.this_block);

            // region.trans 语义：this(src) -> nb(dst)
            pack_to_neighbor_order_scaled(
                fb, it.sb, ncomp, it.r->trans, it.factor, send_buf[is]);

            it.len = (int32_t)send_buf[is].size();
            length[is] = it.len;
        }

        // ---------------- prepare recv buffers ----------------
        for (size_t ir = 0; ir < recvs.size(); ++ir)
        {
            RecvItem &it = recvs[ir];
            const Box3 &rb = it.rb;
            const int32_t n_total =
                (rb.hi.i - rb.lo.i) * (rb.hi.j - rb.lo.j) * (rb.hi.k - rb.lo.k) * ncomp;

            it.len = n_total;
            length_corner_recv[ir] = n_total;

            if ((int32_t)recv_buf[ir].size() < n_total)
                recv_buf[ir].resize(n_total);
        }

        // ---------------- MPI (recv first, then send) ----------------
        for (size_t ir = 0; ir < recvs.size(); ++ir)
        {
            const HaloRegion *r = recvs[ir].r;
            PARALLEL::mpi_data_recv(
                r->neighbor_rank, r->recv_flag,
                recv_buf[ir].data(), length_corner_recv[ir],
                &(req_recv[ir]));
        }

        for (size_t is = 0; is < sends.size(); ++is)
        {
            const HaloRegion *r = sends[is].r;
            PARALLEL::mpi_data_send(
                r->neighbor_rank, r->send_flag,
                send_buf[is].data(), length[is],
                &(req_send[is]));
        }

        if (!sends.empty())
        {
            int nsend = (int)sends.size();
            PARALLEL::mpi_wait(nsend, req_send.data(), stat_send.data());
        }
        if (!recvs.empty())
        {
            int nrecv = (int)recvs.size();
            PARALLEL::mpi_wait(nrecv, req_recv.data(), stat_recv.data());
        }

        // 严格流水线
        PARALLEL::mpi_barrier();

        // ---------------- unpack to coupling buffers ----------------
        for (size_t ir = 0; ir < recvs.size(); ++ir)
        {
            RecvItem &it = recvs[ir];
            HALO_TOOLS::unpack_to_coupling_buffer(*(it.buf), it.rb, ncomp, recv_buf[ir]);
        }

        PARALLEL::mpi_barrier();
    }
}

void Halo::coupling_parallel_edge(const std::string &src, const std::string &dst, std::vector<int32_t> &field_cids)
{
    // 2D+ 才有 edge corner
    const int dim = fld_->grd->dimension;
    if (dim < 2)
        return;

    if (!fld_->has_coupling_pair(src, dst))
        return;

    CouplingBuffersForPair &cb = fld_->coupling_buffers(src, dst);
    const CouplingPairDesc &desc = fld_->coupling_pair(src, dst);

    if (cb.parallel_edge.empty())
        return;

    // 强制串行：按 channel 逐个 stage
    for (size_t cid : field_cids)
    {
        const CouplingChannelSpec &ch = desc.channels[cid];
        const int ncomp = ch.ncomp;

        // ---- 若 topo 根本没有该 pair 的 parallel coupling edge patch，就跳过（避免无谓 fatal）----
        bool has_any_patch = false;
        const auto &parallel_edges = TOPO_VIEW::edge_patches(*topo_, TOPO::PatchKind::Parallel);
        for (const auto &ep : parallel_edges)
        {
            if (!ep.is_coupling)
                continue;
            if ((ep.this_block_name == dst && ep.nb_block_name == src) ||
                (ep.this_block_name == src && ep.nb_block_name == dst))
            {
                has_any_patch = true;
                break;
            }
        }
        if (!has_any_patch)
            continue;

        // ---- pattern 取出（必须已在 build 阶段建好并缓存）----
        CouplingPatternKey pkey{src, dst, ch.location, ch.nghost};

        auto it_send = coupling_parallel_edge_patterns_send.find(pkey);
        auto it_recv = coupling_parallel_edge_patterns_recv.find(pkey);

        if (it_send == coupling_parallel_edge_patterns_send.end() ||
            it_recv == coupling_parallel_edge_patterns_recv.end())
        {
            std::cout << "Fatal: coupling_parallel_edge pattern not built for pair "
                      << src << " -> " << dst << " (loc=" << (int)ch.location
                      << ", nghost=" << ch.nghost << ")\n";
            std::exit(-1);
        }

        const HaloPattern &pat_send = it_send->second;
        const HaloPattern &pat_recv = it_recv->second;

        // ---------------- collect send/recv items ----------------
        struct SendItem
        {
            const HaloRegion *r;
            int fid;
            int this_block;
            Box3 sb;
            double factor = 1.0;
            int32_t len = 0;
        };
        struct RecvItem
        {
            const HaloRegion *r;
            CouplingBufferBlock *buf;
            Box3 rb;
            int32_t len = 0;
        };

        std::vector<SendItem> sends;
        std::vector<RecvItem> recvs;

        std::vector<const HaloPattern *> send_patterns;
        if (coupling_channel_needs_form_transfer_(ch))
        {
            for (int axis = 0; axis < 3; ++axis)
            {
                CouplingPatternKey skey{src, dst,
                                        coupling_form_location_from_axis_(ch.value_kind, axis),
                                        ch.nghost};
                auto sit = coupling_parallel_edge_patterns_send.find(skey);
                if (sit != coupling_parallel_edge_patterns_send.end())
                    send_patterns.push_back(&sit->second);
            }
        }
        else
        {
            send_patterns.push_back(&pat_send);
        }

        // send side regions（src -> dst）
        for (const HaloPattern *src_pat : send_patterns)
        for (const HaloRegion &r : src_pat->regions)
        {
            int src_axis = -1;
            std::string tag = ch.tag;
            double factor = 1.0;
            if (coupling_channel_needs_form_transfer_(ch))
            {
                const int dst_axis = coupling_form_axis_from_location_(ch.location);
                src_axis = coupling_src_axis_from_src_to_dst_transform_(r.trans, dst_axis);
                if (coupling_form_location_from_axis_(ch.value_kind, src_axis) != src_pat->location)
                    continue;
                tag = find_triplet_field_name_(ch.tag, ch.value_kind, src_axis);
                factor = static_cast<double>(
                    coupling_form_orientation_sign_src_to_dst_(ch, r.trans, src_axis));
            }

            // send_box 为空就跳过
            const Box3 &sb = r.send_box;
            const int Ni = sb.hi.i - sb.lo.i;
            const int Nj = sb.hi.j - sb.lo.j;
            const int Nk = sb.hi.k - sb.lo.k;
            if (Ni <= 0 || Nj <= 0 || Nk <= 0)
                continue;

            const int fid = fld_->field_id(tag);

            // 发送侧取 src field 数据；若该块未分配该 field，可选择跳过或 fatal
            if (!field_active_(fid, r.this_block))
                continue;

            SendItem it;
            it.r = &r;
            it.fid = fid;
            it.this_block = r.this_block;
            it.sb = sb;
            it.factor = factor;
            sends.push_back(std::move(it));
        }

        // recv side regions（dst 接收 src）
        for (const HaloRegion &r : pat_recv.regions)
        {
            const Box3 &rb = r.recv_box;
            const int Ni = rb.hi.i - rb.lo.i;
            const int Nj = rb.hi.j - rb.lo.j;
            const int Nk = rb.hi.k - rb.lo.k;
            if (Ni <= 0 || Nj <= 0 || Nk <= 0)
                continue;

            // 在 cb.parallel_edge[cid][*] 里找对应 buffer（用 this_block + box 精确匹配）
            CouplingBufferBlock *found = nullptr;
            if (cid < cb.parallel_edge.size())
            {
                for (auto &blk : cb.parallel_edge[cid])
                {
                    if (!blk.allocated)
                        continue;
                    if (blk.this_block != r.this_block)
                        continue;
                    if (blk.location != ch.location)
                        continue;
                    if (blk.nghost != ch.nghost)
                        continue;
                    if (!HALO_TOOLS::box_equal(blk.box, rb))
                        continue;

                    found = &blk;
                    break;
                }
            }
            if (found == nullptr)
            {
                std::cout << "Fatal: coupling_parallel_edge cannot find matching CouplingBufferBlock "
                          << "for pair " << src << " -> " << dst
                          << " (cid=" << cid << ", loc=" << (int)ch.location
                          << ", nghost=" << ch.nghost << ")\n";
                std::exit(-1);
            }

            RecvItem it;
            it.r = &r;
            it.buf = found;
            it.rb = rb;
            recvs.push_back(std::move(it));
        }

        // ---------------- ensure reusable buffers ----------------
        if (send_buf.size() < sends.size())
            send_buf.resize(sends.size());
        if (req_send.size() < sends.size())
            req_send.resize(sends.size());
        if (stat_send.size() < sends.size())
            stat_send.resize(sends.size());
        if (length.size() < sends.size())
            length.resize(sends.size());

        if (recv_buf.size() < recvs.size())
            recv_buf.resize(recvs.size());
        if (req_recv.size() < recvs.size())
            req_recv.resize(recvs.size());
        if (stat_recv.size() < recvs.size())
            stat_recv.resize(recvs.size());
        if (length_corner_recv.size() < recvs.size())
            length_corner_recv.resize(recvs.size());

        // ---------------- pack sends (neighbor order) ----------------
        for (size_t is = 0; is < sends.size(); ++is)
        {
            SendItem &it = sends[is];
            FieldBlock &fb = fld_->field(it.fid, it.this_block);

            // region.trans 语义：this(src) -> nb(dst)
            pack_to_neighbor_order_scaled(
                fb, it.sb, ncomp, it.r->trans, it.factor, send_buf[is]);

            it.len = (int32_t)send_buf[is].size();
            length[is] = it.len;
        }

        // ---------------- prepare recv buffers ----------------
        for (size_t ir = 0; ir < recvs.size(); ++ir)
        {
            RecvItem &it = recvs[ir];
            const Box3 &rb = it.rb;
            const int32_t n_total =
                (rb.hi.i - rb.lo.i) * (rb.hi.j - rb.lo.j) * (rb.hi.k - rb.lo.k) * ncomp;

            it.len = n_total;
            length_corner_recv[ir] = n_total;

            if ((int32_t)recv_buf[ir].size() < n_total)
                recv_buf[ir].resize(n_total);
        }

        // ---------------- MPI (recv first, then send) ----------------
        for (size_t ir = 0; ir < recvs.size(); ++ir)
        {
            const HaloRegion *r = recvs[ir].r;
            PARALLEL::mpi_data_recv(
                r->neighbor_rank, r->recv_flag,
                recv_buf[ir].data(), length_corner_recv[ir],
                &(req_recv[ir]));
        }

        for (size_t is = 0; is < sends.size(); ++is)
        {
            const HaloRegion *r = sends[is].r;
            PARALLEL::mpi_data_send(
                r->neighbor_rank, r->send_flag,
                send_buf[is].data(), length[is],
                &(req_send[is]));
        }

        if (!sends.empty())
        {
            int nsend = (int)sends.size();
            PARALLEL::mpi_wait(nsend, req_send.data(), stat_send.data());
        }
        if (!recvs.empty())
        {
            int nrecv = (int)recvs.size();
            PARALLEL::mpi_wait(nrecv, req_recv.data(), stat_recv.data());
        }

        // 严格流水线
        PARALLEL::mpi_barrier();

        // ---------------- unpack to coupling buffers ----------------
        for (size_t ir = 0; ir < recvs.size(); ++ir)
        {
            RecvItem &it = recvs[ir];
            HALO_TOOLS::unpack_to_coupling_buffer(*(it.buf), it.rb, ncomp, recv_buf[ir]);
        }

        PARALLEL::mpi_barrier();
    }
}

void Halo::coupling_parallel_vertex(const std::string &src, const std::string &dst)
{
    // 3D 才有 vertex corner
    const int dim = fld_->grd->dimension;
    if (dim < 3)
        return;

    if (!fld_->has_coupling_pair(src, dst))
        return;

    CouplingBuffersForPair &cb = fld_->coupling_buffers(src, dst);
    const CouplingPairDesc &desc = fld_->coupling_pair(src, dst);

    if (cb.parallel_vertex.empty())
        return;

    // 强制串行：按 channel 逐个 stage
    for (size_t cid = 0; cid < desc.channels.size(); ++cid)
    {
        const CouplingChannelSpec &ch = desc.channels[cid];
        const int ncomp = ch.ncomp;

        bool has_any_patch = false;
        const auto &parallel_vertices = TOPO_VIEW::vertex_patches(*topo_, TOPO::PatchKind::Parallel);
        for (const auto &vp : parallel_vertices)
        {
            if (!vp.is_coupling)
                continue;
            if ((vp.this_block_name == dst && vp.nb_block_name == src) ||
                (vp.this_block_name == src && vp.nb_block_name == dst))
            {
                has_any_patch = true;
                break;
            }
        }
        if (!has_any_patch)
            continue;

        CouplingPatternKey pkey{src, dst, ch.location, ch.nghost};

        auto it_send = coupling_parallel_vertex_patterns_send.find(pkey);
        auto it_recv = coupling_parallel_vertex_patterns_recv.find(pkey);

        if (it_send == coupling_parallel_vertex_patterns_send.end() ||
            it_recv == coupling_parallel_vertex_patterns_recv.end())
        {
            std::cout << "Fatal: coupling_parallel_vertex pattern not built for pair "
                      << src << " -> " << dst << " (loc=" << (int)ch.location
                      << ", nghost=" << ch.nghost << ")\n";
            std::exit(-1);
        }

        const HaloPattern &pat_send = it_send->second;
        const HaloPattern &pat_recv = it_recv->second;

        struct SendItem
        {
            const HaloRegion *r;
            int fid;
            int this_block;
            Box3 sb;
            double factor = 1.0;
            int32_t len = 0;
        };
        struct RecvItem
        {
            const HaloRegion *r;
            CouplingBufferBlock *buf;
            Box3 rb;
            int32_t len = 0;
        };

        std::vector<SendItem> sends;
        std::vector<RecvItem> recvs;

        std::vector<const HaloPattern *> send_patterns;
        if (coupling_channel_needs_form_transfer_(ch))
        {
            for (int axis = 0; axis < 3; ++axis)
            {
                CouplingPatternKey skey{src, dst,
                                        coupling_form_location_from_axis_(ch.value_kind, axis),
                                        ch.nghost};
                auto sit = coupling_parallel_vertex_patterns_send.find(skey);
                if (sit != coupling_parallel_vertex_patterns_send.end())
                    send_patterns.push_back(&sit->second);
            }
        }
        else
        {
            send_patterns.push_back(&pat_send);
        }

        // send regions
        for (const HaloPattern *src_pat : send_patterns)
        for (const HaloRegion &r : src_pat->regions)
        {
            int src_axis = -1;
            std::string tag = ch.tag;
            double factor = 1.0;
            if (coupling_channel_needs_form_transfer_(ch))
            {
                const int dst_axis = coupling_form_axis_from_location_(ch.location);
                src_axis = coupling_src_axis_from_src_to_dst_transform_(r.trans, dst_axis);
                if (coupling_form_location_from_axis_(ch.value_kind, src_axis) != src_pat->location)
                    continue;
                tag = find_triplet_field_name_(ch.tag, ch.value_kind, src_axis);
                factor = static_cast<double>(
                    coupling_form_orientation_sign_src_to_dst_(ch, r.trans, src_axis));
            }

            const Box3 &sb = r.send_box;
            const int Ni = sb.hi.i - sb.lo.i;
            const int Nj = sb.hi.j - sb.lo.j;
            const int Nk = sb.hi.k - sb.lo.k;
            if (Ni <= 0 || Nj <= 0 || Nk <= 0)
                continue;

            const int fid = fld_->field_id(tag);

            if (!field_active_(fid, r.this_block))
                continue;

            SendItem it;
            it.r = &r;
            it.fid = fid;
            it.this_block = r.this_block;
            it.sb = sb;
            it.factor = factor;
            sends.push_back(std::move(it));
        }

        // recv regions
        for (const HaloRegion &r : pat_recv.regions)
        {
            const Box3 &rb = r.recv_box;
            const int Ni = rb.hi.i - rb.lo.i;
            const int Nj = rb.hi.j - rb.lo.j;
            const int Nk = rb.hi.k - rb.lo.k;
            if (Ni <= 0 || Nj <= 0 || Nk <= 0)
                continue;

            CouplingBufferBlock *found = nullptr;
            if (cid < cb.parallel_vertex.size())
            {
                for (auto &blk : cb.parallel_vertex[cid])
                {
                    if (!blk.allocated)
                        continue;
                    if (blk.this_block != r.this_block)
                        continue;
                    if (blk.location != ch.location)
                        continue;
                    if (blk.nghost != ch.nghost)
                        continue;
                    if (!HALO_TOOLS::box_equal(blk.box, rb))
                        continue;

                    found = &blk;
                    break;
                }
            }
            if (found == nullptr)
            {
                std::cout << "Fatal: coupling_parallel_vertex cannot find matching CouplingBufferBlock "
                          << "for pair " << src << " -> " << dst
                          << " (cid=" << cid << ", loc=" << (int)ch.location
                          << ", nghost=" << ch.nghost << ")\n";
                std::exit(-1);
            }

            RecvItem it;
            it.r = &r;
            it.buf = found;
            it.rb = rb;
            recvs.push_back(std::move(it));
        }

        // buffers
        if (send_buf.size() < sends.size())
            send_buf.resize(sends.size());
        if (req_send.size() < sends.size())
            req_send.resize(sends.size());
        if (stat_send.size() < sends.size())
            stat_send.resize(sends.size());
        if (length.size() < sends.size())
            length.resize(sends.size());

        if (recv_buf.size() < recvs.size())
            recv_buf.resize(recvs.size());
        if (req_recv.size() < recvs.size())
            req_recv.resize(recvs.size());
        if (stat_recv.size() < recvs.size())
            stat_recv.resize(recvs.size());
        if (length_corner_recv.size() < recvs.size())
            length_corner_recv.resize(recvs.size());

        // pack
        for (size_t is = 0; is < sends.size(); ++is)
        {
            SendItem &it = sends[is];
            FieldBlock &fb = fld_->field(it.fid, it.this_block);

            pack_to_neighbor_order_scaled(
                fb, it.sb, ncomp, it.r->trans, it.factor, send_buf[is]);

            it.len = (int32_t)send_buf[is].size();
            length[is] = it.len;
        }

        // recv size
        for (size_t ir = 0; ir < recvs.size(); ++ir)
        {
            RecvItem &it = recvs[ir];
            const Box3 &rb = it.rb;
            const int32_t n_total =
                (rb.hi.i - rb.lo.i) * (rb.hi.j - rb.lo.j) * (rb.hi.k - rb.lo.k) * ncomp;

            it.len = n_total;
            length_corner_recv[ir] = n_total;

            if ((int32_t)recv_buf[ir].size() < n_total)
                recv_buf[ir].resize(n_total);
        }

        // MPI
        for (size_t ir = 0; ir < recvs.size(); ++ir)
        {
            const HaloRegion *r = recvs[ir].r;
            PARALLEL::mpi_data_recv(
                r->neighbor_rank, r->recv_flag,
                recv_buf[ir].data(), length_corner_recv[ir],
                &(req_recv[ir]));
        }

        for (size_t is = 0; is < sends.size(); ++is)
        {
            const HaloRegion *r = sends[is].r;
            PARALLEL::mpi_data_send(
                r->neighbor_rank, r->send_flag,
                send_buf[is].data(), length[is],
                &(req_send[is]));
        }

        if (!sends.empty())
        {
            int nsend = (int)sends.size();
            PARALLEL::mpi_wait(nsend, req_send.data(), stat_send.data());
        }
        if (!recvs.empty())
        {
            int nrecv = (int)recvs.size();
            PARALLEL::mpi_wait(nrecv, req_recv.data(), stat_recv.data());
        }

        PARALLEL::mpi_barrier();

        // unpack
        for (size_t ir = 0; ir < recvs.size(); ++ir)
        {
            RecvItem &it = recvs[ir];
            HALO_TOOLS::unpack_to_coupling_buffer(*(it.buf), it.rb, ncomp, recv_buf[ir]);
        }

        PARALLEL::mpi_barrier();
    }
}

void Halo::coupling_parallel_vertex(const std::string &src, const std::string &dst, std::vector<int32_t> &field_cids)
{
    // 3D 才有 vertex corner
    const int dim = fld_->grd->dimension;
    if (dim < 3)
        return;

    if (!fld_->has_coupling_pair(src, dst))
        return;

    CouplingBuffersForPair &cb = fld_->coupling_buffers(src, dst);
    const CouplingPairDesc &desc = fld_->coupling_pair(src, dst);

    if (cb.parallel_vertex.empty())
        return;

    // 强制串行：按 channel 逐个 stage
    for (auto cid : field_cids)
    {
        const CouplingChannelSpec &ch = desc.channels[cid];
        const int ncomp = ch.ncomp;

        bool has_any_patch = false;
        const auto &parallel_vertices = TOPO_VIEW::vertex_patches(*topo_, TOPO::PatchKind::Parallel);
        for (const auto &vp : parallel_vertices)
        {
            if (!vp.is_coupling)
                continue;
            if ((vp.this_block_name == dst && vp.nb_block_name == src) ||
                (vp.this_block_name == src && vp.nb_block_name == dst))
            {
                has_any_patch = true;
                break;
            }
        }
        if (!has_any_patch)
            continue;

        CouplingPatternKey pkey{src, dst, ch.location, ch.nghost};

        auto it_send = coupling_parallel_vertex_patterns_send.find(pkey);
        auto it_recv = coupling_parallel_vertex_patterns_recv.find(pkey);

        if (it_send == coupling_parallel_vertex_patterns_send.end() ||
            it_recv == coupling_parallel_vertex_patterns_recv.end())
        {
            std::cout << "Fatal: coupling_parallel_vertex pattern not built for pair "
                      << src << " -> " << dst << " (loc=" << (int)ch.location
                      << ", nghost=" << ch.nghost << ")\n";
            std::exit(-1);
        }

        const HaloPattern &pat_send = it_send->second;
        const HaloPattern &pat_recv = it_recv->second;

        struct SendItem
        {
            const HaloRegion *r;
            int fid;
            int this_block;
            Box3 sb;
            double factor = 1.0;
            int32_t len = 0;
        };
        struct RecvItem
        {
            const HaloRegion *r;
            CouplingBufferBlock *buf;
            Box3 rb;
            int32_t len = 0;
        };

        std::vector<SendItem> sends;
        std::vector<RecvItem> recvs;

        std::vector<const HaloPattern *> send_patterns;
        if (coupling_channel_needs_form_transfer_(ch))
        {
            for (int axis = 0; axis < 3; ++axis)
            {
                CouplingPatternKey skey{src, dst,
                                        coupling_form_location_from_axis_(ch.value_kind, axis),
                                        ch.nghost};
                auto sit = coupling_parallel_vertex_patterns_send.find(skey);
                if (sit != coupling_parallel_vertex_patterns_send.end())
                    send_patterns.push_back(&sit->second);
            }
        }
        else
        {
            send_patterns.push_back(&pat_send);
        }

        // send regions
        for (const HaloPattern *src_pat : send_patterns)
        for (const HaloRegion &r : src_pat->regions)
        {
            int src_axis = -1;
            std::string tag = ch.tag;
            double factor = 1.0;
            if (coupling_channel_needs_form_transfer_(ch))
            {
                const int dst_axis = coupling_form_axis_from_location_(ch.location);
                src_axis = coupling_src_axis_from_src_to_dst_transform_(r.trans, dst_axis);
                if (coupling_form_location_from_axis_(ch.value_kind, src_axis) != src_pat->location)
                    continue;
                tag = find_triplet_field_name_(ch.tag, ch.value_kind, src_axis);
                factor = static_cast<double>(
                    coupling_form_orientation_sign_src_to_dst_(ch, r.trans, src_axis));
            }

            const Box3 &sb = r.send_box;
            const int Ni = sb.hi.i - sb.lo.i;
            const int Nj = sb.hi.j - sb.lo.j;
            const int Nk = sb.hi.k - sb.lo.k;
            if (Ni <= 0 || Nj <= 0 || Nk <= 0)
                continue;

            const int fid = fld_->field_id(tag);

            if (!field_active_(fid, r.this_block))
                continue;

            SendItem it;
            it.r = &r;
            it.fid = fid;
            it.this_block = r.this_block;
            it.sb = sb;
            it.factor = factor;
            sends.push_back(std::move(it));
        }

        // recv regions
        for (const HaloRegion &r : pat_recv.regions)
        {
            const Box3 &rb = r.recv_box;
            const int Ni = rb.hi.i - rb.lo.i;
            const int Nj = rb.hi.j - rb.lo.j;
            const int Nk = rb.hi.k - rb.lo.k;
            if (Ni <= 0 || Nj <= 0 || Nk <= 0)
                continue;

            CouplingBufferBlock *found = nullptr;
            if (cid < cb.parallel_vertex.size())
            {
                for (auto &blk : cb.parallel_vertex[cid])
                {
                    if (!blk.allocated)
                        continue;
                    if (blk.this_block != r.this_block)
                        continue;
                    if (blk.location != ch.location)
                        continue;
                    if (blk.nghost != ch.nghost)
                        continue;
                    if (!HALO_TOOLS::box_equal(blk.box, rb))
                        continue;

                    found = &blk;
                    break;
                }
            }
            if (found == nullptr)
            {
                std::cout << "Fatal: coupling_parallel_vertex cannot find matching CouplingBufferBlock "
                          << "for pair " << src << " -> " << dst
                          << " (cid=" << cid << ", loc=" << (int)ch.location
                          << ", nghost=" << ch.nghost << ")\n";
                std::exit(-1);
            }

            RecvItem it;
            it.r = &r;
            it.buf = found;
            it.rb = rb;
            recvs.push_back(std::move(it));
        }

        // buffers
        if (send_buf.size() < sends.size())
            send_buf.resize(sends.size());
        if (req_send.size() < sends.size())
            req_send.resize(sends.size());
        if (stat_send.size() < sends.size())
            stat_send.resize(sends.size());
        if (length.size() < sends.size())
            length.resize(sends.size());

        if (recv_buf.size() < recvs.size())
            recv_buf.resize(recvs.size());
        if (req_recv.size() < recvs.size())
            req_recv.resize(recvs.size());
        if (stat_recv.size() < recvs.size())
            stat_recv.resize(recvs.size());
        if (length_corner_recv.size() < recvs.size())
            length_corner_recv.resize(recvs.size());

        // pack
        for (size_t is = 0; is < sends.size(); ++is)
        {
            SendItem &it = sends[is];
            FieldBlock &fb = fld_->field(it.fid, it.this_block);

            pack_to_neighbor_order_scaled(
                fb, it.sb, ncomp, it.r->trans, it.factor, send_buf[is]);

            it.len = (int32_t)send_buf[is].size();
            length[is] = it.len;
        }

        // recv size
        for (size_t ir = 0; ir < recvs.size(); ++ir)
        {
            RecvItem &it = recvs[ir];
            const Box3 &rb = it.rb;
            const int32_t n_total =
                (rb.hi.i - rb.lo.i) * (rb.hi.j - rb.lo.j) * (rb.hi.k - rb.lo.k) * ncomp;

            it.len = n_total;
            length_corner_recv[ir] = n_total;

            if ((int32_t)recv_buf[ir].size() < n_total)
                recv_buf[ir].resize(n_total);
        }

        // MPI
        for (size_t ir = 0; ir < recvs.size(); ++ir)
        {
            const HaloRegion *r = recvs[ir].r;
            PARALLEL::mpi_data_recv(
                r->neighbor_rank, r->recv_flag,
                recv_buf[ir].data(), length_corner_recv[ir],
                &(req_recv[ir]));
        }

        for (size_t is = 0; is < sends.size(); ++is)
        {
            const HaloRegion *r = sends[is].r;
            PARALLEL::mpi_data_send(
                r->neighbor_rank, r->send_flag,
                send_buf[is].data(), length[is],
                &(req_send[is]));
        }

        if (!sends.empty())
        {
            int nsend = (int)sends.size();
            PARALLEL::mpi_wait(nsend, req_send.data(), stat_send.data());
        }
        if (!recvs.empty())
        {
            int nrecv = (int)recvs.size();
            PARALLEL::mpi_wait(nrecv, req_recv.data(), stat_recv.data());
        }

        PARALLEL::mpi_barrier();

        // unpack
        for (size_t ir = 0; ir < recvs.size(); ++ir)
        {
            RecvItem &it = recvs[ir];
            HALO_TOOLS::unpack_to_coupling_buffer(*(it.buf), it.rb, ncomp, recv_buf[ir]);
        }

        PARALLEL::mpi_barrier();
    }
}
