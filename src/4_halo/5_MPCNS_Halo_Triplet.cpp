#include "4_halo/1_MPCNS_Halo.h"
#include "2_topology/TopologyView.h"
#include "0_basic/MPI_WRAPPER.h"
#include "4_halo/detail/halo_build_boxmakers.h"
#include "4_halo/detail/halo_build_tools.h"

#include <array>
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>

namespace
{
StaggerLocation edge_loc_from_axis(int axis)
{
    if (axis == 0)
        return StaggerLocation::EdgeXi;
    if (axis == 1)
        return StaggerLocation::EdgeEt;
    return StaggerLocation::EdgeZe;
}

int inverse_perm_axis(const TOPO::IndexTransform &tr, int dst_axis)
{
    for (int a = 0; a < 3; ++a)
        if (tr.perm[a] == dst_axis)
            return a;

    std::cout << "Fatal Error!!! invalid IndexTransform permutation in triplet halo.\n";
    std::exit(-1);
}

void map_index(const TOPO::IndexTransform &tr,
               int i, int j, int k,
               int &io, int &jo, int &ko)
{
    int loc[3] = {i, j, k};
    int tar[3] = {0, 0, 0};
    int offset[3] = {tr.offset.i, tr.offset.j, tr.offset.k};

    for (int a = 0; a < 3; ++a)
        tar[tr.perm[a]] = tr.sign[a] * loc[a] + offset[a];

    io = tar[0];
    jo = tar[1];
    ko = tar[2];
}

std::array<int, 3> edge_triplet_ids(Field *fld, const std::vector<std::string> &fields)
{
    if (fields.size() != 3)
    {
        std::cout << "Fatal Error!!! edge triplet halo expects {xi, eta, zeta}.\n";
        std::exit(-1);
    }

    std::array<int, 3> fid = {
        fld->field_id(fields[0]),
        fld->field_id(fields[1]),
        fld->field_id(fields[2])};

    const StaggerLocation expect[3] = {
        StaggerLocation::EdgeXi,
        StaggerLocation::EdgeEt,
        StaggerLocation::EdgeZe};

    const int ncomp = fld->descriptor(fid[0]).ncomp;
    const int nghost = fld->descriptor(fid[0]).nghost;
    for (int a = 0; a < 3; ++a)
    {
        const auto &d = fld->descriptor(fid[a]);
        if (d.location != expect[a] || d.ncomp != ncomp || d.nghost != nghost)
        {
            std::cout << "Fatal Error!!! invalid edge triplet fields: "
                      << fields[0] << ", " << fields[1] << ", " << fields[2] << "\n";
            std::exit(-1);
        }
    }

    return fid;
}

int triplet_ncomp(Field *fld, const std::array<int, 3> &fid)
{
    return fld->descriptor(fid[0]).ncomp;
}

int triplet_nghost(Field *fld, const std::array<int, 3> &fid)
{
    return fld->descriptor(fid[0]).nghost;
}

struct TripletCornerMeta
{
    int recv_rank;
    int send_rank;
    int recv_block;
    int send_block;

    int recv_axis;
    int send_axis;
    int sign;

    Box3 recv_box;
    Box3 node_box_on_send;
    int dir1_send;
    int dir2_send;
    int dir3_send;

    TOPO::IndexTransform trans_recv_to_send;
    int tag;
};

void exchange_triplet_corner_meta(const std::map<int, std::vector<TripletCornerMeta>> &meta_to_send,
                                  std::vector<TripletCornerMeta> &recv_metas)
{
    int nrank = 1;
    PARALLEL::mpi_size(&nrank);

    std::vector<int> send_counts(nrank, 0), recv_counts(nrank, 0);
    for (const auto &kv : meta_to_send)
        send_counts[kv.first] = static_cast<int>(kv.second.size());

    PARALLEL::mpi_alltoall(send_counts.data(), 1, recv_counts.data(), 1);

    std::vector<int> sdispls(nrank, 0), rdispls(nrank, 0);
    int total_send = 0;
    int total_recv = 0;
    for (int r = 0; r < nrank; ++r)
    {
        sdispls[r] = total_send;
        rdispls[r] = total_recv;
        total_send += send_counts[r];
        total_recv += recv_counts[r];
    }

    std::vector<TripletCornerMeta> send_buf(total_send);
    recv_metas.resize(total_recv);

    for (int r = 0; r < nrank; ++r)
    {
        auto it = meta_to_send.find(r);
        if (it == meta_to_send.end())
            continue;
        std::copy(it->second.begin(), it->second.end(), send_buf.begin() + sdispls[r]);
    }

    const int sz = static_cast<int>(sizeof(TripletCornerMeta));
    std::vector<int> send_counts_bytes(nrank), recv_counts_bytes(nrank);
    std::vector<int> sdispls_bytes(nrank), rdispls_bytes(nrank);
    for (int r = 0; r < nrank; ++r)
    {
        send_counts_bytes[r] = send_counts[r] * sz;
        recv_counts_bytes[r] = recv_counts[r] * sz;
        sdispls_bytes[r] = sdispls[r] * sz;
        rdispls_bytes[r] = rdispls[r] * sz;
    }

    MPI_Alltoallv(reinterpret_cast<const void *>(send_buf.data()),
                  send_counts_bytes.data(),
                  sdispls_bytes.data(),
                  MPI_BYTE,
                  reinterpret_cast<void *>(recv_metas.data()),
                  recv_counts_bytes.data(),
                  rdispls_bytes.data(),
                  MPI_BYTE,
                  MPI_COMM_WORLD);
}

int box_volume(const Box3 &b)
{
    return (b.hi.i - b.lo.i) * (b.hi.j - b.lo.j) * (b.hi.k - b.lo.k);
}

void pack_triplet_corner_send(Field *fld,
                              const std::array<int, 3> &fid,
                              const TripletCornerMeta &m,
                              const Box3 &send_box,
                              int ncomp,
                              std::vector<double> &buf)
{
    FieldBlock &fb = fld->field(fid[m.send_axis], m.send_block);
    if (!fb.is_allocated())
        return;

    const int ni = m.recv_box.hi.i - m.recv_box.lo.i;
    const int nj = m.recv_box.hi.j - m.recv_box.lo.j;
    const int nk = m.recv_box.hi.k - m.recv_box.lo.k;
    const int n_total = ni * nj * nk * ncomp;
    if (static_cast<int>(buf.size()) < n_total)
        buf.resize(n_total);

    const TOPO::IndexTransform send_to_recv =
        HALO_TOOLS::inverse_transform(m.trans_recv_to_send);

    for (int i = send_box.lo.i; i < send_box.hi.i; ++i)
        for (int j = send_box.lo.j; j < send_box.hi.j; ++j)
            for (int k = send_box.lo.k; k < send_box.hi.k; ++k)
            {
                int ir, jr, kr;
                map_index(send_to_recv, i, j, k, ir, jr, kr);

                const int ii = ir - m.recv_box.lo.i;
                const int jj = jr - m.recv_box.lo.j;
                const int kk = kr - m.recv_box.lo.k;
                const int base = ((ii * nj + jj) * nk + kk) * ncomp;

                for (int c = 0; c < ncomp; ++c)
                    buf[base + c] = static_cast<double>(m.sign) * fb(i, j, k, c);
            }
}

void unpack_triplet_corner_recv(Field *fld,
                                const std::array<int, 3> &fid,
                                const TripletCornerMeta &m,
                                int ncomp,
                                const std::vector<double> &buf)
{
    FieldBlock &fb = fld->field(fid[m.recv_axis], m.recv_block);
    if (!fb.is_allocated())
        return;

    const Box3 &rb = m.recv_box;
    const int ni = rb.hi.i - rb.lo.i;
    const int nj = rb.hi.j - rb.lo.j;
    const int nk = rb.hi.k - rb.lo.k;

    for (int ii = 0; ii < ni; ++ii)
        for (int jj = 0; jj < nj; ++jj)
            for (int kk = 0; kk < nk; ++kk)
            {
                const int base = ((ii * nj + jj) * nk + kk) * ncomp;
                for (int c = 0; c < ncomp; ++c)
                    fb(rb.lo.i + ii, rb.lo.j + jj, rb.lo.k + kk, c) = buf[base + c];
            }
}
} // namespace

void Halo::data_trans_edge_1form_triplet(const std::vector<std::string> &fields,
                                         HaloLevel stage)
{
    if (stage == HaloLevel::FaceOnly)
    {
        exchange_inner_face_edge_1form_triplet_(fields);
        exchange_parallel_face_edge_1form_triplet_(fields);
        return;
    }

    if (stage == HaloLevel::Edge)
    {
        exchange_inner_edge_edge_1form_triplet_(fields);
        exchange_parallel_edge_edge_1form_triplet_(fields);
        return;
    }

    exchange_inner_vertex_edge_1form_triplet_(fields);
    exchange_parallel_vertex_edge_1form_triplet_(fields);
}

void Halo::exchange_inner_face_edge_1form_triplet_(const std::vector<std::string> &fields)
{
    const auto fid = edge_triplet_ids(fld_, fields);
    const int ncomp = triplet_ncomp(fld_, fid);
    const int nghost = triplet_nghost(fld_, fid);

    for (int dst_axis = 0; dst_axis < 3; ++dst_axis)
    {
        PatternKey dst_key{edge_loc_from_axis(dst_axis), nghost};
        auto it_dst = inner_patterns_.find(dst_key);
        if (it_dst == inner_patterns_.end())
            continue;

        const auto &dst_regions = it_dst->second.regions;
        for (std::size_t ir = 0; ir < dst_regions.size(); ++ir)
        {
            const HaloRegion &rdst = dst_regions[ir];
            const int src_axis = inverse_perm_axis(rdst.trans, dst_axis);
            PatternKey src_key{edge_loc_from_axis(src_axis), nghost};
            auto it_src = inner_patterns_.find(src_key);
            if (it_src == inner_patterns_.end() || ir >= it_src->second.regions.size())
                continue;

            const HaloRegion &rsrc = it_src->second.regions[ir];
            FieldBlock &fb_send = fld_->field(fid[src_axis], rsrc.this_block);
            FieldBlock &fb_recv = fld_->field(fid[dst_axis], rdst.neighbor_block);
            if (!fb_send.is_allocated() || !fb_recv.is_allocated())
                continue;

            const double factor = static_cast<double>(rdst.trans.sign[src_axis]);
            const Box3 &sb = rsrc.send_box;

            for (int i = sb.lo.i; i < sb.hi.i; ++i)
                for (int j = sb.lo.j; j < sb.hi.j; ++j)
                    for (int k = sb.lo.k; k < sb.hi.k; ++k)
                    {
                        int io, jo, ko;
                        map_index(rdst.trans, i, j, k, io, jo, ko);
                        for (int m = 0; m < ncomp; ++m)
                            fb_recv(io, jo, ko, m) = factor * fb_send(i, j, k, m);
                    }
        }
    }
}

void Halo::exchange_inner_edge_edge_1form_triplet_(const std::vector<std::string> &fields)
{
    const auto fid = edge_triplet_ids(fld_, fields);
    const int ncomp = triplet_ncomp(fld_, fid);
    const int nghost = triplet_nghost(fld_, fid);

    for (int dst_axis = 0; dst_axis < 3; ++dst_axis)
    {
        PatternKey dst_key{edge_loc_from_axis(dst_axis), nghost};
        auto it_dst = inner_edge_patterns_.find(dst_key);
        if (it_dst == inner_edge_patterns_.end())
            continue;

        for (const HaloRegion &r : it_dst->second.regions)
        {
            const int src_axis = r.trans.perm[dst_axis];
            FieldBlock &fb_recv = fld_->field(fid[dst_axis], r.this_block);
            FieldBlock &fb_send = fld_->field(fid[src_axis], r.neighbor_block);
            if (!fb_recv.is_allocated() || !fb_send.is_allocated())
                continue;

            const double factor = static_cast<double>(r.trans.sign[dst_axis]);
            const Box3 &rb = r.recv_box;
            for (int i = rb.lo.i; i < rb.hi.i; ++i)
                for (int j = rb.lo.j; j < rb.hi.j; ++j)
                    for (int k = rb.lo.k; k < rb.hi.k; ++k)
                    {
                        int is, js, ks;
                        map_index(r.trans, i, j, k, is, js, ks);
                        for (int m = 0; m < ncomp; ++m)
                            fb_recv(i, j, k, m) = factor * fb_send(is, js, ks, m);
                    }
        }
    }
}

void Halo::exchange_inner_vertex_edge_1form_triplet_(const std::vector<std::string> &fields)
{
    const auto fid = edge_triplet_ids(fld_, fields);
    const int ncomp = triplet_ncomp(fld_, fid);
    const int nghost = triplet_nghost(fld_, fid);

    for (int dst_axis = 0; dst_axis < 3; ++dst_axis)
    {
        PatternKey dst_key{edge_loc_from_axis(dst_axis), nghost};
        auto it_dst = inner_vertex_patterns_.find(dst_key);
        if (it_dst == inner_vertex_patterns_.end())
            continue;

        for (const HaloRegion &r : it_dst->second.regions)
        {
            const int src_axis = r.trans.perm[dst_axis];
            FieldBlock &fb_recv = fld_->field(fid[dst_axis], r.this_block);
            FieldBlock &fb_send = fld_->field(fid[src_axis], r.neighbor_block);
            if (!fb_recv.is_allocated() || !fb_send.is_allocated())
                continue;

            const double factor = static_cast<double>(r.trans.sign[dst_axis]);
            const Box3 &rb = r.recv_box;
            for (int i = rb.lo.i; i < rb.hi.i; ++i)
                for (int j = rb.lo.j; j < rb.hi.j; ++j)
                    for (int k = rb.lo.k; k < rb.hi.k; ++k)
                    {
                        int is, js, ks;
                        map_index(r.trans, i, j, k, is, js, ks);
                        for (int m = 0; m < ncomp; ++m)
                            fb_recv(i, j, k, m) = factor * fb_send(is, js, ks, m);
                    }
        }
    }
}

void Halo::exchange_parallel_face_edge_1form_triplet_(const std::vector<std::string> &fields)
{
    const auto fid = edge_triplet_ids(fld_, fields);
    const int ncomp = triplet_ncomp(fld_, fid);
    const int nghost = triplet_nghost(fld_, fid);

    for (int dst_axis = 0; dst_axis < 3; ++dst_axis)
    {
        PatternKey dst_key{edge_loc_from_axis(dst_axis), nghost};
        auto it_dst = parallel_patterns_.find(dst_key);
        if (it_dst == parallel_patterns_.end())
            continue;

        const HaloPattern &pat_dst = it_dst->second;
        const int nreg = static_cast<int>(pat_dst.regions.size());

        if (send_buf.size() < nreg)
            send_buf.resize(nreg);
        if (recv_buf.size() < nreg)
            recv_buf.resize(nreg);
        if (req_send.size() < nreg)
            req_send.resize(nreg);
        if (req_recv.size() < nreg)
            req_recv.resize(nreg);
        if (stat_send.size() < nreg)
            stat_send.resize(nreg);
        if (stat_recv.size() < nreg)
            stat_recv.resize(nreg);
        if (length.size() < nreg)
            length.resize(nreg);

        for (int ir = 0; ir < nreg; ++ir)
        {
            const HaloRegion &rdst = pat_dst.regions[ir];
            const int src_axis = inverse_perm_axis(rdst.trans, dst_axis);
            PatternKey src_key{edge_loc_from_axis(src_axis), nghost};
            auto it_src = parallel_patterns_.find(src_key);
            if (it_src == parallel_patterns_.end() || ir >= static_cast<int>(it_src->second.regions.size()))
                continue;

            const HaloRegion &rsrc = it_src->second.regions[ir];
            FieldBlock &fb = fld_->field(fid[src_axis], rsrc.this_block);
            const Box3 &sb = rsrc.send_box;

            const int ni = rdst.recv_box.hi.i - rdst.recv_box.lo.i;
            const int nj = rdst.recv_box.hi.j - rdst.recv_box.lo.j;
            const int nk = rdst.recv_box.hi.k - rdst.recv_box.lo.k;
            const int32_t n_total = ni * nj * nk * ncomp;
            length[ir] = n_total;
            if (send_buf[ir].size() < static_cast<std::size_t>(n_total))
                send_buf[ir].resize(n_total);
            if (recv_buf[ir].size() < static_cast<std::size_t>(n_total))
                recv_buf[ir].resize(n_total);

            const double factor = static_cast<double>(rdst.trans.sign[src_axis]);

            int loc_lo[3] = {sb.lo.i, sb.lo.j, sb.lo.k};
            int loc_hi[3] = {sb.hi.i - 1, sb.hi.j - 1, sb.hi.k - 1};
            int offset[3] = {rdst.trans.offset.i, rdst.trans.offset.j, rdst.trans.offset.k};
            int tar1[3], tar2[3], tar_ref[3];
            for (int a = 0; a < 3; ++a)
                tar1[rdst.trans.perm[a]] = rdst.trans.sign[a] * loc_lo[a] + offset[a];
            for (int a = 0; a < 3; ++a)
                tar2[rdst.trans.perm[a]] = rdst.trans.sign[a] * loc_hi[a] + offset[a];

            for (int a = 0; a < 3; ++a)
                tar_ref[a] = (tar1[a] <= tar2[a]) ? tar1[a] : tar2[a];

            int len_nb[3] = {ni, nj, nk};
            for (int i = sb.lo.i; i < sb.hi.i; ++i)
                for (int j = sb.lo.j; j < sb.hi.j; ++j)
                    for (int k = sb.lo.k; k < sb.hi.k; ++k)
                    {
                        int io, jo, ko;
                        map_index(rdst.trans, i, j, k, io, jo, ko);
                        const int ri = io - tar_ref[0];
                        const int rj = jo - tar_ref[1];
                        const int rk = ko - tar_ref[2];
                        const int32_t base = ((ri * len_nb[1] + rj) * len_nb[2] + rk) * ncomp;
                        for (int m = 0; m < ncomp; ++m)
                            send_buf[ir][base + m] = factor * fb(i, j, k, m);
                    }
        }

        for (int ir = 0; ir < nreg; ++ir)
        {
            const HaloRegion &r = pat_dst.regions[ir];
            PARALLEL::mpi_data_send(r.neighbor_rank, r.send_flag,
                                    send_buf[ir].data(), length[ir], &(req_send[ir]));
            PARALLEL::mpi_data_recv(r.neighbor_rank, r.recv_flag,
                                    recv_buf[ir].data(), length[ir], &(req_recv[ir]));
        }

        int nsend = nreg;
        int nrecv = nreg;
        PARALLEL::mpi_wait(nsend, req_send.data(), stat_send.data());
        PARALLEL::mpi_wait(nrecv, req_recv.data(), stat_recv.data());
        PARALLEL::mpi_barrier();

        for (int ir = 0; ir < nreg; ++ir)
        {
            const HaloRegion &r = pat_dst.regions[ir];
            FieldBlock &fb = fld_->field(fid[dst_axis], r.this_block);
            if (!fb.is_allocated())
                continue;

            const Box3 &rb = r.recv_box;
            const int ni = rb.hi.i - rb.lo.i;
            const int nj = rb.hi.j - rb.lo.j;
            const int nk = rb.hi.k - rb.lo.k;
            for (int ii = 0; ii < ni; ++ii)
                for (int jj = 0; jj < nj; ++jj)
                    for (int kk = 0; kk < nk; ++kk)
                    {
                        const int32_t base = ((ii * nj + jj) * nk + kk) * ncomp;
                        for (int m = 0; m < ncomp; ++m)
                            fb(rb.lo.i + ii, rb.lo.j + jj, rb.lo.k + kk, m) =
                                recv_buf[ir][base + m];
                    }
        }
    }
}

void Halo::exchange_parallel_edge_edge_1form_triplet_(const std::vector<std::string> &fields)
{
    const auto fid = edge_triplet_ids(fld_, fields);
    const int ncomp = triplet_ncomp(fld_, fid);
    const int nghost = triplet_nghost(fld_, fid);

    int myid = 0;
    PARALLEL::mpi_rank(&myid);

    std::map<int, std::vector<TripletCornerMeta>> meta_to_send;
    std::map<int, int> next_tag;
    std::vector<TripletCornerMeta> recv_plan;

    const auto &parallel_edges = TOPO_VIEW::edge_patches(*topo_, TOPO::PatchKind::Parallel);
    for (const TOPO::EdgePatch &ep : parallel_edges)
    {
        if (ep.is_coupling)
            continue;

        Direction this_d1 = HALO_TOOLS::int_to_direction(ep.dir1);
        Direction this_d2 = HALO_TOOLS::int_to_direction(ep.dir2);

        int send_dir1 = HALO_TOOLS::map_dir_to_neighbor(ep.dir1, ep.trans, true);
        int send_dir2 = HALO_TOOLS::map_dir_to_neighbor(ep.dir2, ep.trans, false);
        if (ep.trans.perm[std::abs(ep.dir1) - 1] + 1 != std::abs(send_dir1))
        {
            if (ep.trans.perm[std::abs(ep.dir1) - 1] + 1 == std::abs(send_dir2))
                std::swap(send_dir1, send_dir2);
            else
            {
                std::cout << "Fatal Error!!! triplet parallel edge direction mapping broken.\n";
                std::exit(-1);
            }
        }

        for (int recv_axis = 0; recv_axis < 3; ++recv_axis)
        {
            TripletCornerMeta m{};
            m.recv_rank = ep.this_rank;
            m.send_rank = ep.nb_rank;
            m.recv_block = ep.this_block;
            m.send_block = ep.nb_block;
            m.recv_axis = recv_axis;
            m.send_axis = ep.trans.perm[recv_axis];
            m.sign = ep.trans.sign[recv_axis];
            m.recv_box = HALO_BOX::make_2DCorner_ghost_box(
                edge_loc_from_axis(recv_axis), ep.this_box_node, this_d1, this_d2, nghost);
            m.node_box_on_send = ep.nb_box_node;
            m.dir1_send = send_dir1;
            m.dir2_send = send_dir2;
            m.dir3_send = 0;
            m.trans_recv_to_send = ep.trans;
            m.tag = next_tag[ep.nb_rank]++;

            recv_plan.push_back(m);
            meta_to_send[ep.nb_rank].push_back(m);
        }
    }

    std::vector<TripletCornerMeta> send_plan;
    exchange_triplet_corner_meta(meta_to_send, send_plan);

    const int nsend = static_cast<int>(send_plan.size());
    const int nrecv = static_cast<int>(recv_plan.size());

    std::vector<std::vector<double>> sbuf(nsend);
    std::vector<std::vector<double>> rbuf(nrecv);
    std::vector<MPI_Request> sreq(nsend), rreq(nrecv);
    std::vector<MPI_Status> sstat(nsend), rstat(nrecv);
    std::vector<int> slen(nsend, 0), rlen(nrecv, 0);

    for (int i = 0; i < nsend; ++i)
    {
        const auto &m = send_plan[i];
        if (m.send_rank != myid)
        {
            std::cout << "Fatal Error!!! triplet edge meta delivered to wrong rank.\n";
            std::exit(-1);
        }

        Direction d1 = HALO_TOOLS::int_to_direction(m.dir1_send);
        Direction d2 = HALO_TOOLS::int_to_direction(m.dir2_send);
        const Box3 send_box = HALO_BOX::make_2DCorner_innerghost_box(
            edge_loc_from_axis(m.send_axis), m.node_box_on_send, d1, d2, nghost);

        slen[i] = box_volume(m.recv_box) * ncomp;
        sbuf[i].resize(slen[i]);
        pack_triplet_corner_send(fld_, fid, m, send_box, ncomp, sbuf[i]);
    }

    for (int i = 0; i < nrecv; ++i)
    {
        rlen[i] = box_volume(recv_plan[i].recv_box) * ncomp;
        rbuf[i].resize(rlen[i]);
    }

    for (int i = 0; i < nrecv; ++i)
        PARALLEL::mpi_data_recv(recv_plan[i].send_rank, recv_plan[i].tag,
                                rbuf[i].data(), rlen[i], &rreq[i]);
    for (int i = 0; i < nsend; ++i)
        PARALLEL::mpi_data_send(send_plan[i].recv_rank, send_plan[i].tag,
                                sbuf[i].data(), slen[i], &sreq[i]);

    int wait_recv = nrecv;
    int wait_send = nsend;
    if (wait_recv > 0)
        PARALLEL::mpi_wait(wait_recv, rreq.data(), rstat.data());
    if (wait_send > 0)
        PARALLEL::mpi_wait(wait_send, sreq.data(), sstat.data());
    PARALLEL::mpi_barrier();

    for (int i = 0; i < nrecv; ++i)
        unpack_triplet_corner_recv(fld_, fid, recv_plan[i], ncomp, rbuf[i]);
}

void Halo::exchange_parallel_vertex_edge_1form_triplet_(const std::vector<std::string> &fields)
{
    const auto fid = edge_triplet_ids(fld_, fields);
    const int ncomp = triplet_ncomp(fld_, fid);
    const int nghost = triplet_nghost(fld_, fid);

    int myid = 0;
    PARALLEL::mpi_rank(&myid);

    std::map<int, std::vector<TripletCornerMeta>> meta_to_send;
    std::map<int, int> next_tag;
    std::vector<TripletCornerMeta> recv_plan;

    const auto &parallel_vertices = TOPO_VIEW::vertex_patches(*topo_, TOPO::PatchKind::Parallel);
    for (const TOPO::VertexPatch &vp : parallel_vertices)
    {
        if (vp.is_coupling)
            continue;

        Direction this_d1 = HALO_TOOLS::int_to_direction(vp.dir1);
        Direction this_d2 = HALO_TOOLS::int_to_direction(vp.dir2);
        Direction this_d3 = HALO_TOOLS::int_to_direction(vp.dir3);

        int send_dir1 = HALO_TOOLS::map_dir_to_neighbor(vp.dir1, vp.trans, true);
        int send_dir2 = HALO_TOOLS::map_dir_to_neighbor(vp.dir2, vp.trans, false);
        int send_dir3 = HALO_TOOLS::map_dir_to_neighbor(vp.dir3, vp.trans, false);
        if (vp.trans.perm[std::abs(vp.dir1) - 1] + 1 != std::abs(send_dir1))
        {
            if (vp.trans.perm[std::abs(vp.dir1) - 1] + 1 == std::abs(send_dir2))
                std::swap(send_dir1, send_dir2);
            else if (vp.trans.perm[std::abs(vp.dir1) - 1] + 1 == std::abs(send_dir3))
                std::swap(send_dir1, send_dir3);
            else
            {
                std::cout << "Fatal Error!!! triplet parallel vertex direction mapping broken.\n";
                std::exit(-1);
            }
        }

        for (int recv_axis = 0; recv_axis < 3; ++recv_axis)
        {
            TripletCornerMeta m{};
            m.recv_rank = vp.this_rank;
            m.send_rank = vp.nb_rank;
            m.recv_block = vp.this_block;
            m.send_block = vp.nb_block;
            m.recv_axis = recv_axis;
            m.send_axis = vp.trans.perm[recv_axis];
            m.sign = vp.trans.sign[recv_axis];
            m.recv_box = HALO_BOX::make_3DCorner_ghost_box(
                edge_loc_from_axis(recv_axis), vp.this_box_node, this_d1, this_d2, this_d3, nghost);
            m.node_box_on_send = vp.nb_box_node;
            m.dir1_send = send_dir1;
            m.dir2_send = send_dir2;
            m.dir3_send = send_dir3;
            m.trans_recv_to_send = vp.trans;
            m.tag = next_tag[vp.nb_rank]++;

            recv_plan.push_back(m);
            meta_to_send[vp.nb_rank].push_back(m);
        }
    }

    std::vector<TripletCornerMeta> send_plan;
    exchange_triplet_corner_meta(meta_to_send, send_plan);

    const int nsend = static_cast<int>(send_plan.size());
    const int nrecv = static_cast<int>(recv_plan.size());

    std::vector<std::vector<double>> sbuf(nsend);
    std::vector<std::vector<double>> rbuf(nrecv);
    std::vector<MPI_Request> sreq(nsend), rreq(nrecv);
    std::vector<MPI_Status> sstat(nsend), rstat(nrecv);
    std::vector<int> slen(nsend, 0), rlen(nrecv, 0);

    for (int i = 0; i < nsend; ++i)
    {
        const auto &m = send_plan[i];
        if (m.send_rank != myid)
        {
            std::cout << "Fatal Error!!! triplet vertex meta delivered to wrong rank.\n";
            std::exit(-1);
        }

        Direction d1 = HALO_TOOLS::int_to_direction(m.dir1_send);
        Direction d2 = HALO_TOOLS::int_to_direction(m.dir2_send);
        Direction d3 = HALO_TOOLS::int_to_direction(m.dir3_send);
        const Box3 send_box = HALO_BOX::make_3DCorner_innerghost_box(
            edge_loc_from_axis(m.send_axis), m.node_box_on_send, d1, d2, d3, nghost);

        slen[i] = box_volume(m.recv_box) * ncomp;
        sbuf[i].resize(slen[i]);
        pack_triplet_corner_send(fld_, fid, m, send_box, ncomp, sbuf[i]);
    }

    for (int i = 0; i < nrecv; ++i)
    {
        rlen[i] = box_volume(recv_plan[i].recv_box) * ncomp;
        rbuf[i].resize(rlen[i]);
    }

    for (int i = 0; i < nrecv; ++i)
        PARALLEL::mpi_data_recv(recv_plan[i].send_rank, recv_plan[i].tag,
                                rbuf[i].data(), rlen[i], &rreq[i]);
    for (int i = 0; i < nsend; ++i)
        PARALLEL::mpi_data_send(send_plan[i].recv_rank, send_plan[i].tag,
                                sbuf[i].data(), slen[i], &sreq[i]);

    int wait_recv = nrecv;
    int wait_send = nsend;
    if (wait_recv > 0)
        PARALLEL::mpi_wait(wait_recv, rreq.data(), rstat.data());
    if (wait_send > 0)
        PARALLEL::mpi_wait(wait_send, sreq.data(), sstat.data());
    PARALLEL::mpi_barrier();

    for (int i = 0; i < nrecv; ++i)
        unpack_triplet_corner_recv(fld_, fid, recv_plan[i], ncomp, rbuf[i]);
}
