#include "4_halo/1_MPCNS_Halo.h"

#include <array>
#include <cstdlib>
#include <iostream>

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

        // Parallel corner patterns are currently keyed by one destination
        // field. Keep the existing transport there until the corner metadata
        // also carries source/destination axes.
        for (const auto &fn : fields)
        {
            std::string name = fn;
            exchange_parallel_edge(name);
        }
        return;
    }

    exchange_inner_vertex_edge_1form_triplet_(fields);
    for (const auto &fn : fields)
    {
        std::string name = fn;
        exchange_parallel_vertex(name);
    }
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
