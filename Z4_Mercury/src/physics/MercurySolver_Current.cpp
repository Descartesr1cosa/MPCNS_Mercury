#include "MercurySolver.h"

namespace
{
// Keep the shared current operator consistent with Lunar.  Extraordinary
// edges use the full dual polygon; ordinary seams use a deferred correction.
constexpr double regular_seam_curl_correction = 0.25;
}

void MercurySolver::ReduceEdgeAliasCandidatesToOwners_(const IdTriplet &fid_edge)
{
    if(!topo_) return;
    const int myid=par_->GetInt("myid");
    if(!edge_alias_reduce_cache_ready_)
    {
        edge_alias_reduce_classes_.clear();
        edge_alias_reduce_classes_.reserve(topo_->edges.classes.size());
        for(const auto &cls:topo_->edges.classes)
            if(cls.global_id>=0 && cls.members.size()>1)
                edge_alias_reduce_classes_.push_back(&cls);
        std::sort(edge_alias_reduce_classes_.begin(),edge_alias_reduce_classes_.end(),
                  [](const auto *a,const auto *b)
                  { return a->global_id<b->global_id; });
        edge_alias_reduce_local_sum_.resize(edge_alias_reduce_classes_.size());
        edge_alias_reduce_global_sum_.resize(edge_alias_reduce_classes_.size());

        edge_alias_reduce_local_terms_.clear();
        edge_alias_reduce_owner_writes_.clear();
        for(std::size_t n=0;n<edge_alias_reduce_classes_.size();++n)
        {
            const auto &cls=*edge_alias_reduce_classes_[n];
            for(const auto &m:cls.members)
                if(m.entity.rank==myid)
                    edge_alias_reduce_local_terms_.push_back(
                        {n,m.entity,m.orient_sign});
            if(cls.owner.entity.rank==myid)
                edge_alias_reduce_owner_writes_.push_back(
                    {n,cls.owner.entity,cls.owner.orient_sign,
                     static_cast<double>(cls.members.size())});
        }
        edge_alias_reduce_cache_ready_=true;
    }

    const auto &classes=edge_alias_reduce_classes_;
    const std::size_t nclass=classes.size();
    if(nclass==0) return;
    auto &local_sum=edge_alias_reduce_local_sum_;
    auto &global_sum=edge_alias_reduce_global_sum_;
    std::fill(local_sum.begin(),local_sum.end(),0.0);
    for(const auto &term:edge_alias_reduce_local_terms_)
    {
        const auto &e=term.edge;
        FieldBlock &value=fld_->field(fid_edge.at(static_cast<int>(e.axis)+1),e.block);
        local_sum[term.class_index]+=term.orient_sign*value(e.i,e.j,e.k,0);
    }
    MPI_Allreduce(local_sum.data(),global_sum.data(),static_cast<int>(nclass),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    for(const auto &write:edge_alias_reduce_owner_writes_)
    {
        const auto &e=write.edge;
        FieldBlock &owner=fld_->field(fid_edge.at(static_cast<int>(e.axis)+1),e.block);
        owner(e.i,e.j,e.k,0)=write.orient_sign*
            global_sum[write.class_index]/write.member_count;
    }
}

void MercurySolver::AssembleSingularEdgeCurrent_(const IdTriplet &fid_Bface,
                                                 const IdTriplet &fid_Jedge)
{
    if (!singular_edges_ || singular_edges_->empty()) return;
    (void)fid_Bface;
    singular_edges_->assemble_cell_vector_circulation_to_local_owners(
        *fld_, "Bind_cell",
        {fld_->descriptor(fid_Jedge.xi).name,
         fld_->descriptor(fid_Jedge.eta).name,
         fld_->descriptor(fid_Jedge.zeta).name},
        regular_seam_curl_correction);
}

void MercurySolver::Calc_J_Edge()
{
    for (int iblk = 0; iblk < fld_->num_blocks(); ++iblk)
    {
        // B_xi/B_eta/B_zeta are the Bind face-flux DOFs.  Badd is prescribed
        // analytically and must not contribute to the solved current.
        auto &Bxi = fld_->field(fid_.fid_B.xi, iblk);
        auto &Beta = fld_->field(fid_.fid_B.eta, iblk);
        auto &Bzeta = fld_->field(fid_.fid_B.zeta, iblk);

        auto &Jxi = fld_->field(fid_.fid_J.xi, iblk);
        auto &Jeta = fld_->field(fid_.fid_J.eta, iblk);
        auto &Jzeta = fld_->field(fid_.fid_J.zeta, iblk);

        auto &beta_xi = fld_->field(fid_.Face_beta.xi, iblk); // Hodge *: 2-form face -> 1-form face
        auto &beta_eta = fld_->field(fid_.Face_beta.eta, iblk);
        auto &beta_zeta = fld_->field(fid_.Face_beta.zeta, iblk);

        auto &alpha_xi = fld_->field(fid_.Edge_alpha.xi, iblk); // Hodge *: 2-form edge -> 1-form edge
        auto &alpha_eta = fld_->field(fid_.Edge_alpha.eta, iblk);
        auto &alpha_zeta = fld_->field(fid_.Edge_alpha.zeta, iblk);

        // compute J (edge 1-form) from face B (2-form)
        // multiper 用 +1.0 J =curl B。
        CTOperators::CurlAdjFaceToEdge(iblk,
                                       Bxi, Beta, Bzeta,
                                       beta_xi, beta_eta, beta_zeta,
                                       Jxi, Jeta, Jzeta,
                                       /*multiper=*/1.0);

        //  J_edge^(1-form) = (alpha_edge / mu0) * Jcirc_edge
        //     alpha = |e|/|S*|  (⋆1^{-1})
        {
            // Edge_xi
            Int3 lo = Jxi.inner_lo(), hi = Jxi.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        Jxi(i, j, k, 0) *= (alpha_xi(i, j, k, 0));
        }
        {
            // Edge_eta
            Int3 lo = Jeta.inner_lo(), hi = Jeta.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        Jeta(i, j, k, 0) *= (alpha_eta(i, j, k, 0));
        }
        {
            // Edge_zeta
            Int3 lo = Jzeta.inner_lo(), hi = Jzeta.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        Jzeta(i, j, k, 0) *= (alpha_zeta(i, j, k, 0));
        }
    }

    ReduceEdgeAliasCandidatesToOwners_(fid_.fid_J);
    AssembleSingularEdgeCurrent_(fid_.fid_B, fid_.fid_J);
    mercury_bound_.Sync("Jedge"); // ApplyBC_EdgeJ_();
}

void MercurySolver::calc_Jcell()
{
    const int nblock = fld_->num_blocks();

    for (int ib = 0; ib < nblock; ++ib)
    {
        auto &Jcell = fld_->field(fid_.fid_Jcell, ib);

        auto &Jxi = fld_->field(fid_.fid_J.xi, ib);
        auto &Jeta = fld_->field(fid_.fid_J.eta, ib);
        auto &Jzeta = fld_->field(fid_.fid_J.zeta, ib);

        auto &W = fld_->field(fid_.fid_Jcell_from_Jedge_w, ib);

        if (!Jcell.is_allocated() || !Jxi.is_allocated() ||
            !Jeta.is_allocated() || !Jzeta.is_allocated())
            continue;

        Int3 lo = Jcell.inner_lo();
        Int3 hi = Jcell.inner_hi();

        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                {
                    const double s0 = Jxi(i, j, k, 0);
                    const double s1 = Jxi(i, j + 1, k, 0);
                    const double s2 = Jxi(i, j, k + 1, 0);
                    const double s3 = Jxi(i, j + 1, k + 1, 0);

                    const double s4 = Jeta(i, j, k, 0);
                    const double s5 = Jeta(i + 1, j, k, 0);
                    const double s6 = Jeta(i, j, k + 1, 0);
                    const double s7 = Jeta(i + 1, j, k + 1, 0);

                    const double s8 = Jzeta(i, j, k, 0);
                    const double s9 = Jzeta(i + 1, j, k, 0);
                    const double s10 = Jzeta(i, j + 1, k, 0);
                    const double s11 = Jzeta(i + 1, j + 1, k, 0);

                    const double s[12] = {s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11};

                    double Jx = 0.0, Jy = 0.0, Jz = 0.0;
                    for (int n = 0; n < 12; ++n)
                    {
                        const double sn = s[n];
                        Jx += W(i, j, k, n) * sn;
                        Jy += W(i, j, k, 12 + n) * sn;
                        Jz += W(i, j, k, 24 + n) * sn;
                    }

                    Jcell(i, j, k, 0) = Jx;
                    Jcell(i, j, k, 1) = Jy;
                    Jcell(i, j, k, 2) = Jz;
                }
    }

    if (topo_)
    {
        auto solve3x3 = [](double A[3][3], double b[3], double x[3]) -> bool
        {
            const double trA = A[0][0] + A[1][1] + A[2][2];
            const double reg = 1.0e-12 * std::max(1.0, trA);

            A[0][0] += reg;
            A[1][1] += reg;
            A[2][2] += reg;

            double M[3][4] = {
                {A[0][0], A[0][1], A[0][2], b[0]},
                {A[1][0], A[1][1], A[1][2], b[1]},
                {A[2][0], A[2][1], A[2][2], b[2]}};

            for (int c = 0; c < 3; ++c)
            {
                int piv = c;
                double amax = std::abs(M[c][c]);

                for (int r0 = c + 1; r0 < 3; ++r0)
                {
                    const double av = std::abs(M[r0][c]);
                    if (av > amax)
                    {
                        amax = av;
                        piv = r0;
                    }
                }

                if (amax < 1.0e-300)
                    return false;

                if (piv != c)
                {
                    for (int q = c; q < 4; ++q)
                        std::swap(M[c][q], M[piv][q]);
                }

                const double inv = 1.0 / M[c][c];
                for (int q = c; q < 4; ++q)
                    M[c][q] *= inv;

                for (int r0 = 0; r0 < 3; ++r0)
                {
                    if (r0 == c)
                        continue;

                    const double fac = M[r0][c];
                    for (int q = c; q < 4; ++q)
                        M[r0][q] -= fac * M[c][q];
                }
            }

            x[0] = M[0][3];
            x[1] = M[1][3];
            x[2] = M[2][3];
            return true;
        };

        constexpr double eps_len = 1.0e-14;

        for (const auto &p : topo_->physical_patches)
        {
            if (p.bc_name != "Pole")
                continue;

            const int dir = std::abs(p.direction);
            if (dir != 1 && dir != 2)
                continue;

            const int ib = p.this_block;
            if (ib < 0 || ib >= nblock)
                continue;

            auto &Jcell = fld_->field(fid_.fid_Jcell, ib);
            auto &Jxi = fld_->field(fid_.fid_J.xi, ib);
            auto &Jeta = fld_->field(fid_.fid_J.eta, ib);
            auto &Jzeta = fld_->field(fid_.fid_J.zeta, ib);
            auto &dr_xi = fld_->field(fid_.Edge_dr.xi, ib);
            auto &dr_eta = fld_->field(fid_.Edge_dr.eta, ib);
            auto &dr_zeta = fld_->field(fid_.Edge_dr.zeta, ib);

            if (!Jcell.is_allocated() || !Jxi.is_allocated() ||
                !Jeta.is_allocated() || !Jzeta.is_allocated() ||
                !dr_xi.is_allocated() || !dr_eta.is_allocated() ||
                !dr_zeta.is_allocated())
                continue;

            const bool high_side = (p.direction > 0);
            const Box3 &node = p.this_box_node;

            auto push_edge = [&](FieldBlock &Je, FieldBlock &dr,
                                 int i, int j, int k,
                                 double A[3][3], double b[3], int &cnt)
            {
                const double dx = dr(i, j, k, 0);
                const double dy = dr(i, j, k, 1);
                const double dz = dr(i, j, k, 2);
                const double L = std::sqrt(dx * dx + dy * dy + dz * dz);

                if (L < eps_len)
                    return;

                const double tau[3] = {dx / L, dy / L, dz / L};
                const double y = Je(i, j, k, 0) / L;

                for (int a = 0; a < 3; ++a)
                {
                    b[a] += y * tau[a];
                    for (int c = 0; c < 3; ++c)
                        A[a][c] += tau[a] * tau[c];
                }

                ++cnt;
            };

            auto write_ring = [&](int ic, int jc,
                                  double A[3][3], double b[3], int cnt,
                                  int klo, int khi)
            {
                if (cnt <= 0)
                    return;

                double Jp[3] = {0.0, 0.0, 0.0};
                if (!solve3x3(A, b, Jp))
                    return;

                for (int k = klo; k < khi; ++k)
                {
                    Jcell(ic, jc, k, 0) = Jp[0];
                    Jcell(ic, jc, k, 1) = Jp[1];
                    Jcell(ic, jc, k, 2) = Jp[2];
                }
            };

            if (dir == 1)
            {
                const int ic = high_side ? (node.lo.i - 1) : node.lo.i;
                const int it = high_side ? ic : (ic + 1);

                for (int j = node.lo.j; j < node.hi.j - 1; ++j)
                {
                    double A[3][3] = {
                        {0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0}};
                    double b[3] = {0.0, 0.0, 0.0};
                    int cnt = 0;

                    for (int k = node.lo.k; k < node.hi.k - 1; ++k)
                    {
                        // norm edges on the Pole layer
                        push_edge(Jxi, dr_xi, ic, j, k, A, b, cnt);
                        push_edge(Jxi, dr_xi, ic, j + 1, k, A, b, cnt);
                        push_edge(Jxi, dr_xi, ic, j, k + 1, A, b, cnt);
                        push_edge(Jxi, dr_xi, ic, j + 1, k + 1, A, b, cnt);

                        // axis/rotate edges one layer away from the Pole
                        push_edge(Jeta, dr_eta, it, j, k, A, b, cnt);
                        push_edge(Jeta, dr_eta, it, j, k + 1, A, b, cnt);
                        push_edge(Jzeta, dr_zeta, it, j, k, A, b, cnt);
                        push_edge(Jzeta, dr_zeta, it, j + 1, k, A, b, cnt);
                    }

                    write_ring(ic, j, A, b, cnt, node.lo.k, node.hi.k - 1);
                }
            }
            else
            {
                const int jc = high_side ? (node.lo.j - 1) : node.lo.j;
                const int jt = high_side ? jc : (jc + 1);

                for (int i = node.lo.i; i < node.hi.i - 1; ++i)
                {
                    double A[3][3] = {
                        {0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0}};
                    double b[3] = {0.0, 0.0, 0.0};
                    int cnt = 0;

                    for (int k = node.lo.k; k < node.hi.k - 1; ++k)
                    {
                        // norm edges on the Pole layer
                        push_edge(Jeta, dr_eta, i, jc, k, A, b, cnt);
                        push_edge(Jeta, dr_eta, i + 1, jc, k, A, b, cnt);
                        push_edge(Jeta, dr_eta, i, jc, k + 1, A, b, cnt);
                        push_edge(Jeta, dr_eta, i + 1, jc, k + 1, A, b, cnt);

                        // axis/rotate edges one layer away from the Pole
                        push_edge(Jxi, dr_xi, i, jt, k, A, b, cnt);
                        push_edge(Jxi, dr_xi, i, jt, k + 1, A, b, cnt);
                        push_edge(Jzeta, dr_zeta, i, jt, k, A, b, cnt);
                        push_edge(Jzeta, dr_zeta, i + 1, jt, k, A, b, cnt);
                    }

                    write_ring(i, jc, A, b, cnt, node.lo.k, node.hi.k - 1);
                }
            }
        }
    }

    mercury_bound_.Sync("J_cell");

}
