#include "4_Hall_Implicit.h"
#if HALL_IMPLICIT == 1

void ImplicitHallSolver::ClearEdgeTriplet_(const IdTriplet &fid_triplet)
{
    for (int ib = 0; ib < fld_->num_blocks(); ++ib)
    {
        for (int fid : {fid_triplet.xi, fid_triplet.eta, fid_triplet.zeta})
        {
            auto &F = fld_->field(fid, ib);
            if (!F.is_allocated())
                continue;
            Int3 lo = F.get_lo(), hi = F.get_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        for (int m = 0; m < F.descriptor().ncomp; ++m)
                            F(i, j, k, m) = 0.0;
        }
    }
}

void ImplicitHallSolver::ClearFaceTriplet_(const IdTriplet &fid_triplet)
{
    for (int ib = 0; ib < fld_->num_blocks(); ++ib)
    {
        for (int fid : {fid_triplet.xi, fid_triplet.eta, fid_triplet.zeta})
        {
            auto &F = fld_->field(fid, ib);
            if (!F.is_allocated())
                continue;
            Int3 lo = F.inner_lo(), hi = F.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        for (int m = 0; m < F.descriptor().ncomp; ++m)
                            F(i, j, k, m) = 0.0;
        }
    }
}

void ImplicitHallSolver::CopyEhallToE_()
{
    const int nb = fld_->num_blocks();
    for (int ib = 0; ib < nb; ++ib)
    {
        auto &EHx = fld_->field(fid_.fid_Ehall.xi, ib);
        auto &EHy = fld_->field(fid_.fid_Ehall.eta, ib);
        auto &EHz = fld_->field(fid_.fid_Ehall.zeta, ib);

        auto &Ex = fld_->field(fid_.fid_E.xi, ib);
        auto &Ey = fld_->field(fid_.fid_E.eta, ib);
        auto &Ez = fld_->field(fid_.fid_E.zeta, ib);

        if (!Ex.is_allocated())
            continue;

        auto copy0 = [](FieldBlock &src, FieldBlock &dst)
        {
            Int3 lo = dst.get_lo(), hi = dst.get_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        dst(i, j, k, 0) = src(i, j, k, 0);
        };

        copy0(EHx, Ex);
        copy0(EHy, Ey);
        copy0(EHz, Ez);
    }
}

void ImplicitHallSolver::PackFaceInner_(int fid, std::vector<Scalar> &buf)
{
    for (int ib = 0; ib < fld_->num_blocks(); ++ib)
    {
        auto &F = fld_->field(fid, ib);
        if (!F.is_allocated())
            continue;

        Int3 lo = F.inner_lo(), hi = F.inner_hi();

        auto &v = buf[ib];

        size_t t = 0;
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                    v(i, j, k) = F(i, j, k, 0);
    }
}

void ImplicitHallSolver::RestoreFaceInner_(int fid, std::vector<Scalar> &buf)
{
    for (int ib = 0; ib < fld_->num_blocks(); ++ib)
    {
        auto &F = fld_->field(fid, ib);
        if (!F.is_allocated())
            continue;

        Int3 lo = F.inner_lo(), hi = F.inner_hi();
        auto &v = buf[static_cast<size_t>(ib)];

        size_t t = 0;
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                    F(i, j, k, 0) = v(i, j, k);
    }
}

void ImplicitHallSolver::SnapshotCurrentBface_()
{
    PackFaceInner_(fid_.fid_B.xi, Bstar_xi_);
    PackFaceInner_(fid_.fid_B.eta, Bstar_eta_);
    PackFaceInner_(fid_.fid_B.zeta, Bstar_ze_);
}

void ImplicitHallSolver::RestoreCurrentBfaceFromSnapshot_()
{
    RestoreFaceInner_(fid_.fid_B.xi, Bstar_xi_);
    RestoreFaceInner_(fid_.fid_B.eta, Bstar_eta_);
    RestoreFaceInner_(fid_.fid_B.zeta, Bstar_ze_);
}

void ImplicitHallSolver::AddFaceInnerFromRHS_(int fid_B, int fid_RHS, double factor)
{
    for (int ib = 0; ib < fld_->num_blocks(); ++ib)
    {
        auto &B = fld_->field(fid_B, ib);
        auto &R = fld_->field(fid_RHS, ib);
        if (!B.is_allocated())
            continue;

        Int3 lo = B.inner_lo(), hi = B.inner_hi();
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                    B(i, j, k, 0) += factor * R(i, j, k, 0);
    }
}

#endif
