
#include "LunarSolver.h"

void LunarSolver::ZeroRHS_()
{
    const int nb = fld_->num_blocks();
    for (int ib = 0; ib < nb; ++ib)
    {
        FieldBlock &Jac = fld_->field(fid_.fid_Jac, ib);
        FieldBlock &RHSH = fld_->field(fid_.fid_RHS_H, ib);
        FieldBlock &RHSB_xi = fld_->field(fid_.fid_RHS_b.xi, ib);
        FieldBlock &RHSB_et = fld_->field(fid_.fid_RHS_b.eta, ib);
        FieldBlock &RHSB_ze = fld_->field(fid_.fid_RHS_b.zeta, ib);

        Int3 lo; // = Jac.inner_lo();
        Int3 hi; // = Jac.inner_hi();

        if (RHSH.is_allocated())
        {
            lo = RHSH.inner_lo();
            hi = RHSH.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        for (int m = 0; m < 5; ++m)
                            RHSH(i, j, k, m) = 0.0;
        }

        if (RHSB_xi.is_allocated())
        {
            lo = RHSB_xi.inner_lo();
            hi = RHSB_xi.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        RHSB_xi(i, j, k, 0) = 0.0;
        }

        if (RHSB_et.is_allocated())
        {
            lo = RHSB_et.inner_lo();
            hi = RHSB_et.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        RHSB_et(i, j, k, 0) = 0.0;
        }

        if (RHSB_ze.is_allocated())
        {
            lo = RHSB_ze.inner_lo();
            hi = RHSB_ze.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        RHSB_ze(i, j, k, 0) = 0.0;
        }
    }
}

void LunarSolver::ApplyUpdate_Euler_()
{
    const int nb = fld_->num_blocks();
    for (int ib = 0; ib < nb; ++ib)
    {
        FieldBlock &Jac = fld_->field(fid_.fid_Jac, ib);

        FieldBlock &UH = fld_->field(fid_.fid_U_H, ib);
        FieldBlock &Ub_xi = fld_->field(fid_.fid_B.xi, ib);
        FieldBlock &Ub_eta = fld_->field(fid_.fid_B.eta, ib);
        FieldBlock &Ub_zeta = fld_->field(fid_.fid_B.zeta, ib);

        FieldBlock &RHSH = fld_->field(fid_.fid_RHS_H, ib);
        FieldBlock &RHSB_xi = fld_->field(fid_.fid_RHS_b.xi, ib);
        FieldBlock &RHSB_eta = fld_->field(fid_.fid_RHS_b.eta, ib);
        FieldBlock &RHSB_zeta = fld_->field(fid_.fid_RHS_b.zeta, ib);

        Int3 lo; //= Jac.inner_lo();
        Int3 hi; //= Jac.inner_hi();

        if (RHSH.is_allocated())
        {
            lo = UH.inner_lo();
            hi = UH.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        for (int m = 0; m < 5; ++m)
                            UH(i, j, k, m) += dt * RHSH(i, j, k, m);
        }

        if (RHSB_xi.is_allocated())
        {
            lo = Ub_xi.inner_lo();
            hi = Ub_xi.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        Ub_xi(i, j, k, 0) += dt * RHSB_xi(i, j, k, 0);
        }

        if (RHSB_eta.is_allocated())
        {
            lo = Ub_eta.inner_lo();
            hi = Ub_eta.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        Ub_eta(i, j, k, 0) += dt * RHSB_eta(i, j, k, 0);
        }

        if (RHSB_zeta.is_allocated())
        {
            lo = Ub_zeta.inner_lo();
            hi = Ub_zeta.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        Ub_zeta(i, j, k, 0) += dt * RHSB_zeta(i, j, k, 0);
        }
    }
}
