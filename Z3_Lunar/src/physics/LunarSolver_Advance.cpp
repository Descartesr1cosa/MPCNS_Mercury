
#include "LunarSolver.h"

#include <algorithm>
#include <cmath>

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

void LunarSolver::ApplyDensityFloor_()
{
    if (!(density_floor_ > 0.0) || !std::isfinite(density_floor_))
        return;

    constexpr double p_floor = 1.0e-12;
    const double gm1 = std::max(gamma_ - 1.0, 1.0e-12);

    for (int ib = 0; ib < fld_->num_blocks(); ++ib)
    {
        FieldBlock &U = fld_->field(fid_.fid_U_H, ib);
        if (!U.is_allocated())
            continue;

        const Int3 lo = U.inner_lo();
        const Int3 hi = U.inner_hi();
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                {
                    const double rho_old = U(i, j, k, 0);
                    if (rho_old >= density_floor_)
                        continue;

                    const double mx = U(i, j, k, 1);
                    const double my = U(i, j, k, 2);
                    const double mz = U(i, j, k, 3);
                    const double energy_old = U(i, j, k, 4);

                    if (rho_old > 0.0 && std::isfinite(rho_old) &&
                        std::isfinite(mx) && std::isfinite(my) &&
                        std::isfinite(mz) && std::isfinite(energy_old))
                    {
                        // Preserve the local velocity and thermal pressure.
                        // Added floor mass therefore carries the same velocity;
                        // momentum and total energy must both be reconstructed.
                        const double ux = mx / rho_old;
                        const double uy = my / rho_old;
                        const double uz = mz / rho_old;
                        const double speed2 = ux * ux + uy * uy + uz * uz;
                        const double kinetic_old = 0.5 * rho_old * speed2;
                        const double internal_raw = energy_old - kinetic_old;

                        if (!std::isfinite(speed2) ||
                            !std::isfinite(internal_raw))
                        {
                            U(i, j, k, 0) = density_floor_;
                            U(i, j, k, 1) = 0.0;
                            U(i, j, k, 2) = 0.0;
                            U(i, j, k, 3) = 0.0;
                            U(i, j, k, 4) = p_floor / gm1;
                            continue;
                        }

                        const double internal = std::max(
                            internal_raw, p_floor / gm1);

                        U(i, j, k, 0) = density_floor_;
                        U(i, j, k, 1) = density_floor_ * ux;
                        U(i, j, k, 2) = density_floor_ * uy;
                        U(i, j, k, 3) = density_floor_ * uz;
                        U(i, j, k, 4) = internal +
                            0.5 * density_floor_ * speed2;
                    }
                    else
                    {
                        // A non-positive/non-finite density has no meaningful
                        // primitive velocity. Replace it with the least-energetic
                        // admissible floor state so the next CFL scan can fail
                        // neither through division by zero nor through NaN.
                        U(i, j, k, 0) = density_floor_;
                        U(i, j, k, 1) = 0.0;
                        U(i, j, k, 2) = 0.0;
                        U(i, j, k, 3) = 0.0;
                        U(i, j, k, 4) = p_floor / gm1;
                    }
                }
    }
}
