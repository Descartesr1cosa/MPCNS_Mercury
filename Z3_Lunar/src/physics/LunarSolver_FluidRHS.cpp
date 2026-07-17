#include "LunarSolver.h"

void LunarSolver::Scheme_U_()
{
    auto calc_Jac_radius_GCL = [&](double &out,
                                   const double rho, const double u, const double v, const double w, const double p,
                                   const double k1, const double k2, const double k3)
    {
        // k1,k2,k3 are "Jac * contravariant metric components" (GCL form)
        const double uvw = k1 * u + k2 * v + k3 * w;
        const double kk = (k1 * k1 + k2 * k2 + k3 * k3);
        const double c = std::sqrt(std::max(0.0, gamma_ * p / std::max(rho, 1e-30)));
        const double a = std::abs(uvw) + c * std::sqrt(kk);
        out = a;
    };

    auto calc_Jac_Flux_GCL = [&](double *flux,
                                 const double rho, const double u, const double v, const double w, const double p,
                                 const double k1, const double k2, const double k3)
    {
        const double uvw = k1 * u + k2 * v + k3 * w;
        const double rhoe = p / (gamma_ - 1.0) + 0.5 * rho * (u * u + v * v + w * w);

        flux[0] = rho * uvw;
        flux[1] = rho * uvw * u + k1 * p;
        flux[2] = rho * uvw * v + k2 * p;
        flux[3] = rho * uvw * w + k3 * p;
        flux[4] = uvw * (rhoe + p);
    };

    // One-face Rusanov (piecewise constant; you can upgrade to MUSCL later)
    auto Reconstruction_Rusanov = [&](const double *metric, int direction,
                                      FieldBlock &PV, FieldBlock &U,
                                      int i, int j, int k,
                                      double *out_flux)
    {
        const double k1 = metric[0];
        const double k2 = metric[1];
        const double k3 = metric[2];

        // left/right states at the face: (i-1,i) or (j-1,j) or (k-1,k)
        int iL = i, jL = j, kL = k;
        int iR = i, jR = j, kR = k;

        if (direction == 0)
        {
            iL = i - 1;
            iR = i;
        }
        if (direction == 1)
        {
            jL = j - 1;
            jR = j;
        }
        if (direction == 2)
        {
            kL = k - 1;
            kR = k;
        }

        // Left primitive
        const double rhoL = U(iL, jL, kL, 0);
        const double uL = PV(iL, jL, kL, 0);
        const double vL = PV(iL, jL, kL, 1);
        const double wL = PV(iL, jL, kL, 2);
        const double pL = PV(iL, jL, kL, 3);

        // Right primitive
        const double rhoR = U(iR, jR, kR, 0);
        const double uR = PV(iR, jR, kR, 0);
        const double vR = PV(iR, jR, kR, 1);
        const double wR = PV(iR, jR, kR, 2);
        const double pR = PV(iR, jR, kR, 3);

        // Conservative UL/UR
        double UL[5], UR[5];
        UL[0] = rhoL;
        UL[1] = rhoL * uL;
        UL[2] = rhoL * vL;
        UL[3] = rhoL * wL;
        UL[4] = pL / (gamma_ - 1.0) + 0.5 * rhoL * (uL * uL + vL * vL + wL * wL);

        UR[0] = rhoR;
        UR[1] = rhoR * uR;
        UR[2] = rhoR * vR;
        UR[3] = rhoR * wR;
        UR[4] = pR / (gamma_ - 1.0) + 0.5 * rhoR * (uR * uR + vR * vR + wR * wR);

        // Spectral radius
        double radL = 0.0, radR = 0.0;
        calc_Jac_radius_GCL(radL, rhoL, uL, vL, wL, pL, k1, k2, k3);
        calc_Jac_radius_GCL(radR, rhoR, uR, vR, wR, pR, k1, k2, k3);
        const double rad = std::max(radL, radR);

        // Fluxes
        double FL[5], FR[5];
        calc_Jac_Flux_GCL(FL, rhoL, uL, vL, wL, pL, k1, k2, k3);
        calc_Jac_Flux_GCL(FR, rhoR, uR, vR, wR, pR, k1, k2, k3);

        for (int m = 0; m < 5; ++m)
            out_flux[m] = 0.5 * (FL[m] + FR[m]) - 0.5 * rad * (UR[m] - UL[m]);
    };

    auto do_one_species = [&](int fidU, int fidPV, int fidRHS)
    {
        for (int iblk = 0; iblk < fld_->num_blocks(); ++iblk)
        {
            FieldBlock &U = fld_->field(fidU, iblk);
            FieldBlock &PV = fld_->field(fidPV, iblk);
            FieldBlock &RHS = fld_->field(fidRHS, iblk);

            FieldBlock &Jac = fld_->field(fid_.fid_Jac, iblk);
            FieldBlock &XI = fld_->field(fid_.fid_metric.xi, iblk);   // FaceXi,3
            FieldBlock &ET = fld_->field(fid_.fid_metric.eta, iblk);  // FaceEt,3
            FieldBlock &ZE = fld_->field(fid_.fid_metric.zeta, iblk); // FaceZe,3

            if (!U.is_allocated() || !PV.is_allocated() || !RHS.is_allocated())
                continue;

            // temp flux fields (shared is OK; we overwrite per species)
            FieldBlock &Fxi = fld_->field("F_xi", iblk);
            FieldBlock &Fet = fld_->field("F_eta", iblk);
            FieldBlock &Fze = fld_->field("F_zeta", iblk);

            double metric[3];
            double Flux[5];

            // ---- xi faces ----
            {
                Int3 sub = Fxi.inner_lo();
                Int3 sup = Fxi.inner_hi();
                for (int i = sub.i; i < sup.i; ++i)
                    for (int j = sub.j; j < sup.j; ++j)
                        for (int k = sub.k; k < sup.k; ++k)
                        {
                            metric[0] = XI(i, j, k, 0);
                            metric[1] = XI(i, j, k, 1);
                            metric[2] = XI(i, j, k, 2);
                            Reconstruction_Rusanov(metric, 0, PV, U, i, j, k, Flux);
                            for (int m = 0; m < 5; ++m)
                                Fxi(i, j, k, m) = Flux[m];
                        }
            }

            // ---- eta faces ----
            {
                Int3 sub = Fet.inner_lo();
                Int3 sup = Fet.inner_hi();
                for (int i = sub.i; i < sup.i; ++i)
                    for (int j = sub.j; j < sup.j; ++j)
                        for (int k = sub.k; k < sup.k; ++k)
                        {
                            metric[0] = ET(i, j, k, 0);
                            metric[1] = ET(i, j, k, 1);
                            metric[2] = ET(i, j, k, 2);
                            Reconstruction_Rusanov(metric, 1, PV, U, i, j, k, Flux);
                            for (int m = 0; m < 5; ++m)
                                Fet(i, j, k, m) = Flux[m];
                        }
            }

            // ---- zeta faces ----
            {
                Int3 sub = Fze.inner_lo();
                Int3 sup = Fze.inner_hi();
                for (int i = sub.i; i < sup.i; ++i)
                    for (int j = sub.j; j < sup.j; ++j)
                        for (int k = sub.k; k < sup.k; ++k)
                        {
                            metric[0] = ZE(i, j, k, 0);
                            metric[1] = ZE(i, j, k, 1);
                            metric[2] = ZE(i, j, k, 2);
                            Reconstruction_Rusanov(metric, 2, PV, U, i, j, k, Flux);
                            for (int m = 0; m < 5; ++m)
                                Fze(i, j, k, m) = Flux[m];
                        }
            }

            // ---- divergence to RHS ----
            {
                Int3 sub = RHS.inner_lo();
                Int3 sup = RHS.inner_hi();
                for (int i = sub.i; i < sup.i; ++i)
                    for (int j = sub.j; j < sup.j; ++j)
                        for (int k = sub.k; k < sup.k; ++k)
                        {
                            const double invJ = 1.0 / Jac(i, j, k, 0);
                            for (int m = 0; m < 5; ++m)
                            {
                                RHS(i, j, k, m) -= (Fxi(i + 1, j, k, m) - Fxi(i, j, k, m)) * invJ;
                                RHS(i, j, k, m) -= (Fet(i, j + 1, k, m) - Fet(i, j, k, m)) * invJ;
                                RHS(i, j, k, m) -= (Fze(i, j, k + 1, m) - Fze(i, j, k, m)) * invJ;
                            }
                        }
            }
        }
    };

    // Single H+ conservative system.
    do_one_species(fid_.fid_U_H, fid_.fid_PV_H, fid_.fid_RHS_H);
}

// Electromagnetic momentum/energy coupling for the single H+ fluid.
// This package has no second-species creation, drag, or conservation sources.
void LunarSolver::AddSourceToRHS_Fluid()
{
    // Match the cell-centred electron-pressure gradient used by Mercury.
    auto dd = [](double a, double b, double c,
                 double fp_i, double fm_i,
                 double fp_j, double fm_j,
                 double fp_k, double fm_k) -> double
    {
        return 0.5 * (a * (fp_i - fm_i) +
                      b * (fp_j - fm_j) +
                      c * (fp_k - fm_k));
    };

    auto derv_at = [](FieldBlock &Jac, FieldBlock &Axi,
                      FieldBlock &Aet, FieldBlock &Aze,
                      int i, int j, int k,
                      double &ax, double &ay, double &az,
                      double &bx, double &by, double &bz,
                      double &cx, double &cy, double &cz)
    {
        const double V = std::abs(Jac(i, j, k, 0));
        if (V <= 0.0)
        {
            ax = ay = az = bx = by = bz = cx = cy = cz = 0.0;
            return;
        }

        ax = 0.5 * (Axi(i - 1, j, k, 0) + Axi(i, j, k, 0)) / V;
        ay = 0.5 * (Axi(i - 1, j, k, 1) + Axi(i, j, k, 1)) / V;
        az = 0.5 * (Axi(i - 1, j, k, 2) + Axi(i, j, k, 2)) / V;

        bx = 0.5 * (Aet(i, j - 1, k, 0) + Aet(i, j, k, 0)) / V;
        by = 0.5 * (Aet(i, j - 1, k, 1) + Aet(i, j, k, 1)) / V;
        bz = 0.5 * (Aet(i, j - 1, k, 2) + Aet(i, j, k, 2)) / V;

        cx = 0.5 * (Aze(i, j, k - 1, 0) + Aze(i, j, k, 0)) / V;
        cy = 0.5 * (Aze(i, j, k - 1, 1) + Aze(i, j, k, 1)) / V;
        cz = 0.5 * (Aze(i, j, k - 1, 2) + Aze(i, j, k, 2)) / V;
    };

    for (int ib = 0; ib < fld_->num_blocks(); ++ib)
    {
        FieldBlock &Jac = fld_->field(fid_.fid_Jac, ib);
        FieldBlock &Axi = fld_->field(fid_.fid_metric.xi, ib);
        FieldBlock &Aet = fld_->field(fid_.fid_metric.eta, ib);
        FieldBlock &Aze = fld_->field(fid_.fid_metric.zeta, ib);
        FieldBlock &UH = fld_->field(fid_.fid_U_H, ib);
        FieldBlock &PVH = fld_->field(fid_.fid_PV_H, ib);
        FieldBlock &B = fld_->field(fid_.fid_Bcell, ib);
        FieldBlock &J = fld_->field(fid_.fid_Jcell, ib);
        FieldBlock &RHS = fld_->field(fid_.fid_RHS_H, ib);
        const bool missing_ambipolar_geometry =
            ambipolar_control.enabled &&
            (!Jac.is_allocated() || !Axi.is_allocated() ||
             !Aet.is_allocated() || !Aze.is_allocated());
        if (missing_ambipolar_geometry ||
            !UH.is_allocated() || !PVH.is_allocated() ||
            !B.is_allocated() || !J.is_allocated() || !RHS.is_allocated())
            continue;

        const Int3 lo = RHS.inner_lo();
        const Int3 hi = RHS.inner_hi();
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                {
                    const double jxb_x = J(i,j,k,1)*B(i,j,k,2) - J(i,j,k,2)*B(i,j,k,1);
                    const double jxb_y = J(i,j,k,2)*B(i,j,k,0) - J(i,j,k,0)*B(i,j,k,2);
                    const double jxb_z = J(i,j,k,0)*B(i,j,k,1) - J(i,j,k,1)*B(i,j,k,0);
                    const double weight = Hall_Num_Limiter(UH(i,j,k,0)).wH_mhd;
                    const double sx = momentum_hall_coeff * weight * jxb_x;
                    const double sy = momentum_hall_coeff * weight * jxb_y;
                    const double sz = momentum_hall_coeff * weight * jxb_z;
                    RHS(i,j,k,1) += sx;
                    RHS(i,j,k,2) += sy;
                    RHS(i,j,k,3) += sz;
                    RHS(i,j,k,4) += sx*PVH(i,j,k,0) + sy*PVH(i,j,k,1) + sz*PVH(i,j,k,2);

                    if (ambipolar_control.enabled)
                    {
                        // For Lunar's single H+ fluid, pe is the same p_H used
                        // by AddAmbipolarEdgeEMF_.
                        double ax, ay, az, bx, by, bz, cx, cy, cz;
                        derv_at(Jac, Axi, Aet, Aze, i, j, k,
                                ax, ay, az, bx, by, bz, cx, cy, cz);

                        const double pe_p_i = PVH(i + 1, j, k, 3);
                        const double pe_m_i = PVH(i - 1, j, k, 3);
                        const double pe_p_j = PVH(i, j + 1, k, 3);
                        const double pe_m_j = PVH(i, j - 1, k, 3);
                        const double pe_p_k = PVH(i, j, k + 1, 3);
                        const double pe_m_k = PVH(i, j, k - 1, 3);

                        const double dpex = dd(ax, bx, cx, pe_p_i, pe_m_i,
                                               pe_p_j, pe_m_j, pe_p_k, pe_m_k);
                        const double dpey = dd(ay, by, cy, pe_p_i, pe_m_i,
                                               pe_p_j, pe_m_j, pe_p_k, pe_m_k);
                        const double dpez = dd(az, bz, cz, pe_p_i, pe_m_i,
                                               pe_p_j, pe_m_j, pe_p_k, pe_m_k);

                        RHS(i,j,k,1) -= weight * dpex;
                        RHS(i,j,k,2) -= weight * dpey;
                        RHS(i,j,k,3) -= weight * dpez;
                        RHS(i,j,k,4) -= weight *
                            (dpex*PVH(i,j,k,0) +
                             dpey*PVH(i,j,k,1) +
                             dpez*PVH(i,j,k,2));
                    }
                }
    }
}
