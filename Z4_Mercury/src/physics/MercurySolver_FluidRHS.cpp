#include "MercurySolver.h"

void MercurySolver::Scheme_U_()
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

    // H+ and Na+ (no split, still one function)
    do_one_species(fid_.fid_U_H, fid_.fid_PV_H, fid_.fid_RHS_H);
    do_one_species(fid_.fid_U_Na, fid_.fid_PV_Na, fid_.fid_RHS_Na);
}

// 添加 Fortran source_species 中对流体 dq 的“生效项”
// 依赖字段：U_H/U_Na, PV_H/PV_Na, U_plus, B_cell, J_cell, Na(neutral), Photo_rate, Jac, metric(Axi/Aet/Aze)
void MercurySolver::AddSourceToRHS_Fluid()
{
    // ---------- constants  ----------
    const double Tn0 = 185.0;  // K
    const double sk1 = 5.0e-5; // 1/s  (day side)
    const double sk2 = 1.0e-5; // 1/s  (night side)

    // coefficients matching Fortran structure (see explanation in my previous message)

    const double a1_Na = (L_ref / U_ref) * (m_Na / rho_ref) * 1e6; // mass source coeff for Na+

    // a4 = L_ref / U_ref(秒)，Fortran : a4 = sl8 / v8
    const double a4 = (L_ref / U_ref);

    const double a5 = (3.0 * L_ref * k_Boltz * n_ref) / (rho_ref * U_ref * U_ref * U_ref);

    const double a6 = (1.0e6 * L_ref * k_Boltz) / (rho_ref * U_ref * U_ref * U_ref * (gamma_ - 1.0));

    // dd: Fortran 0.5*(a*(f(i+1)-f(i-1)) + b*(f(j+1)-f(j-1)) + c*(f(k+1)-f(k-1)))
    auto dd = [](double a, double b, double c,
                 double fp_i, double fm_i,
                 double fp_j, double fm_j,
                 double fp_k, double fm_k) -> double
    {
        return 0.5 * (a * (fp_i - fm_i) + b * (fp_j - fm_j) + c * (fp_k - fm_k));
    };

    // derv_at: 从 (Axi/Aet/Aze)/Jac 还原 derv(1..3,1..3) 的 ax..cz（和 1_Solver_Scheme_B.cpp 一致）
    auto derv_at = [&](FieldBlock &Jac, FieldBlock &Axi, FieldBlock &Aet, FieldBlock &Aze,
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

    const int nb = fld_->num_blocks();
    for (int ib = 0; ib < nb; ++ib)
    {
        FieldBlock &Jac = fld_->field(fid_.fid_Jac, ib);
        FieldBlock &Axi = fld_->field(fid_.fid_metric.xi, ib);
        FieldBlock &Aet = fld_->field(fid_.fid_metric.eta, ib);
        FieldBlock &Aze = fld_->field(fid_.fid_metric.zeta, ib);

        FieldBlock &UH = fld_->field(fid_.fid_U_H, ib);
        FieldBlock &UNa = fld_->field(fid_.fid_U_Na, ib);
        FieldBlock &PVH = fld_->field(fid_.fid_PV_H, ib);
        FieldBlock &PVN = fld_->field(fid_.fid_PV_Na, ib);

        FieldBlock &Up = fld_->field(fid_.fid_U_plus, ib);
        FieldBlock &Bt = fld_->field(fid_.fid_Bcell, ib);
        auto &Jc = fld_->field(fid_.fid_Jcell, ib);
        FieldBlock &NaNeu = fld_->field(fld_->field_id("Na"), ib);
        FieldBlock &Photo = fld_->field(fld_->field_id("Photo_rate"), ib);

        FieldBlock &RHS_H = fld_->field(fid_.fid_RHS_H, ib);
        FieldBlock &RHS_Na = fld_->field(fid_.fid_RHS_Na, ib);

        auto &cx = grd_->grids(ib).dual_x;
        auto &cy = grd_->grids(ib).dual_y;
        auto &cz = grd_->grids(ib).dual_z;

        auto radius = [&](int i, int j, int k) -> double
        {
            double xx = cx(i + 1, j + 1, k + 1);
            double yy = cy(i + 1, j + 1, k + 1);
            double zz = cz(i + 1, j + 1, k + 1);
            return std::sqrt(xx * xx + yy * yy + zz * zz);
        };

        if (!Jac.is_allocated() || !Axi.is_allocated())
            continue;
        if (!UH.is_allocated() || !UNa.is_allocated())
            continue;
        if (!PVH.is_allocated() || !PVN.is_allocated())
            continue;
        if (!Up.is_allocated() || !Bt.is_allocated())
            continue;
        if (!RHS_H.is_allocated() || !RHS_Na.is_allocated())
            continue;
        if (!NaNeu.is_allocated() || !Photo.is_allocated())
            continue;

        Int3 lo = Jac.inner_lo();
        Int3 hi = Jac.inner_hi();

        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                {
                    // ---------- common state ----------
                    const double upx = Up(i, j, k, 0);
                    const double upy = Up(i, j, k, 1);
                    const double upz = Up(i, j, k, 2);

                    const double Bx = Bt(i, j, k, 0);
                    const double By = Bt(i, j, k, 1);
                    const double Bz = Bt(i, j, k, 2);

                    const double Jx = Jc(i, j, k, 0);
                    const double Jy = Jc(i, j, k, 1);
                    const double Jz = Jc(i, j, k, 2);

                    const double sjbx = Jy * Bz - Jz * By;
                    const double sjby = Jz * Bx - Jx * Bz;
                    const double sjbz = Jx * By - Jy * Bx;

                    // double num[3];
                    // Hall_Num_Limiter(UH(i, j, k, 0), UNa(i, j, k, 0), num);
                    // const double nH = num[0];
                    // const double nNa = num[1];
                    // const double ne = num[2];

                    NumInfo num = Hall_Num_Limiter(UH(i, j, k, 0), UNa(i, j, k, 0));
                    const double nH = num.nH_true;
                    const double nNa = num.nNa_true;
                    const double chi_H = num.wH_mhd;   // chiH;
                    const double chi_Na = num.wNa_mhd; // chi_Na;

                    // // ---------- Hall coefficient (same convention as induction) ----------
                    // const double nH_hall = UH(i, j, k, 0) / M_H;
                    // const double nNa_hall = UNa(i, j, k, 0) / M_Na;
                    // const double ne_hall = nH_hall + nNa_hall;

                    // const double rloc = radius(i, j, k);
                    // const double alphaH = HallAlpha_Coeffient(ne_hall, 1.5);
                    // // const double alphaH = HallAlpha_Coeffient(ne_hall, rloc);

                    // // ---------- physical number densities for induce term ----------
                    // const double rhoH_nd = std::max(UH(i, j, k, 0), 0.0);
                    // const double rhoNa_nd = std::max(UNa(i, j, k, 0), 0.0);

                    // const double nH_m = (rho_ref * rhoH_nd) / m_H;
                    // const double nNa_m = (rho_ref * rhoNa_nd) / m_Na;

                    // -------------------- metric derv(ax..cz) --------------------
                    // double ax, ay, az, bx, by, bz, cx, cy, cz;
                    // derv_at(Jac, Axi, Aet, Aze, i, j, k, ax, ay, az, bx, by, bz, cx, cy, cz);
                    // -------------------- grad(pe) --------------------
                    // double dpex, dpey, dpez;
                    // {
                    //     const double pe_p_i = PVH(i + 1, j, k, 3) + PVN(i + 1, j, k, 3);
                    //     const double pe_m_i = PVH(i - 1, j, k, 3) + PVN(i - 1, j, k, 3);

                    //     double pe_p_j = 0, pe_m_j = 0, pe_p_k = 0, pe_m_k = 0;

                    //     pe_p_j = PVH(i, j + 1, k, 3) + PVN(i, j + 1, k, 3);
                    //     pe_m_j = PVH(i, j - 1, k, 3) + PVN(i, j - 1, k, 3);

                    //     pe_p_k = PVH(i, j, k + 1, 3) + PVN(i, j, k + 1, 3);
                    //     pe_m_k = PVH(i, j, k - 1, 3) + PVN(i, j, k - 1, 3);

                    //     // dpex = dd(ax, bx, cx, pe_p_i, pe_m_i, pe_p_j, pe_m_j, pe_p_k, pe_m_k);
                    //     // dpey = dd(ay, by, cy, pe_p_i, pe_m_i, pe_p_j, pe_m_j, pe_p_k, pe_m_k);
                    //     // dpez = dd(az, bz, cz, pe_p_i, pe_m_i, pe_p_j, pe_m_j, pe_p_k, pe_m_k);

                    //     dpex = 0.0;
                    //     dpey = 0.0;
                    //     dpez = 0.0;
                    // }

                    // ---------- ionization source sss (cm^-3 s^-1), only for Na+ ----------
                    const double sss = Photo(i, j, k, 0); // cm^-3 s^-1

                    // -------------------- drag frequency vst (1/s) --------------------
                    // Fortran: if(x>=0) vst=sk1 else vst=sk2  （x 为无量纲坐标）
                    // C++ 目前没有直接取物理坐标，这里用 Jac/metric 无法得到 x；
                    // 因此：如果你要严格复现 Fortran 的 x 判定，请从几何场里取坐标场（例如 fid_.fid_coord）再判断。
                    // 这里给一个保守默认：统一用 sk1（你也可以改成 sk2 或按你已有的坐标场实现）
                    // const double sk1 = 5E-5;
                    // const double sk2 = 1E-5;
                    const double vst = (cx(i, j, k) >= 0) ? sk1 : sk2;

                    // b2 = (sm2/(sm1+sm2))*vst ;  b1 = (Tn0 - Ts0)*sm1/(sm1+sm2)

                    // sse = qm1 (Fortran), here Photo is already (cm^-3 s^-1) For electrics
                    const double sse = Photo(i, j, k, 0);

                    // =====================
                    // species H+  (ls=1)
                    // =====================
                    {
                        const double rho = std::max(UH(i, j, k, 0), 0.0);
                        const double u = PVH(i, j, k, 0);
                        const double v = PVH(i, j, k, 1);
                        const double w = PVH(i, j, k, 2);
                        const double u2 = u * u + v * v + w * w;

                        const double subx = (v - upy) * Bz - (w - upz) * By;
                        const double suby = (w - upz) * Bx - (u - upx) * Bz;
                        const double subz = (u - upx) * By - (v - upy) * Bx;

                        const double subu = subx * u + suby * v + subz * w;
                        const double sjbu = sjbx * u + sjby * v + sjbz * w;
                        // const double dpeu = dpex * u + dpey * v + dpez * w;

                        // Ts0 in Kelvin
                        const double Ts0 = PVH(i, j, k, 4) * T_ref;
                        const double b1 = vst * (Tn0 - Ts0) * (M_H / (M_H + M_Na));
                        const double b2 = (M_Na / (M_H + M_Na)) * vst;

                        //=============================================================================================
                        // Electromagnetic Source Terms
                        RHS_H(i, j, k, 1) += momentum_induce_coeff * nH * subx + momentum_hall_coeff * chi_H * sjbx;
                        RHS_H(i, j, k, 2) += momentum_induce_coeff * nH * suby + momentum_hall_coeff * chi_H * sjby;
                        RHS_H(i, j, k, 3) += momentum_induce_coeff * nH * subz + momentum_hall_coeff * chi_H * sjbz;
                        RHS_H(i, j, k, 4) += momentum_induce_coeff * nH * subu + momentum_hall_coeff * chi_H * sjbu;
                        // RHS_H(i, j, k, 1) -= sns0 * (dpex / ne_cm);
                        // RHS_H(i, j, k, 2) -= sns0 * (dpey / ne_cm);
                        // RHS_H(i, j, k, 3) -= sns0 * (dpez / ne_cm);
                        // RHS_H(i, j, k, 4) -= sns0 * (dpeu / ne_cm);
                        //=============================================================================================

                        //=============================================================================================
                        // Photoionization Source Terms
                        // no source terms for H+
                        // RHS_H(i, j, k, 4) += a6 * sse * Tn0 * chi_H;
                        // This is the electron energy that has been divided into corresponding species
                        // 光致电离产生电子的那一部分电子能量，按照组分数密度分配给对应组分
                        //=============================================================================================

                        //=============================================================================================
                        // Drag Related Source Terms
                        RHS_H(i, j, k, 1) -= a4 * rho * u * vst;
                        RHS_H(i, j, k, 2) -= a4 * rho * v * vst;
                        RHS_H(i, j, k, 3) -= a4 * rho * w * vst;
                        RHS_H(i, j, k, 4) -= a4 * rho * u2 * vst;
                        //=============================================================================================

                        //=============================================================================================
                        // Temperature relation & Thermal Drift
                        RHS_H(i, j, k, 4) += a5 * nH * b1;       // Temperature relation
                        RHS_H(i, j, k, 4) += a4 * rho * u2 * b2; //  Thermal Drift
                        //=============================================================================================
                    }

                    // =====================
                    // species Na+ (ls=2)
                    // =====================
                    {
                        const double rho = std::max(UNa(i, j, k, 0), 0.0);
                        const double u = PVN(i, j, k, 0);
                        const double v = PVN(i, j, k, 1);
                        const double w = PVN(i, j, k, 2);
                        const double u2 = u * u + v * v + w * w;

                        const double subx = (v - upy) * Bz - (w - upz) * By;
                        const double suby = (w - upz) * Bx - (u - upx) * Bz;
                        const double subz = (u - upx) * By - (v - upy) * Bx;

                        const double subu = subx * u + suby * v + subz * w;
                        const double sjbu = sjbx * u + sjby * v + sjbz * w;
                        // const double dpeu = dpex * u + dpey * v + dpez * w;

                        // Ts0 in Kelvin
                        const double Ts0 = PVN(i, j, k, 4) * T_ref;
                        // const double us2 = uN * uN + vN * vN + wN * wN;
                        const double b1 = vst * (Tn0 - Ts0) * (M_Na / (M_Na + M_Na));
                        const double b2 = (M_Na / (M_Na + M_Na)) * vst;

                        //=============================================================================================
                        // Electromagnetic Source Terms
                        RHS_Na(i, j, k, 1) += momentum_induce_coeff * nNa * subx + momentum_hall_coeff * chi_Na * sjbx;
                        RHS_Na(i, j, k, 2) += momentum_induce_coeff * nNa * suby + momentum_hall_coeff * chi_Na * sjby;
                        RHS_Na(i, j, k, 3) += momentum_induce_coeff * nNa * subz + momentum_hall_coeff * chi_Na * sjbz;
                        RHS_Na(i, j, k, 4) += momentum_induce_coeff * nNa * subu + momentum_hall_coeff * chi_Na * sjbu;
                        // RHS_Na(i, j, k, 1) -= sns0 * (dpex / ne_cm);
                        // RHS_Na(i, j, k, 2) -= sns0 * (dpey / ne_cm);
                        // RHS_Na(i, j, k, 3) -= sns0 * (dpez / ne_cm);
                        // RHS_Na(i, j, k, 4) -= sns0 * (dpeu / ne_cm);
                        //=============================================================================================

                        //=============================================================================================
                        // Photoionization Source Terms
                        RHS_Na(i, j, k, 0) += a1_Na * sss;    // Na+ mass creation
                                                              // no Photoionization related source term for momentum eqs
                        RHS_Na(i, j, k, 4) += a6 * sss * Tn0; // Photoionization energy (internal) pump into this species
                        // RHS_Na(i, j, k, 4) += a6 * sse * Tn0 * chi_Na;
                        // This is the electron energy that has been divided into corresponding species
                        // 光致电离产生电子的那一部分电子能量，按照组分数密度分配给对应组分
                        //=============================================================================================

                        //=============================================================================================
                        // Drag Related Source Terms
                        RHS_Na(i, j, k, 1) -= a4 * rho * u * vst;
                        RHS_Na(i, j, k, 2) -= a4 * rho * v * vst;
                        RHS_Na(i, j, k, 3) -= a4 * rho * w * vst;
                        RHS_Na(i, j, k, 4) -= a4 * rho * u2 * vst;
                        //=============================================================================================

                        //=============================================================================================
                        // Temperature relation & Thermal Drift
                        RHS_Na(i, j, k, 4) += a5 * nNa * b1;      // Temperature relation
                        RHS_Na(i, j, k, 4) += a4 * rho * u2 * b2; //  Thermal Drift
                        //=============================================================================================
                    }
                }
    }
}
