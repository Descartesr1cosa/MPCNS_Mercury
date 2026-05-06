#include "MercurySolver.h"
void MercurySolver::AddHallEdgeEMF_()
{
    BuildHallFaceEMF_Rusanov_diff_();

    // // 0) 清零 Ehall_edge / Ehall_face
    // // for (int ib = 0; ib < fld_->num_blocks(); ++ib)
    // // {
    // //     auto zero_scalar_edge = [&](FieldBlock &F)
    // //     {
    // //         if (!F.is_allocated())
    // //             return;
    // //         Int3 lo = F.inner_lo(), hi = F.inner_hi();
    // //         for (int i = lo.i; i < hi.i; ++i)
    // //             for (int j = lo.j; j < hi.j; ++j)
    // //                 for (int k = lo.k; k < hi.k; ++k)
    // //                     F(i, j, k, 0) = 0.0;
    // //     };

    // //     auto zero_vec_face = [&](FieldBlock &F)
    // //     {
    // //         if (!F.is_allocated())
    // //             return;
    // //         Int3 lo = F.inner_lo(), hi = F.inner_hi();
    // //         for (int i = lo.i; i < hi.i; ++i)
    // //             for (int j = lo.j; j < hi.j; ++j)
    // //                 for (int k = lo.k; k < hi.k; ++k)
    // //                     for (int m = 0; m < 3; ++m)
    // //                         F(i, j, k, m) = 0.0;
    // //     };

    // //     zero_scalar_edge(fld_->field(fid_.fid_Ehall.xi, ib));
    // //     zero_scalar_edge(fld_->field(fid_.fid_Ehall.eta, ib));
    // //     zero_scalar_edge(fld_->field(fid_.fid_Ehall.zeta, ib));

    // //     zero_vec_face(fld_->field(fid_.fid_Eface.xi, ib));
    // //     zero_vec_face(fld_->field(fid_.fid_Eface.eta, ib));
    // //     zero_vec_face(fld_->field(fid_.fid_Eface.zeta, ib));
    // // }

    // // 2) face 上做 Hall-Rusanov
    // BuildHallFaceEMF_Rusanov_();

    // // 3) sync Hall face field
    // mercury_bound_.Sync("Eface");

    // // 4) face -> edge
    // AssembleEdgeEMF_FromFaceE_Hall_();

    // // 5) sync edge Hall EMF
    // mercury_bound_.Sync("Ehall");

    // // 6) 加到总 E 上
}

void MercurySolver::BuildHallFaceEMF_Rusanov_diff_()
{
    constexpr double eps = 1e-14;
    const double Cwh = 0.1; // 先沿用你原来的 whistler-LLF 系数

    auto avg4 = [](double a, double b, double c, double d) -> double
    {
        return 0.25 * (a + b + c + d);
    };

    auto max4 = [](double a, double b, double c, double d) -> double
    {
        return std::max(std::max(a, b), std::max(c, d));
    };

    auto norm3 = [](double x, double y, double z) -> double
    {
        return std::sqrt(x * x + y * y + z * z);
    };

    const int nb = fld_->num_blocks();

    for (int ib = 0; ib < nb; ++ib)
    {
        // ---------- cell-centered total B / J / plasma ----------
        auto &Bc = fld_->field(fid_.fid_Bcell, ib); // total B on cell (Bind + Badd)
        auto &Jc = fld_->field(fid_.fid_Jcell, ib);

        auto &UH = fld_->field(fid_.fid_U_H, ib);
        auto &UNa = fld_->field(fid_.fid_U_Na, ib);

        // ---------- face 2-form magnetic flux DOFs (induced only) ----------
        auto &Bxi = fld_->field(fid_.fid_B.xi, ib);
        auto &Bet = fld_->field(fid_.fid_B.eta, ib);
        auto &Bze = fld_->field(fid_.fid_B.zeta, ib);

        // ---------- edge 1-form Hall EMF DOFs ----------
        auto &Exi = fld_->field(fid_.fid_Ehall.xi, ib);
        auto &Eet = fld_->field(fid_.fid_Ehall.eta, ib);
        auto &Eze = fld_->field(fid_.fid_Ehall.zeta, ib);

        // ---------- static geometry fields ----------
        auto &Axi = fld_->field(fid_.Face_Area.xi, ib);
        auto &Aet = fld_->field(fid_.Face_Area.eta, ib);
        auto &Aze = fld_->field(fid_.Face_Area.zeta, ib);

        auto &dl_xi = fld_->field(fid_.Edge_dl.xi, ib);
        auto &dl_et = fld_->field(fid_.Edge_dl.eta, ib);
        auto &dl_ze = fld_->field(fid_.Edge_dl.zeta, ib);

        auto &dr_xi = fld_->field(fid_.Edge_dr.xi, ib);
        auto &dr_et = fld_->field(fid_.Edge_dr.eta, ib);
        auto &dr_ze = fld_->field(fid_.Edge_dr.zeta, ib);

        if (!Bc.is_allocated() || !Jc.is_allocated() ||
            !UH.is_allocated() || !UNa.is_allocated() ||
            !Bxi.is_allocated() || !Bet.is_allocated() || !Bze.is_allocated() ||
            !Exi.is_allocated() || !Eet.is_allocated() || !Eze.is_allocated() ||
            !Axi.is_allocated() || !Aet.is_allocated() || !Aze.is_allocated() ||
            !dl_xi.is_allocated() || !dl_et.is_allocated() || !dl_ze.is_allocated())
        {
            continue;
        }

        auto &buf = hall_face_scratch_[ib];
        auto &Ehc = buf.Ehc;   // cell-centered physical Hall E = alpha (J x B)
        auto &beta = buf.beta; // cell-centered Hall diffusivity scale = |alpha| |B|

        // ============================================================
        // 1) cell 上预计算:
        //    Ehc = alpha * (J x B), beta = |alpha| |B|
        // ============================================================
        {
            const Int3 clo = Bc.get_lo();
            const Int3 chi = Bc.get_hi();

            for (int i = clo.i; i < chi.i; ++i)
                for (int j = clo.j; j < chi.j; ++j)
                    for (int k = clo.k; k < chi.k; ++k)
                    {
                        // double num[3];
                        // Hall_Num_Limiter(UH(i, j, k, 0), UNa(i, j, k, 0), num);
                        // const double ne = std::max(num[2], eps);

                        NumInfo num = Hall_Num_Limiter(UH(i, j, k, 0), UNa(i, j, k, 0));
                        const double ne = num.ne_eff;

                        const double alpha = hall_coef / ne;

                        const double Jx = Jc(i, j, k, 0);
                        const double Jy = Jc(i, j, k, 1);
                        const double Jz = Jc(i, j, k, 2);

                        const double Bx = Bc(i, j, k, 0);
                        const double By = Bc(i, j, k, 1);
                        const double Bz = Bc(i, j, k, 2);

                        Ehc(i, j, k, 0) = alpha * (Jy * Bz - Jz * By);
                        Ehc(i, j, k, 1) = alpha * (Jz * Bx - Jx * Bz);
                        Ehc(i, j, k, 2) = alpha * (Jx * By - Jy * Bx);

                        beta(i, j, k) = std::abs(alpha) * norm3(Bx, By, Bz);
                    }
        }

        // ============================================================
        // 2) xi-edge:
        //
        // Exi(i,j,k) = <Ehc>_4cell · dr_xi
        //            - 0.5 * mu_eta  * (Bzeta(i,j,k)-Bzeta(i,j-1,k))
        //            + 0.5 * mu_zeta * (Beta (i,j,k)-Beta (i,j,k-1))
        //
        // surrounding 4 cells:
        // (i, j-1, k-1), (i, j-1, k), (i, j, k-1), (i, j, k)
        // ============================================================
        {
            Int3 lo = Exi.inner_lo();
            Int3 hi = Exi.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        // ---- central Hall EMF: average 4 cells then line-integrate ----
                        const double Ecx = avg4(Ehc(i, j - 1, k - 1, 0),
                                                Ehc(i, j - 1, k, 0),
                                                Ehc(i, j, k - 1, 0),
                                                Ehc(i, j, k, 0));

                        const double Ecy = avg4(Ehc(i, j - 1, k - 1, 1),
                                                Ehc(i, j - 1, k, 1),
                                                Ehc(i, j, k - 1, 1),
                                                Ehc(i, j, k, 1));

                        const double Ecz = avg4(Ehc(i, j - 1, k - 1, 2),
                                                Ehc(i, j - 1, k, 2),
                                                Ehc(i, j, k - 1, 2),
                                                Ehc(i, j, k, 2));

                        const double Ecen =
                            Ecx * dr_xi(i, j, k, 0) +
                            Ecy * dr_xi(i, j, k, 1) +
                            Ecz * dr_xi(i, j, k, 2);

                        // ---- edge-local Hall diffusivity scale ----
                        const double beta_e = max4(beta(i, j - 1, k - 1),
                                                   beta(i, j - 1, k),
                                                   beta(i, j, k - 1),
                                                   beta(i, j, k));

                        const double L = std::max(dl_xi(i, j, k, 0), eps);

                        // h_eta 由周围 zeta-face area / L_xi 给出
                        const double Abar_zeta = 0.5 * (Aze(i, j, k, 0) +
                                                        Aze(i, j - 1, k, 0));
                        const double inv_h_eta = L / std::max(Abar_zeta, eps);
                        const double mu_eta = Cwh * beta_e * inv_h_eta * inv_h_eta;

                        // h_zeta 由周围 eta-face area / L_xi 给出
                        const double Abar_eta = 0.5 * (Aet(i, j, k, 0) +
                                                       Aet(i, j, k - 1, 0));
                        const double inv_h_zeta = L / std::max(Abar_eta, eps);
                        const double mu_zeta = Cwh * beta_e * inv_h_zeta * inv_h_zeta;

                        // ---- jumps of face 2-form magnetic flux DOFs ----
                        const double dBzeta_eta = Bze(i, j, k, 0) - Bze(i, j - 1, k, 0);
                        const double dBeta_zeta = Bet(i, j, k, 0) - Bet(i, j, k - 1, 0);

                        Exi(i, j, k, 0) =
                            Ecen - 0.5 * mu_eta * dBzeta_eta + 0.5 * mu_zeta * dBeta_zeta;
                    }
        }

        // ============================================================
        // 3) eta-edge:
        //
        // Eeta(i,j,k) = <Ehc>_4cell · dr_eta
        //             - 0.5 * mu_zeta * (Bxi(i,j,k)-Bxi(i,j,k-1))
        //             + 0.5 * mu_xi   * (Bzeta(i,j,k)-Bzeta(i-1,j,k))
        //
        // surrounding 4 cells:
        // (i-1, j, k-1), (i-1, j, k), (i, j, k-1), (i, j, k)
        // ============================================================
        {
            Int3 lo = Eet.inner_lo();
            Int3 hi = Eet.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const double Ecx = avg4(Ehc(i - 1, j, k - 1, 0),
                                                Ehc(i - 1, j, k, 0),
                                                Ehc(i, j, k - 1, 0),
                                                Ehc(i, j, k, 0));

                        const double Ecy = avg4(Ehc(i - 1, j, k - 1, 1),
                                                Ehc(i - 1, j, k, 1),
                                                Ehc(i, j, k - 1, 1),
                                                Ehc(i, j, k, 1));

                        const double Ecz = avg4(Ehc(i - 1, j, k - 1, 2),
                                                Ehc(i - 1, j, k, 2),
                                                Ehc(i, j, k - 1, 2),
                                                Ehc(i, j, k, 2));

                        const double Ecen =
                            Ecx * dr_et(i, j, k, 0) +
                            Ecy * dr_et(i, j, k, 1) +
                            Ecz * dr_et(i, j, k, 2);

                        const double beta_e = max4(beta(i - 1, j, k - 1),
                                                   beta(i - 1, j, k),
                                                   beta(i, j, k - 1),
                                                   beta(i, j, k));

                        const double L = std::max(dl_et(i, j, k, 0), eps);

                        // h_zeta from xi-face area / L_eta
                        const double Abar_xi = 0.5 * (Axi(i, j, k, 0) +
                                                      Axi(i, j, k - 1, 0));
                        const double inv_h_zeta = L / std::max(Abar_xi, eps);
                        const double mu_zeta = Cwh * beta_e * inv_h_zeta * inv_h_zeta;

                        // h_xi from zeta-face area / L_eta
                        const double Abar_zeta = 0.5 * (Aze(i, j, k, 0) +
                                                        Aze(i - 1, j, k, 0));
                        const double inv_h_xi = L / std::max(Abar_zeta, eps);
                        const double mu_xi = Cwh * beta_e * inv_h_xi * inv_h_xi;

                        const double dBxi_zeta = Bxi(i, j, k, 0) - Bxi(i, j, k - 1, 0);
                        const double dBzeta_xi = Bze(i, j, k, 0) - Bze(i - 1, j, k, 0);

                        Eet(i, j, k, 0) =
                            Ecen - 0.5 * mu_zeta * dBxi_zeta + 0.5 * mu_xi * dBzeta_xi;
                    }
        }

        // ============================================================
        // 4) zeta-edge:
        //
        // Ezeta(i,j,k) = <Ehc>_4cell · dr_zeta
        //              - 0.5 * mu_xi  * (Beta(i,j,k)-Beta(i-1,j,k))
        //              + 0.5 * mu_eta * (Bxi (i,j,k)-Bxi (i,j-1,k))
        //
        // surrounding 4 cells:
        // (i-1, j-1, k), (i-1, j, k), (i, j-1, k), (i, j, k)
        // ============================================================
        {
            Int3 lo = Eze.inner_lo();
            Int3 hi = Eze.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const double Ecx = avg4(Ehc(i - 1, j - 1, k, 0),
                                                Ehc(i - 1, j, k, 0),
                                                Ehc(i, j - 1, k, 0),
                                                Ehc(i, j, k, 0));

                        const double Ecy = avg4(Ehc(i - 1, j - 1, k, 1),
                                                Ehc(i - 1, j, k, 1),
                                                Ehc(i, j - 1, k, 1),
                                                Ehc(i, j, k, 1));

                        const double Ecz = avg4(Ehc(i - 1, j - 1, k, 2),
                                                Ehc(i - 1, j, k, 2),
                                                Ehc(i, j - 1, k, 2),
                                                Ehc(i, j, k, 2));

                        const double Ecen =
                            Ecx * dr_ze(i, j, k, 0) +
                            Ecy * dr_ze(i, j, k, 1) +
                            Ecz * dr_ze(i, j, k, 2);

                        const double beta_e = max4(beta(i - 1, j - 1, k),
                                                   beta(i - 1, j, k),
                                                   beta(i, j - 1, k),
                                                   beta(i, j, k));

                        const double L = std::max(dl_ze(i, j, k, 0), eps);

                        // h_xi from eta-face area / L_zeta
                        const double Abar_eta = 0.5 * (Aet(i, j, k, 0) +
                                                       Aet(i - 1, j, k, 0));
                        const double inv_h_xi = L / std::max(Abar_eta, eps);
                        const double mu_xi = Cwh * beta_e * inv_h_xi * inv_h_xi;

                        // h_eta from xi-face area / L_zeta
                        const double Abar_xi = 0.5 * (Axi(i, j, k, 0) +
                                                      Axi(i, j - 1, k, 0));
                        const double inv_h_eta = L / std::max(Abar_xi, eps);
                        const double mu_eta = Cwh * beta_e * inv_h_eta * inv_h_eta;

                        const double dBeta_xi = Bet(i, j, k, 0) - Bet(i - 1, j, k, 0);
                        const double dBxi_eta = Bxi(i, j, k, 0) - Bxi(i, j - 1, k, 0);

                        Eze(i, j, k, 0) =
                            Ecen - 0.5 * mu_xi * dBeta_xi + 0.5 * mu_eta * dBxi_eta;
                    }
        }

        // ------------------------------------------------------------
        // hard zero Hall edge EMF on exact Pole edges
        // place this AFTER Exi/Eet/Eze have all been assembled for block ib
        // ------------------------------------------------------------
        auto zero_edge_box = [](FieldBlock &E, StaggerLocation loc, const Box3 &node_box)
        {
            if (!E.is_allocated())
                return;

            Int3 lo = node_box.lo;
            Int3 hi = node_box.hi;

            // node-box -> edge dof box
            if (loc == StaggerLocation::EdgeXi)
            {
                hi.i -= 1;
            }
            else if (loc == StaggerLocation::EdgeEt)
            {
                hi.j -= 1;
            }
            else if (loc == StaggerLocation::EdgeZe)
            {
                hi.k -= 1;
            }

            Int3 elo = E.inner_lo();
            Int3 ehi = E.inner_hi();

            lo.i = std::max(lo.i, elo.i);
            lo.j = std::max(lo.j, elo.j);
            lo.k = std::max(lo.k, elo.k);

            hi.i = std::min(hi.i, ehi.i);
            hi.j = std::min(hi.j, ehi.j);
            hi.k = std::min(hi.k, ehi.k);

            if (!(lo.i < hi.i && lo.j < hi.j && lo.k < hi.k))
                return;

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        E(i, j, k, 0) = 0.0;
        };

        auto zero_edge_box_inward = [](FieldBlock &E,
                                       StaggerLocation loc,
                                       const Box3 &node_box,
                                       int direction,
                                       int nlayer)
        {
            if (!E.is_allocated())
                return;

            Int3 lo = node_box.lo;
            Int3 hi = node_box.hi;

            const int ax = std::abs(direction) - 1;

            // 沿法向向内扩展 nlayer 层（node box）
            if (direction < 0)
            {
                if (ax == 0)
                    hi.i += nlayer;
                if (ax == 1)
                    hi.j += nlayer;
                if (ax == 2)
                    hi.k += nlayer;
            }
            else
            {
                if (ax == 0)
                    lo.i -= nlayer;
                if (ax == 1)
                    lo.j -= nlayer;
                if (ax == 2)
                    lo.k -= nlayer;
            }

            // node-box -> edge dof box
            if (loc == StaggerLocation::EdgeXi)
            {
                hi.i -= 1;
            }
            else if (loc == StaggerLocation::EdgeEt)
            {
                hi.j -= 1;
            }
            else if (loc == StaggerLocation::EdgeZe)
            {
                hi.k -= 1;
            }

            Int3 elo = E.inner_lo();
            Int3 ehi = E.inner_hi();

            lo.i = std::max(lo.i, elo.i);
            lo.j = std::max(lo.j, elo.j);
            lo.k = std::max(lo.k, elo.k);

            hi.i = std::min(hi.i, ehi.i);
            hi.j = std::min(hi.j, ehi.j);
            hi.k = std::min(hi.k, ehi.k);

            if (!(lo.i < hi.i && lo.j < hi.j && lo.k < hi.k))
                return;

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        E(i, j, k, 0) = 0.0;
        };

        for (const auto &p : topo_->physical_patches)
        {
            if (p.this_block != ib)
                continue;

            if (p.bc_name != "Solid_Surface" && p.bc_name != "Coupled-Solid") // 换成你的壁面 patch 名
                continue;

            const int nlayer = 20; // 先试 1；不够再试 2

            zero_edge_box_inward(Exi, StaggerLocation::EdgeXi, p.this_box_node, p.direction, nlayer);
            zero_edge_box_inward(Eet, StaggerLocation::EdgeEt, p.this_box_node, p.direction, nlayer);
            zero_edge_box_inward(Eze, StaggerLocation::EdgeZe, p.this_box_node, p.direction, nlayer);
        }

        for (const auto &p : topo_->physical_patches)
        {
            if (p.this_block != ib)
                continue;
            if (p.bc_name != "Pole")
                continue;

            // Pole normal = xi  -> zero eta-edge electric field
            if (std::abs(p.direction) == 1)
            {
                zero_edge_box(Eet, StaggerLocation::EdgeEt, p.this_box_node);
            }

            // Pole normal = eta -> zero xi-edge electric field
            if (std::abs(p.direction) == 2)
            {
                zero_edge_box(Exi, StaggerLocation::EdgeXi, p.this_box_node);
            }

            // 如果你想更狠一点，也可以把 Pole 上的 zeta-edge 一并关掉：
            zero_edge_box(Eze, StaggerLocation::EdgeZe, p.this_box_node);
        }
    }
}

void MercurySolver::BuildHallFaceEMF_Rusanov_()
{
    constexpr double eps = 1e-14;
    const double Cwh = 0.5;
    const double hall_coef_abs = std::abs(hall_coef);

    auto norm3 = [](double x, double y, double z) -> double
    {
        return std::sqrt(x * x + y * y + z * z);
    };

    const int nb = fld_->num_blocks();

    for (int ib = 0; ib < nb; ++ib)
    {
        auto &Bc = fld_->field(fid_.fid_Bcell, ib); // total B on cells
        auto &Binduce = fld_->field(fid_.fid_Bindcell, ib);
        auto &Jc = fld_->field(fid_.fid_Jcell, ib);

        auto &UH = fld_->field(fid_.fid_U_H, ib);
        auto &UNa = fld_->field(fid_.fid_U_Na, ib);

        auto &Efxi = fld_->field(fid_.fid_Eface.xi, ib);
        auto &Efet = fld_->field(fid_.fid_Eface.eta, ib);
        auto &Efze = fld_->field(fid_.fid_Eface.zeta, ib);

        auto &JDxi = fld_->field(fid_.fid_metric.xi, ib);
        auto &JDet = fld_->field(fid_.fid_metric.eta, ib);
        auto &JDze = fld_->field(fid_.fid_metric.zeta, ib);

        auto &dlst_xi = fld_->field("dlstar_xi", ib);
        auto &dlst_et = fld_->field("dlstar_eta", ib);
        auto &dlst_ze = fld_->field("dlstar_zeta", ib);

        if (!Bc.is_allocated() || !Binduce.is_allocated() || !Jc.is_allocated() ||
            !UH.is_allocated() || !UNa.is_allocated() ||
            !Efxi.is_allocated() || !Efet.is_allocated() || !Efze.is_allocated())
        {
            continue;
        }

        auto &buf = hall_face_scratch_[ib];
        auto &Ehc = buf.Ehc;        // Ehc(i,j,k,comp)
        auto &beta_hall = buf.beta; // beta_hall(i,j,k)

        const Int3 clo = Bc.get_lo();
        const Int3 chi = Bc.get_hi();

        // ============================================================
        // 1) cell 上预计算 Ehc = alpha * (J x B), beta_hall
        // ============================================================
        for (int i = clo.i; i < chi.i; ++i)
            for (int j = clo.j; j < chi.j; ++j)
                for (int k = clo.k; k < chi.k; ++k)
                {
                    // double num[3];
                    // Hall_Num_Limiter(UH(i, j, k, 0), UNa(i, j, k, 0), num);
                    // const double ne = num[2];

                    NumInfo num = Hall_Num_Limiter(UH(i, j, k, 0), UNa(i, j, k, 0));
                    const double ne = num.ne_eff;

                    const double alpha = hall_coef / ne;

                    const double Jx = Jc(i, j, k, 0);
                    const double Jy = Jc(i, j, k, 1);
                    const double Jz = Jc(i, j, k, 2);

                    const double Bx = Bc(i, j, k, 0);
                    const double By = Bc(i, j, k, 1);
                    const double Bz = Bc(i, j, k, 2);

                    Ehc(i, j, k, 0) = alpha * (Jy * Bz - Jz * By);
                    Ehc(i, j, k, 1) = alpha * (Jz * Bx - Jx * Bz);
                    Ehc(i, j, k, 2) = alpha * (Jx * By - Jy * Bx);

                    beta_hall(i, j, k) = hall_coef_abs * norm3(Bx, By, Bz) / ne;
                }

        // ============================================================
        // 2) xi-face
        // ============================================================
        {
            Int3 lo = Efxi.inner_lo();
            Int3 hi = Efxi.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const int iL = i - 1;
                        const int iR = i;

                        double ELx = Ehc(iL, j, k, 0);
                        double ELy = Ehc(iL, j, k, 1);
                        double ELz = Ehc(iL, j, k, 2);

                        double ERx = Ehc(iR, j, k, 0);
                        double ERy = Ehc(iR, j, k, 1);
                        double ERz = Ehc(iR, j, k, 2);

                        const double Sx = JDxi(i, j, k, 0);
                        const double Sy = JDxi(i, j, k, 1);
                        const double Sz = JDxi(i, j, k, 2);

                        const double Smag = norm3(Sx, Sy, Sz) + eps;
                        const double nx = Sx / Smag;
                        const double ny = Sy / Smag;
                        const double nz = Sz / Smag;

                        {
                            const double En = ELx * nx + ELy * ny + ELz * nz;
                            ELx -= En * nx;
                            ELy -= En * ny;
                            ELz -= En * nz;
                        }
                        {
                            const double En = ERx * nx + ERy * ny + ERz * nz;
                            ERx -= En * nx;
                            ERy -= En * ny;
                            ERz -= En * nz;
                        }

                        const double h_n = std::max(dlst_xi(i, j, k, 0), eps);
                        const double sH = Cwh * std::max(beta_hall(iL, j, k), beta_hall(iR, j, k)) / h_n;

                        const double BLx = Binduce(iL, j, k, 0);
                        const double BLy = Binduce(iL, j, k, 1);
                        const double BLz = Binduce(iL, j, k, 2);

                        const double BRx = Binduce(iR, j, k, 0);
                        const double BRy = Binduce(iR, j, k, 1);
                        const double BRz = Binduce(iR, j, k, 2);

                        const double dBx = BRx - BLx;
                        const double dBy = BRy - BLy;
                        const double dBz = BRz - BLz;

                        const double cx = ny * dBz - nz * dBy;
                        const double cy = nz * dBx - nx * dBz;
                        const double cz = nx * dBy - ny * dBx;

                        Efxi(i, j, k, 0) = 0.5 * (ELx + ERx) + 0.5 * sH * cx;
                        Efxi(i, j, k, 1) = 0.5 * (ELy + ERy) + 0.5 * sH * cy;
                        Efxi(i, j, k, 2) = 0.5 * (ELz + ERz) + 0.5 * sH * cz;
                    }
        }

        // ============================================================
        // 3) eta-face
        // ============================================================
        {
            Int3 lo = Efet.inner_lo();
            Int3 hi = Efet.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const int jL = j - 1;
                        const int jR = j;

                        double ELx = Ehc(i, jL, k, 0);
                        double ELy = Ehc(i, jL, k, 1);
                        double ELz = Ehc(i, jL, k, 2);

                        double ERx = Ehc(i, jR, k, 0);
                        double ERy = Ehc(i, jR, k, 1);
                        double ERz = Ehc(i, jR, k, 2);

                        const double Sx = JDet(i, j, k, 0);
                        const double Sy = JDet(i, j, k, 1);
                        const double Sz = JDet(i, j, k, 2);

                        const double Smag = norm3(Sx, Sy, Sz) + eps;
                        const double nx = Sx / Smag;
                        const double ny = Sy / Smag;
                        const double nz = Sz / Smag;

                        {
                            const double En = ELx * nx + ELy * ny + ELz * nz;
                            ELx -= En * nx;
                            ELy -= En * ny;
                            ELz -= En * nz;
                        }
                        {
                            const double En = ERx * nx + ERy * ny + ERz * nz;
                            ERx -= En * nx;
                            ERy -= En * ny;
                            ERz -= En * nz;
                        }

                        const double h_n = std::max(dlst_et(i, j, k, 0), eps);
                        const double sH = Cwh * std::max(beta_hall(i, jL, k), beta_hall(i, jR, k)) / h_n;

                        const double BLx = Binduce(i, jL, k, 0);
                        const double BLy = Binduce(i, jL, k, 1);
                        const double BLz = Binduce(i, jL, k, 2);

                        const double BRx = Binduce(i, jR, k, 0);
                        const double BRy = Binduce(i, jR, k, 1);
                        const double BRz = Binduce(i, jR, k, 2);

                        const double dBx = BRx - BLx;
                        const double dBy = BRy - BLy;
                        const double dBz = BRz - BLz;

                        const double cx = ny * dBz - nz * dBy;
                        const double cy = nz * dBx - nx * dBz;
                        const double cz = nx * dBy - ny * dBx;

                        Efet(i, j, k, 0) = 0.5 * (ELx + ERx) + 0.5 * sH * cx;
                        Efet(i, j, k, 1) = 0.5 * (ELy + ERy) + 0.5 * sH * cy;
                        Efet(i, j, k, 2) = 0.5 * (ELz + ERz) + 0.5 * sH * cz;
                    }
        }

        // ============================================================
        // 4) zeta-face
        // ============================================================
        {
            Int3 lo = Efze.inner_lo();
            Int3 hi = Efze.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const int kL = k - 1;
                        const int kR = k;

                        double ELx = Ehc(i, j, kL, 0);
                        double ELy = Ehc(i, j, kL, 1);
                        double ELz = Ehc(i, j, kL, 2);

                        double ERx = Ehc(i, j, kR, 0);
                        double ERy = Ehc(i, j, kR, 1);
                        double ERz = Ehc(i, j, kR, 2);

                        const double Sx = JDze(i, j, k, 0);
                        const double Sy = JDze(i, j, k, 1);
                        const double Sz = JDze(i, j, k, 2);

                        const double Smag = norm3(Sx, Sy, Sz) + eps;
                        const double nx = Sx / Smag;
                        const double ny = Sy / Smag;
                        const double nz = Sz / Smag;

                        {
                            const double En = ELx * nx + ELy * ny + ELz * nz;
                            ELx -= En * nx;
                            ELy -= En * ny;
                            ELz -= En * nz;
                        }
                        {
                            const double En = ERx * nx + ERy * ny + ERz * nz;
                            ERx -= En * nx;
                            ERy -= En * ny;
                            ERz -= En * nz;
                        }

                        const double h_n = std::max(dlst_ze(i, j, k, 0), eps);
                        const double sH = Cwh * std::max(beta_hall(i, j, kL), beta_hall(i, j, kR)) / h_n;

                        const double BLx = Binduce(i, j, kL, 0);
                        const double BLy = Binduce(i, j, kL, 1);
                        const double BLz = Binduce(i, j, kL, 2);

                        const double BRx = Binduce(i, j, kR, 0);
                        const double BRy = Binduce(i, j, kR, 1);
                        const double BRz = Binduce(i, j, kR, 2);

                        const double dBx = BRx - BLx;
                        const double dBy = BRy - BLy;
                        const double dBz = BRz - BLz;

                        const double cx = ny * dBz - nz * dBy;
                        const double cy = nz * dBx - nx * dBz;
                        const double cz = nx * dBy - ny * dBx;

                        Efze(i, j, k, 0) = 0.5 * (ELx + ERx) + 0.5 * sH * cx;
                        Efze(i, j, k, 1) = 0.5 * (ELy + ERy) + 0.5 * sH * cy;
                        Efze(i, j, k, 2) = 0.5 * (ELz + ERz) + 0.5 * sH * cz;
                    }
        }
    }
}

void MercurySolver::AssembleEdgeEMF_FromFaceE_Hall_()
{
    Int3 sub, sup;

    for (int iblk = 0; iblk < fld_->num_blocks(); ++iblk)
    {
        double3D &x = fld_->grd->grids(iblk).x;
        double3D &y = fld_->grd->grids(iblk).y;
        double3D &z = fld_->grd->grids(iblk).z;

        Vec3 E, dr;

        // -----------------------------------------
        // Exi(edge) from face-eta and face-zeta
        // -----------------------------------------
        {
            auto &Exi = fld_->field(fid_.fid_Ehall.xi, iblk);
            auto &Ef_eta = fld_->field(fid_.fid_Eface.eta, iblk);
            auto &Ef_ze = fld_->field(fid_.fid_Eface.zeta, iblk);

            sub = Exi.inner_lo();
            sup = Exi.inner_hi();

            for (int i = sub.i; i < sup.i; ++i)
                for (int j = sub.j; j < sup.j; ++j)
                    for (int k = sub.k; k < sup.k; ++k)
                    {
                        E.vec[0] = 0.25 * (Ef_eta(i, j, k, 0) + Ef_eta(i, j, k - 1, 0) + Ef_ze(i, j, k, 0) + Ef_ze(i, j - 1, k, 0));
                        E.vec[1] = 0.25 * (Ef_eta(i, j, k, 1) + Ef_eta(i, j, k - 1, 1) + Ef_ze(i, j, k, 1) + Ef_ze(i, j - 1, k, 1));
                        E.vec[2] = 0.25 * (Ef_eta(i, j, k, 2) + Ef_eta(i, j, k - 1, 2) + Ef_ze(i, j, k, 2) + Ef_ze(i, j - 1, k, 2));

                        dr = {
                            x(i + 1, j, k) - x(i, j, k),
                            y(i + 1, j, k) - y(i, j, k),
                            z(i + 1, j, k) - z(i, j, k)};

                        Exi(i, j, k, 0) = E * dr;
                    }
        }

        // -----------------------------------------
        // Eeta(edge) from face-xi and face-zeta
        // -----------------------------------------
        {
            auto &Eeta = fld_->field(fid_.fid_Ehall.eta, iblk);
            auto &Ef_xi = fld_->field(fid_.fid_Eface.xi, iblk);
            auto &Ef_ze = fld_->field(fid_.fid_Eface.zeta, iblk);

            sub = Eeta.inner_lo();
            sup = Eeta.inner_hi();

            for (int i = sub.i; i < sup.i; ++i)
                for (int j = sub.j; j < sup.j; ++j)
                    for (int k = sub.k; k < sup.k; ++k)
                    {
                        E.vec[0] = 0.25 * (Ef_xi(i, j, k, 0) + Ef_xi(i, j, k - 1, 0) + Ef_ze(i, j, k, 0) + Ef_ze(i - 1, j, k, 0));
                        E.vec[1] = 0.25 * (Ef_xi(i, j, k, 1) + Ef_xi(i, j, k - 1, 1) + Ef_ze(i, j, k, 1) + Ef_ze(i - 1, j, k, 1));
                        E.vec[2] = 0.25 * (Ef_xi(i, j, k, 2) + Ef_xi(i, j, k - 1, 2) + Ef_ze(i, j, k, 2) + Ef_ze(i - 1, j, k, 2));

                        dr = {
                            x(i, j + 1, k) - x(i, j, k),
                            y(i, j + 1, k) - y(i, j, k),
                            z(i, j + 1, k) - z(i, j, k)};

                        Eeta(i, j, k, 0) = E * dr;
                    }
        }

        // -----------------------------------------
        // Ezeta(edge) from face-xi and face-eta
        // -----------------------------------------
        {
            auto &Eze = fld_->field(fid_.fid_Ehall.zeta, iblk);
            auto &Ef_xi = fld_->field(fid_.fid_Eface.xi, iblk);
            auto &Ef_et = fld_->field(fid_.fid_Eface.eta, iblk);

            sub = Eze.inner_lo();
            sup = Eze.inner_hi();

            for (int i = sub.i; i < sup.i; ++i)
                for (int j = sub.j; j < sup.j; ++j)
                    for (int k = sub.k; k < sup.k; ++k)
                    {
                        E.vec[0] = 0.25 * (Ef_xi(i, j, k, 0) + Ef_xi(i, j - 1, k, 0) + Ef_et(i, j, k, 0) + Ef_et(i - 1, j, k, 0));
                        E.vec[1] = 0.25 * (Ef_xi(i, j, k, 1) + Ef_xi(i, j - 1, k, 1) + Ef_et(i, j, k, 1) + Ef_et(i - 1, j, k, 1));
                        E.vec[2] = 0.25 * (Ef_xi(i, j, k, 2) + Ef_xi(i, j - 1, k, 2) + Ef_et(i, j, k, 2) + Ef_et(i - 1, j, k, 2));

                        dr = {
                            x(i, j, k + 1) - x(i, j, k),
                            y(i, j, k + 1) - y(i, j, k),
                            z(i, j, k + 1) - z(i, j, k)};

                        Eze(i, j, k, 0) = E * dr;
                    }
        }
    }
}

void MercurySolver::ApplyUpdate_Euler_BfaceOnly_(double dt_sub)
{
    const int nb = fld_->num_blocks();
    for (int ib = 0; ib < nb; ++ib)
    {
        auto &Ub_xi = fld_->field(fid_.fid_B.xi, ib);
        auto &Ub_eta = fld_->field(fid_.fid_B.eta, ib);
        auto &Ub_zeta = fld_->field(fid_.fid_B.zeta, ib);

        auto &RHSB_xi = fld_->field(fid_.fid_RHS_b.xi, ib);
        auto &RHSB_eta = fld_->field(fid_.fid_RHS_b.eta, ib);
        auto &RHSB_zeta = fld_->field(fid_.fid_RHS_b.zeta, ib);

        if (!Ub_xi.is_allocated())
            continue;

        {
            Int3 lo = Ub_xi.inner_lo(), hi = Ub_xi.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        Ub_xi(i, j, k, 0) += dt_sub * RHSB_xi(i, j, k, 0);
        }
        {
            Int3 lo = Ub_eta.inner_lo(), hi = Ub_eta.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        Ub_eta(i, j, k, 0) += dt_sub * RHSB_eta(i, j, k, 0);
        }
        {
            Int3 lo = Ub_zeta.inner_lo(), hi = Ub_zeta.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        Ub_zeta(i, j, k, 0) += dt_sub * RHSB_zeta(i, j, k, 0);
        }
    }
}
