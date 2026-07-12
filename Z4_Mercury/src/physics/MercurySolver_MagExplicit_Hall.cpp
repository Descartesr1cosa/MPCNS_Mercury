#include "MercurySolver.h"
#include "1_grid/1_MPCNS_Grid.h"

#include <cmath>

double MercurySolver::HallRadialTaper_(double x, double y, double z)
{
    const double r = std::sqrt(x*x + y*y + z*z);
    const double width = hall_taper_r_max - hall_taper_r_min;
    if (width <= 0.0) return (r >= hall_taper_r_min) ? 1.0 : 0.0;
    const double mid = 0.5 * (hall_taper_r_min + hall_taper_r_max);
    return 0.5 * (1.0 + std::tanh(2.0 * (r - mid) / width));
}

double MercurySolver::HallRadialTaperEdge_(int ib, StaggerLocation loc, int i, int j, int k)
{
    auto &b = grd_->grids(ib);
    int ip=i, jp=j, kp=k;
    if (loc == StaggerLocation::EdgeXi) ++ip;
    else if (loc == StaggerLocation::EdgeEt) ++jp;
    else ++kp;
    return HallRadialTaper_(0.5*(b.x(i,j,k)+b.x(ip,jp,kp)),
                            0.5*(b.y(i,j,k)+b.y(ip,jp,kp)),
                            0.5*(b.z(i,j,k)+b.z(ip,jp,kp)));
}
void MercurySolver::AddHallEdgeEMF_()
{
#if HALL_IMPLICIT == 1
    auto zero_one = [](FieldBlock &F)
    {
        if (!F.is_allocated())
            return;

        Int3 lo = F.get_lo();
        Int3 hi = F.get_hi();
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                    F(i, j, k, 0) = 0.0;
    };

    for (int ib = 0; ib < fld_->num_blocks(); ++ib)
    {
        zero_one(fld_->field(fid_.fid_Ehall.xi, ib));
        zero_one(fld_->field(fid_.fid_Ehall.eta, ib));
        zero_one(fld_->field(fid_.fid_Ehall.zeta, ib));
    }
#endif

    BuildHallFaceEMF_Rusanov_diff_();
#if HALL_IMPLICIT == 0
    AssembleSingularEdgeEMF_HallExplicit_();
#endif
}

void MercurySolver::AssembleSingularEdgeEMF_HallExplicit_()
{
    if (!singular_edges_ || singular_edges_->empty()) return;
    auto contribution=[this](const METRIC::SingularPhysicalEdge &edge,
                             const METRIC::WeightedIncidentEntity &inc)->double
    {
        const auto &c=inc.entity;
        FieldBlock &bc=fld_->field(fid_.fid_Bcell,c.block);
        FieldBlock &jc=fld_->field(fid_.fid_Jcell,c.block);
        FieldBlock &uh=fld_->field(fid_.fid_U_H,c.block);
        FieldBlock &un=fld_->field(fid_.fid_U_Na,c.block);
        if(!bc.is_allocated()||!jc.is_allocated()||!uh.is_allocated()||!un.is_allocated()) return 0.0;
        const NumInfo num=Hall_Num_Limiter(uh(c.i,c.j,c.k,0),un(c.i,c.j,c.k,0));
        const double alpha=hall_coef/num.ne_eff;
        const double ex=alpha*(jc(c.i,c.j,c.k,1)*bc(c.i,c.j,c.k,2)-jc(c.i,c.j,c.k,2)*bc(c.i,c.j,c.k,1));
        const double ey=alpha*(jc(c.i,c.j,c.k,2)*bc(c.i,c.j,c.k,0)-jc(c.i,c.j,c.k,0)*bc(c.i,c.j,c.k,2));
        const double ez=alpha*(jc(c.i,c.j,c.k,0)*bc(c.i,c.j,c.k,1)-jc(c.i,c.j,c.k,1)*bc(c.i,c.j,c.k,0));
        auto &cx=grd_->grids(c.block).dual_x; auto &cy=grd_->grids(c.block).dual_y; auto &cz=grd_->grids(c.block).dual_z;
        const double taper=HallRadialTaper_(cx(c.i+1,c.j+1,c.k+1),cy(c.i+1,c.j+1,c.k+1),cz(c.i+1,c.j+1,c.k+1));
        return taper*(ex*edge.canonical_edge_vector[0]+ey*edge.canonical_edge_vector[1]+ez*edge.canonical_edge_vector[2]);
    };
    singular_edges_->assemble_cell_field_to_local_owners(*fld_,"Ehall_xi",contribution);
    singular_edges_->assemble_cell_field_to_local_owners(*fld_,"Ehall_eta",contribution);
    singular_edges_->assemble_cell_field_to_local_owners(*fld_,"Ehall_zeta",contribution);
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
        //            + 0.5 * mu_eta  * (Bzeta(i,j,k)-Bzeta(i,j-1,k))
        //            - 0.5 * mu_zeta * (Beta (i,j,k)-Beta (i,j,k-1))
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
                            HallRadialTaperEdge_(ib, StaggerLocation::EdgeXi, i,j,k) *
                            (Ecen + 0.5 * mu_eta * dBzeta_eta - 0.5 * mu_zeta * dBeta_zeta);
                    }
        }

        // ============================================================
        // 3) eta-edge:
        //
        // Eeta(i,j,k) = <Ehc>_4cell · dr_eta
        //             + 0.5 * mu_zeta * (Bxi(i,j,k)-Bxi(i,j,k-1))
        //             - 0.5 * mu_xi   * (Bzeta(i,j,k)-Bzeta(i-1,j,k))
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
                            HallRadialTaperEdge_(ib, StaggerLocation::EdgeEt, i,j,k) *
                            (Ecen + 0.5 * mu_zeta * dBxi_zeta - 0.5 * mu_xi * dBzeta_xi);
                    }
        }

        // ============================================================
        // 4) zeta-edge:
        //
        // Ezeta(i,j,k) = <Ehc>_4cell · dr_zeta
        //              + 0.5 * mu_xi  * (Beta(i,j,k)-Beta(i-1,j,k))
        //              - 0.5 * mu_eta * (Bxi (i,j,k)-Bxi (i,j-1,k))
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
                            HallRadialTaperEdge_(ib, StaggerLocation::EdgeZe, i,j,k) *
                            (Ecen + 0.5 * mu_xi * dBeta_xi - 0.5 * mu_eta * dBxi_eta);
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

void MercurySolver::ApplyUpdate_Euler_BfaceOnly_(double dt_sub, const IdTriplet &fid_RHSB)
{
    const int nb = fld_->num_blocks();
    for (int ib = 0; ib < nb; ++ib)
    {
        auto &Ub_xi = fld_->field(fid_.fid_B.xi, ib);
        auto &Ub_eta = fld_->field(fid_.fid_B.eta, ib);
        auto &Ub_zeta = fld_->field(fid_.fid_B.zeta, ib);

        auto &RHSB_xi = fld_->field(fid_RHSB.xi, ib);
        auto &RHSB_eta = fld_->field(fid_RHSB.eta, ib);
        auto &RHSB_zeta = fld_->field(fid_RHSB.zeta, ib);

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
