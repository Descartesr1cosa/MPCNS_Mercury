#include "MercurySolver.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

void MercurySolver::AddIdealEdgeEMF_()
{
    for (int iblk = 0; iblk < fld_->num_blocks(); iblk++)
    {
        auto &Uplus = fld_->field(fid_.fid_U_plus, iblk);
        auto &UH = fld_->field(fid_.fid_U_H, iblk);
        auto &UN = fld_->field(fid_.fid_U_Na, iblk);
        if (!Uplus.is_allocated() || !UH.is_allocated() || !UN.is_allocated())
            continue;

        auto &Bxi = fld_->field(fid_.fid_B.xi, iblk);
        auto &Beta = fld_->field(fid_.fid_B.eta, iblk);
        auto &Bzeta = fld_->field(fid_.fid_B.zeta, iblk);
        auto &Badd_xi = fld_->field(fid_.fid_Badd.xi, iblk);
        auto &Badd_eta = fld_->field(fid_.fid_Badd.eta, iblk);
        auto &Badd_zeta = fld_->field(fid_.fid_Badd.zeta, iblk);

        auto &Jac = fld_->field(fid_.fid_Jac, iblk);
        auto &JDxi = fld_->field(fid_.fid_metric.xi, iblk);
        auto &JDet = fld_->field(fid_.fid_metric.eta, iblk);
        auto &JDze = fld_->field(fid_.fid_metric.zeta, iblk);

        AssembleOneDirectionEMF_(1, fld_->field(fid_.fid_Eface.xi, iblk),
                                 Bxi, Beta, Bzeta,
                                 Badd_xi, Badd_eta, Badd_zeta,
                                 Jac, JDxi, JDet, JDze, Uplus);
        AssembleOneDirectionEMF_(2, fld_->field(fid_.fid_Eface.eta, iblk),
                                 Bxi, Beta, Bzeta,
                                 Badd_xi, Badd_eta, Badd_zeta,
                                 Jac, JDxi, JDet, JDze, Uplus);
        AssembleOneDirectionEMF_(3, fld_->field(fid_.fid_Eface.zeta, iblk),
                                 Bxi, Beta, Bzeta,
                                 Badd_xi, Badd_eta, Badd_zeta,
                                 Jac, JDxi, JDet, JDze, Uplus);
    }

    mercury_bound_.Sync("Eface");

    AssembleEdgeEMF_FromFaceE_Ideal_();
}

void MercurySolver::AssembleSingularEdgeEMF_NonHall_()
{
    if (!singular_edges_ || singular_edges_->empty()) return;
    BuildStationaryWallSingularEdgeSet_();

    auto contribution = [this](const METRIC::SingularPhysicalEdge &edge,
                               const METRIC::WeightedIncidentEntity &inc) -> double
    {
        const auto &c=inc.entity;
        FieldBlock &u=fld_->field(fid_.fid_U_plus,c.block);
        FieldBlock &b=fld_->field(fid_.fid_Bcell,c.block);
        if (!u.is_allocated() || !b.is_allocated()) return 0.0;

        double emf=0.0;
        // The singular owner assembly replaces the block-local candidate.
        // Preserve a stationary-wall tangential Ideal constraint here rather
        // than relying on call order; non-Ideal contributions below remain.
        if (stationary_wall_singular_edge_gids_.find(edge.global_id) ==
            stationary_wall_singular_edge_gids_.end())
        {
            const double ex=-(u(c.i,c.j,c.k,1)*b(c.i,c.j,c.k,2)-
                              u(c.i,c.j,c.k,2)*b(c.i,c.j,c.k,1));
            const double ey=-(u(c.i,c.j,c.k,2)*b(c.i,c.j,c.k,0)-
                              u(c.i,c.j,c.k,0)*b(c.i,c.j,c.k,2));
            const double ez=-(u(c.i,c.j,c.k,0)*b(c.i,c.j,c.k,1)-
                              u(c.i,c.j,c.k,1)*b(c.i,c.j,c.k,0));
            emf=ex*edge.canonical_edge_vector[0]+
                ey*edge.canonical_edge_vector[1]+
                ez*edge.canonical_edge_vector[2];
        }

        // Electron-pressure (ambipolar) contribution.  Along the singular
        // line only the edge-tangent derivative is needed; transverse
        // collapsed ghosts never enter this centered cell derivative.
        if (ambipolar_control.enabled)
        {
            FieldBlock &pH=fld_->field(fid_.fid_PV_H,c.block);
            FieldBlock &pN=fld_->field(fid_.fid_PV_Na,c.block);
            FieldBlock &uH=fld_->field(fid_.fid_U_H,c.block);
            FieldBlock &uN=fld_->field(fid_.fid_U_Na,c.block);
            int im=c.i,jm=c.j,km=c.k, ip=c.i,jp=c.j,kp=c.k;
            const int ax=static_cast<int>(inc.source_alias.axis);
            if (ax==0) { --im; ++ip; } else if (ax==1) { --jm; ++jp; } else { --km; ++kp; }
            const double dp=0.5*((pH(ip,jp,kp,3)+pN(ip,jp,kp,3))-
                                 (pH(im,jm,km,3)+pN(im,jm,km,3)));
            const NumInfo num=Hall_Num_Limiter(uH(c.i,c.j,c.k,0),uN(c.i,c.j,c.k,0));
            const double pressure_coef=(rho_ref*U_ref)/(q_e*L_ref*B_ref*n_ref);
            emf -= pressure_coef*inc.source_orientation*dp/num.ne_eff;
        }

        // Artificial resistivity is reconstructed from physical cell J and
        // the canonical edge vector, rather than a block-local edge alias.
        FieldBlock &jc=fld_->field(fid_.fid_Jcell,c.block);
        const double jedge=jc(c.i,c.j,c.k,0)*edge.canonical_edge_vector[0]+
                           jc(c.i,c.j,c.k,1)*edge.canonical_edge_vector[1]+
                           jc(c.i,c.j,c.k,2)*edge.canonical_edge_vector[2];
        const double jmag=std::sqrt(jc(c.i,c.j,c.k,0)*jc(c.i,c.j,c.k,0)+
                                    jc(c.i,c.j,c.k,1)*jc(c.i,c.j,c.k,1)+
                                    jc(c.i,c.j,c.k,2)*jc(c.i,c.j,c.k,2));
        if (arti_resist_control.eta_max>0.0)
        {
            double q=(jmag-arti_resist_control.J_range_start)/
                     std::max(arti_resist_control.J_range_on-arti_resist_control.J_range_start,1.e-30);
            q=std::max(0.0,std::min(1.0,q));
            emf += arti_resist_control.eta_max*q*q*(3.0-2.0*q)*jedge;
        }
        if (arti_resist_control.local_enabled && arti_resist_control.local_eta_max>0.0)
        {
            auto &cx=grd_->grids(c.block).dual_x; auto &cy=grd_->grids(c.block).dual_y; auto &cz=grd_->grids(c.block).dual_z;
            const double dx=cx(c.i+1,c.j+1,c.k+1)-arti_resist_control.local_center[0];
            const double dy=cy(c.i+1,c.j+1,c.k+1)-arti_resist_control.local_center[1];
            const double dz=cz(c.i+1,c.j+1,c.k+1)-arti_resist_control.local_center[2];
            const double rr=std::sqrt(dx*dx+dy*dy+dz*dz);
            if (arti_resist_control.local_r_cutoff<=0.0 || rr<arti_resist_control.local_r_cutoff)
                emf += arti_resist_control.local_eta_max*std::exp(-rr/arti_resist_control.local_r_decay)*jedge;
        }
        return emf;
    };

    singular_edges_->assemble_cell_field_to_local_owners(*fld_,"E_xi",contribution);
    singular_edges_->assemble_cell_field_to_local_owners(*fld_,"E_eta",contribution);
    singular_edges_->assemble_cell_field_to_local_owners(*fld_,"E_zeta",contribution);
}

//=========================================================================
// Face candidates to the unique shared edge EMF.
// Eface_* stores candidate components (e_xi,e_eta,e_zeta), not Cartesian E.
//
// This version adds a simple UCT transverse upwind correction.
// It keeps the CT structure because the final B update is still curl(E_edge).
void MercurySolver::AssembleEdgeEMF_FromFaceE_Ideal_()
{
    constexpr double jac_floor = 1.0e-30;

    // UCT correction strength.
    // Start with 0.5. If still too weak, try 1.0.
    // If it becomes noisy, reduce to 0.25.
    constexpr double Cuct = 0.5;

    auto shift_cell = [](int axis, int side, int &i, int &j, int &k)
    {
        if (axis == 0)
            i += side;
        else if (axis == 1)
            j += side;
        else
            k += side;
    };

    Int3 sub, sup;

    for (int iblk = 0; iblk < fld_->num_blocks(); iblk++)
    {
        auto &Uplus = fld_->field(fid_.fid_U_plus, iblk);

        auto &Bxi = fld_->field(fid_.fid_B.xi, iblk);
        auto &Beta = fld_->field(fid_.fid_B.eta, iblk);
        auto &Bzeta = fld_->field(fid_.fid_B.zeta, iblk);

        auto &Jac = fld_->field(fid_.fid_Jac, iblk);
        auto &JDxi = fld_->field(fid_.fid_metric.xi, iblk);
        auto &JDet = fld_->field(fid_.fid_metric.eta, iblk);
        auto &JDze = fld_->field(fid_.fid_metric.zeta, iblk);

        if (!Uplus.is_allocated())
            continue;

        auto metric_component = [&](int axis, int i, int j, int k, int comp) -> double
        {
            if (axis == 0)
                return JDxi(i, j, k, comp);
            if (axis == 1)
                return JDet(i, j, k, comp);
            return JDze(i, j, k, comp);
        };

        auto cell_u_contra = [&](int axis, int i, int j, int k) -> double
        {
            int ip = i, jp = j, kp = k;
            shift_cell(axis, 1, ip, jp, kp);

            const double kx = 0.5 * (metric_component(axis, i, j, k, 0) +
                                     metric_component(axis, ip, jp, kp, 0));

            const double ky = 0.5 * (metric_component(axis, i, j, k, 1) +
                                     metric_component(axis, ip, jp, kp, 1));

            const double kz = 0.5 * (metric_component(axis, i, j, k, 2) +
                                     metric_component(axis, ip, jp, kp, 2));

            const double inv_jac = 1.0 / (std::abs(Jac(i, j, k, 0)) + jac_floor);

            return (kx * Uplus(i, j, k, 0) +
                    ky * Uplus(i, j, k, 1) +
                    kz * Uplus(i, j, k, 2)) *
                   inv_jac;
        };

        auto max_abs_u4 = [&](int axis,
                              int i0, int j0, int k0,
                              int i1, int j1, int k1,
                              int i2, int j2, int k2,
                              int i3, int j3, int k3) -> double
        {
            double a = 0.0;
            a = std::max(a, std::abs(cell_u_contra(axis, i0, j0, k0)));
            a = std::max(a, std::abs(cell_u_contra(axis, i1, j1, k1)));
            a = std::max(a, std::abs(cell_u_contra(axis, i2, j2, k2)));
            a = std::max(a, std::abs(cell_u_contra(axis, i3, j3, k3)));
            return a;
        };

        //=============================================================
        // E_xi edge
        //
        // E_xi = central
        //      + 0.5 * a_eta  * d_beta_zeta / d_eta
        //      - 0.5 * a_zeta * d_beta_eta  / d_zeta
        //
        // Here beta means the evolved face 2-form B, not Badd.
        {
            auto &Exi = fld_->field(fid_.fid_E.xi, iblk);
            auto &E_face_eta = fld_->field(fid_.fid_Eface.eta, iblk);
            auto &E_face_zeta = fld_->field(fid_.fid_Eface.zeta, iblk);

            sub = Exi.inner_lo();
            sup = Exi.inner_hi();

            for (int i = sub.i; i < sup.i; i++)
                for (int j = sub.j; j < sup.j; j++)
                    for (int k = sub.k; k < sup.k; k++)
                    {
                        const double Ec = 0.25 * (E_face_eta(i, j, k, 0) +
                                                  E_face_eta(i, j, k - 1, 0) +
                                                  E_face_zeta(i, j, k, 0) +
                                                  E_face_zeta(i, j - 1, k, 0));

                        const double a_eta = max_abs_u4(
                            1,
                            i, j, k,
                            i, j - 1, k,
                            i, j, k - 1,
                            i, j - 1, k - 1);

                        const double a_zeta = max_abs_u4(
                            2,
                            i, j, k,
                            i, j - 1, k,
                            i, j, k - 1,
                            i, j - 1, k - 1);

                        const double dBzeta_eta =
                            Bzeta(i, j, k, 0) -
                            Bzeta(i, j - 1, k, 0);

                        const double dBeta_zeta =
                            Beta(i, j, k, 0) -
                            Beta(i, j, k - 1, 0);

                        Exi(i, j, k, 0) =
                            Ec + Cuct * (0.5 * a_eta * dBzeta_eta - 0.5 * a_zeta * dBeta_zeta);
                    }
        }

        //=============================================================
        // E_eta edge
        //
        // E_eta = central
        //       + 0.5 * a_zeta * d_beta_xi   / d_zeta
        //       - 0.5 * a_xi   * d_beta_zeta / d_xi
        {
            auto &Eeta = fld_->field(fid_.fid_E.eta, iblk);
            auto &E_face_xi = fld_->field(fid_.fid_Eface.xi, iblk);
            auto &E_face_zeta = fld_->field(fid_.fid_Eface.zeta, iblk);

            sub = Eeta.inner_lo();
            sup = Eeta.inner_hi();

            for (int i = sub.i; i < sup.i; i++)
                for (int j = sub.j; j < sup.j; j++)
                    for (int k = sub.k; k < sup.k; k++)
                    {
                        const double Ec = 0.25 * (E_face_xi(i, j, k, 1) +
                                                  E_face_xi(i, j, k - 1, 1) +
                                                  E_face_zeta(i, j, k, 1) +
                                                  E_face_zeta(i - 1, j, k, 1));

                        const double a_zeta = max_abs_u4(
                            2,
                            i, j, k,
                            i - 1, j, k,
                            i, j, k - 1,
                            i - 1, j, k - 1);

                        const double a_xi = max_abs_u4(
                            0,
                            i, j, k,
                            i - 1, j, k,
                            i, j, k - 1,
                            i - 1, j, k - 1);

                        const double dBxi_zeta =
                            Bxi(i, j, k, 0) -
                            Bxi(i, j, k - 1, 0);

                        const double dBzeta_xi =
                            Bzeta(i, j, k, 0) -
                            Bzeta(i - 1, j, k, 0);

                        Eeta(i, j, k, 0) =
                            Ec + Cuct * (0.5 * a_zeta * dBxi_zeta - 0.5 * a_xi * dBzeta_xi);
                    }
        }

        //=============================================================
        // E_zeta edge
        //
        // E_zeta = central
        //        + 0.5 * a_xi  * d_beta_eta / d_xi
        //        - 0.5 * a_eta * d_beta_xi  / d_eta
        {
            auto &Ezeta = fld_->field(fid_.fid_E.zeta, iblk);
            auto &E_face_xi = fld_->field(fid_.fid_Eface.xi, iblk);
            auto &E_face_eta = fld_->field(fid_.fid_Eface.eta, iblk);

            sub = Ezeta.inner_lo();
            sup = Ezeta.inner_hi();

            for (int i = sub.i; i < sup.i; i++)
                for (int j = sub.j; j < sup.j; j++)
                    for (int k = sub.k; k < sup.k; k++)
                    {
                        const double Ec = 0.25 * (E_face_xi(i, j, k, 2) +
                                                  E_face_xi(i, j - 1, k, 2) +
                                                  E_face_eta(i, j, k, 2) +
                                                  E_face_eta(i - 1, j, k, 2));

                        const double a_xi = max_abs_u4(
                            0,
                            i, j, k,
                            i - 1, j, k,
                            i, j - 1, k,
                            i - 1, j - 1, k);

                        const double a_eta = max_abs_u4(
                            1,
                            i, j, k,
                            i - 1, j, k,
                            i, j - 1, k,
                            i - 1, j - 1, k);

                        const double dBeta_xi =
                            Beta(i, j, k, 0) -
                            Beta(i - 1, j, k, 0);

                        const double dBxi_eta =
                            Bxi(i, j, k, 0) -
                            Bxi(i, j - 1, k, 0);

                        Ezeta(i, j, k, 0) =
                            Ec + Cuct * (0.5 * a_xi * dBeta_xi - 0.5 * a_eta * dBxi_eta);
                    }
        }
    }
}

//=========================================================================
// Face candidates to the unique shared edge EMF.
// Eface_* stores candidate components (e_xi,e_eta,e_zeta), not Cartesian E.
// void MercurySolver::AssembleEdgeEMF_FromFaceE_Ideal_()
// {
//     Int3 sub, sup;
//     for (int iblk = 0; iblk < fld_->num_blocks(); iblk++)
//     {
//         {
//             auto &Exi = fld_->field(fid_.fid_E.xi, iblk);
//             auto &E_face_eta = fld_->field(fid_.fid_Eface.eta, iblk);
//             auto &E_face_zeta = fld_->field(fid_.fid_Eface.zeta, iblk);
//             sub = Exi.inner_lo();
//             sup = Exi.inner_hi();
//             for (int i = sub.i; i < sup.i; i++)
//                 for (int j = sub.j; j < sup.j; j++)
//                     for (int k = sub.k; k < sup.k; k++)
//                     {
//                         Exi(i, j, k, 0) = 0.25 * (
//                             E_face_eta(i, j, k, 0) + E_face_eta(i, j, k - 1, 0) +
//                             E_face_zeta(i, j, k, 0) + E_face_zeta(i, j - 1, k, 0));
//                     }
//         }

//         {
//             auto &Eeta = fld_->field(fid_.fid_E.eta, iblk);
//             auto &E_face_xi = fld_->field(fid_.fid_Eface.xi, iblk);
//             auto &E_face_zeta = fld_->field(fid_.fid_Eface.zeta, iblk);
//             sub = Eeta.inner_lo();
//             sup = Eeta.inner_hi();
//             for (int i = sub.i; i < sup.i; i++)
//                 for (int j = sub.j; j < sup.j; j++)
//                     for (int k = sub.k; k < sup.k; k++)
//                     {
//                         Eeta(i, j, k, 0) = 0.25 * (
//                             E_face_xi(i, j, k, 1) + E_face_xi(i, j, k - 1, 1) +
//                             E_face_zeta(i, j, k, 1) + E_face_zeta(i - 1, j, k, 1));
//                     }
//         }

//         {
//             auto &Ezeta = fld_->field(fid_.fid_E.zeta, iblk);
//             auto &E_face_xi = fld_->field(fid_.fid_Eface.xi, iblk);
//             auto &E_face_eta = fld_->field(fid_.fid_Eface.eta, iblk);
//             sub = Ezeta.inner_lo();
//             sup = Ezeta.inner_hi();
//             for (int i = sub.i; i < sup.i; i++)
//                 for (int j = sub.j; j < sup.j; j++)
//                     for (int k = sub.k; k < sup.k; k++)
//                     {
//                         Ezeta(i, j, k, 0) = 0.25 * (
//                             E_face_xi(i, j, k, 2) + E_face_xi(i, j - 1, k, 2) +
//                             E_face_eta(i, j, k, 2) + E_face_eta(i - 1, j, k, 2));
//                     }
//         }
//     }
// }

//=========================================================================
void MercurySolver::AssembleOneDirectionEMF_(
    int dir, FieldBlock &E_face,
    FieldBlock &Bxi, FieldBlock &Beta, FieldBlock &Bzeta,
    FieldBlock &Badd_xi, FieldBlock &Badd_eta, FieldBlock &Badd_zeta,
    FieldBlock &Jac,
    FieldBlock &JDxi, FieldBlock &JDet, FieldBlock &JDze,
    FieldBlock &Uplus)
{
    constexpr double jac_floor = 1.0e-30;

    auto shift_cell = [](int axis, int side, int &i, int &j, int &k)
    {
        if (axis == 0)
            i += side;
        else if (axis == 1)
            j += side;
        else
            k += side;
    };

    auto face_beta_ind = [&](int axis, int i, int j, int k) -> double
    {
        if (axis == 0)
            return Bxi(i, j, k, 0);
        if (axis == 1)
            return Beta(i, j, k, 0);
        return Bzeta(i, j, k, 0);
    };

    auto face_beta_add = [&](int axis, int i, int j, int k) -> double
    {
        if (axis == 0)
            return Badd_xi(i, j, k, 0);
        if (axis == 1)
            return Badd_eta(i, j, k, 0);
        return Badd_zeta(i, j, k, 0);
    };

    auto face_beta_total = [&](int axis, int i, int j, int k) -> double
    {
        return face_beta_ind(axis, i, j, k) + face_beta_add(axis, i, j, k);
    };

    auto cell_beta_ind = [&](int axis, int i, int j, int k) -> double
    {
        int ip = i, jp = j, kp = k;
        shift_cell(axis, 1, ip, jp, kp);
        return 0.5 * (face_beta_ind(axis, i, j, k) + face_beta_ind(axis, ip, jp, kp));
    };

    auto cell_beta_total = [&](int axis, int i, int j, int k) -> double
    {
        int ip = i, jp = j, kp = k;
        shift_cell(axis, 1, ip, jp, kp);
        return 0.5 * (face_beta_total(axis, i, j, k) + face_beta_total(axis, ip, jp, kp));
    };

    auto metric_component = [&](int axis, int i, int j, int k, int comp) -> double
    {
        if (axis == 0)
            return JDxi(i, j, k, comp);
        if (axis == 1)
            return JDet(i, j, k, comp);
        return JDze(i, j, k, comp);
    };

    auto cell_u_contra = [&](int axis, int i, int j, int k) -> double
    {
        int ip = i, jp = j, kp = k;
        shift_cell(axis, 1, ip, jp, kp);

        const double kx = 0.5 * (metric_component(axis, i, j, k, 0) + metric_component(axis, ip, jp, kp, 0));
        const double ky = 0.5 * (metric_component(axis, i, j, k, 1) + metric_component(axis, ip, jp, kp, 1));
        const double kz = 0.5 * (metric_component(axis, i, j, k, 2) + metric_component(axis, ip, jp, kp, 2));
        const double inv_jac = 1.0 / (std::abs(Jac(i, j, k, 0)) + jac_floor);

        return (kx * Uplus(i, j, k, 0) +
                ky * Uplus(i, j, k, 1) +
                kz * Uplus(i, j, k, 2)) *
               inv_jac;
    };

    auto flux = [&](int beta_axis, int flux_axis,
                    int iL, int jL, int kL,
                    int iR, int jR, int kR,
                    double beta_flux_face) -> double
    {
        const double betaTotalL = cell_beta_total(beta_axis, iL, jL, kL);
        const double betaTotalR = cell_beta_total(beta_axis, iR, jR, kR);
        const double betaIndL = cell_beta_ind(beta_axis, iL, jL, kL);
        const double betaIndR = cell_beta_ind(beta_axis, iR, jR, kR);

        const double uFluxL = cell_u_contra(flux_axis, iL, jL, kL);
        const double uFluxR = cell_u_contra(flux_axis, iR, jR, kR);
        const double uBetaL = cell_u_contra(beta_axis, iL, jL, kL);
        const double uBetaR = cell_u_contra(beta_axis, iR, jR, kR);

        const double GL = uFluxL * betaTotalL - uBetaL * beta_flux_face;
        const double GR = uFluxR * betaTotalR - uBetaR * beta_flux_face;
        const double radius = std::max(std::abs(uFluxL), std::abs(uFluxR));

        return 0.5 * (GL + GR) - 0.5 * radius * (betaIndR - betaIndL);
    };

    const int flux_axis = dir - 1;
    if (flux_axis < 0 || flux_axis > 2)
        throw std::runtime_error("AssembleOneDirectionEMF_: invalid direction");

    Int3 sub = E_face.inner_lo();
    Int3 sup = E_face.inner_hi();

    for (int i = sub.i; i < sup.i; ++i)
        for (int j = sub.j; j < sup.j; ++j)
            for (int k = sub.k; k < sup.k; ++k)
            {
                int iL = i, jL = j, kL = k;
                int iR = i, jR = j, kR = k;
                shift_cell(flux_axis, -1, iL, jL, kL);

                const double beta_flux_face = face_beta_total(flux_axis, i, j, k);
                double e[3] = {0.0, 0.0, 0.0};

                if (flux_axis == 0)
                {
                    const double G_eta_xi = flux(1, 0, iL, jL, kL, iR, jR, kR, beta_flux_face);
                    const double G_zeta_xi = flux(2, 0, iL, jL, kL, iR, jR, kR, beta_flux_face);
                    e[1] = G_zeta_xi;
                    e[2] = -G_eta_xi;
                }
                else if (flux_axis == 1)
                {
                    const double G_xi_eta = flux(0, 1, iL, jL, kL, iR, jR, kR, beta_flux_face);
                    const double G_zeta_eta = flux(2, 1, iL, jL, kL, iR, jR, kR, beta_flux_face);
                    e[0] = -G_zeta_eta;
                    e[2] = G_xi_eta;
                }
                else
                {
                    const double G_xi_zeta = flux(0, 2, iL, jL, kL, iR, jR, kR, beta_flux_face);
                    const double G_eta_zeta = flux(1, 2, iL, jL, kL, iR, jR, kR, beta_flux_face);
                    e[0] = G_eta_zeta;
                    e[1] = -G_xi_zeta;
                }

                E_face(i, j, k, 0) = e[0];
                E_face(i, j, k, 1) = e[1];
                E_face(i, j, k, 2) = e[2];
            }
}
