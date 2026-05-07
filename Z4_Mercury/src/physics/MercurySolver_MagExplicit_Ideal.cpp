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

//=========================================================================
// Face candidates to the unique shared edge EMF.
// Eface_* stores candidate components (e_xi,e_eta,e_zeta), not Cartesian E.
void MercurySolver::AssembleEdgeEMF_FromFaceE_Ideal_()
{
    Int3 sub, sup;
    for (int iblk = 0; iblk < fld_->num_blocks(); iblk++)
    {
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
                        Exi(i, j, k, 0) = 0.25 * (
                            E_face_eta(i, j, k, 0) + E_face_eta(i, j, k - 1, 0) +
                            E_face_zeta(i, j, k, 0) + E_face_zeta(i, j - 1, k, 0));
                    }
        }

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
                        Eeta(i, j, k, 0) = 0.25 * (
                            E_face_xi(i, j, k, 1) + E_face_xi(i, j, k - 1, 1) +
                            E_face_zeta(i, j, k, 1) + E_face_zeta(i - 1, j, k, 1));
                    }
        }

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
                        Ezeta(i, j, k, 0) = 0.25 * (
                            E_face_xi(i, j, k, 2) + E_face_xi(i, j - 1, k, 2) +
                            E_face_eta(i, j, k, 2) + E_face_eta(i - 1, j, k, 2));
                    }
        }
    }
}

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
