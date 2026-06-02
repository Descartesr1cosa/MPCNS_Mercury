#include "MercurySolver.h"
#include <algorithm>
#include <vector>

void MercurySolver::AssembleRHS_Induction_CT_()
{
    const int nb = fld_->num_blocks();

    // 0) E_edge 清零
    for (int ib = 0; ib < nb; ++ib)
    {
        auto &Exi = fld_->field(fid_.fid_E.xi, ib);
        auto &Eeta = fld_->field(fid_.fid_E.eta, ib);
        auto &Eze = fld_->field(fid_.fid_E.zeta, ib);
        if (!Exi.is_allocated())
            continue;

        auto zero_electric = [&](FieldBlock &E)
        {
            Int3 lo = E.get_lo(), hi = E.get_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        E(i, j, k, 0) = 0.0;
        };
        zero_electric(Exi);
        zero_electric(Eeta);
        zero_electric(Eze);
    }

    AddIdealEdgeEMF_();
    AddAmbipolarEdgeEMF_();
    AddArtificialResistivityToEdgeEMF_();

    mercury_bound_.Sync("Eedge");

    for (int ib = 0; ib < nb; ++ib)
    {
        auto &Exi = fld_->field(fid_.fid_E.xi, ib);
        auto &Eeta = fld_->field(fid_.fid_E.eta, ib);
        auto &Eze = fld_->field(fid_.fid_E.zeta, ib);

        auto &RHSBxi = fld_->field(fid_.fid_RHS_b.xi, ib);
        auto &RHSBeta = fld_->field(fid_.fid_RHS_b.eta, ib);
        auto &RHSBze = fld_->field(fid_.fid_RHS_b.zeta, ib);
        if (!Exi.is_allocated())
            continue;

        CTOperators::CurlEdgeToFace(ib, Exi, Eeta, Eze, RHSBxi, RHSBeta, RHSBze, /*multiper=*/-1.0);
    }
}

void MercurySolver::ResistiveDiffusionSubcycles_()
{
    if (!resist_control.is_Mercury_resistance)
        return;

    const int nsub = std::max(1, resist_control.n_subcycles);
    const double dt_res = dt / static_cast<double>(nsub);
    const int nb = fld_->num_blocks();

    auto zero_one = [](FieldBlock &F)
    {
        if (!F.is_allocated())
            return;

        Int3 lo = F.get_lo();
        Int3 hi = F.get_hi();
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                    for (int m = 0; m < F.descriptor().ncomp; ++m)
                        F(i, j, k, m) = 0.0;
    };

    for (int isub = 0; isub < nsub; ++isub)
    {
        mercury_bound_.Sync("Bface");
        UpdateMagneticDerivedFields_();

        for (int ib = 0; ib < nb; ++ib)
        {
            zero_one(fld_->field(fid_.fid_Eres.xi, ib));
            zero_one(fld_->field(fid_.fid_Eres.eta, ib));
            zero_one(fld_->field(fid_.fid_Eres.zeta, ib));
            zero_one(fld_->field(fid_.fid_RHS_b_res.xi, ib));
            zero_one(fld_->field(fid_.fid_RHS_b_res.eta, ib));
            zero_one(fld_->field(fid_.fid_RHS_b_res.zeta, ib));
        }

        AddResistiveEdgeEMF_To_(fid_.fid_Eres);
        mercury_bound_.Sync("Eres");

        for (int ib = 0; ib < nb; ++ib)
        {
            auto &Exi = fld_->field(fid_.fid_Eres.xi, ib);
            auto &Eeta = fld_->field(fid_.fid_Eres.eta, ib);
            auto &Eze = fld_->field(fid_.fid_Eres.zeta, ib);

            auto &RHSBxi = fld_->field(fid_.fid_RHS_b_res.xi, ib);
            auto &RHSBeta = fld_->field(fid_.fid_RHS_b_res.eta, ib);
            auto &RHSBze = fld_->field(fid_.fid_RHS_b_res.zeta, ib);

            if (!Exi.is_allocated())
                continue;

            CTOperators::CurlEdgeToFace(ib, Exi, Eeta, Eze,
                                        RHSBxi, RHSBeta, RHSBze,
                                        /*multiper=*/-1.0);
        }

        ApplyUpdate_Euler_BfaceOnly_(dt_res, fid_.fid_RHS_b_res);
        mercury_bound_.Sync("Bface");
    }
}
