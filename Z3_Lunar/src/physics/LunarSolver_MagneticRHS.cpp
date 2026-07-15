#include "LunarSolver.h"
#include "0_basic/LayoutTraits.h"
#include <algorithm>

void LunarSolver::AssembleRHS_Induction_CT_()
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
    AddLocalArtificialResistivityToEdgeEMF_();

    // A regular physical edge on a rotated block interface has one complete
    // block-local UCT candidate per alias.  Those candidates use different
    // computational bases and are not bitwise identical even when U, B and
    // the interface geometry are conformal.  Owner-only synchronization would
    // arbitrarily keep one panel's candidate and creates an antisymmetric seam
    // forcing.  Reconcile all oriented alias candidates first, just as for
    // Jedge.  Singular physical edges are overwritten immediately below from
    // their variable incident-entity assembly, so they do not use this regular
    // edge average.
    ReduceEdgeAliasCandidatesToOwners_(fid_.fid_E);
    AssembleSingularEdgeEMF_NonHall_();

    AddHallEdgeEMF_();
    lunar_bound_.Sync("Ehall");
    for (int ib=0; ib<nb; ++ib)
    {
        auto add_one=[](FieldBlock &e,FieldBlock &h)
        {
            if (!e.is_allocated() || !h.is_allocated()) return;
            const Int3 lo=e.inner_lo(), hi=e.inner_hi();
            for(int i=lo.i;i<hi.i;++i) for(int j=lo.j;j<hi.j;++j) for(int k=lo.k;k<hi.k;++k)
                e(i,j,k,0)+=h(i,j,k,0);
        };
        add_one(fld_->field(fid_.fid_E.xi,ib),fld_->field(fid_.fid_Ehall.xi,ib));
        add_one(fld_->field(fid_.fid_E.eta,ib),fld_->field(fid_.fid_Ehall.eta,ib));
        add_one(fld_->field(fid_.fid_E.zeta,ib),fld_->field(fid_.fid_Ehall.zeta,ib));
    }

    // Apply the material-interface condition after every non-Lunar-
    // resistive contribution (ideal, ambipolar, artificial and Hall) has been
    // assembled.  Applying it immediately after the ideal term would allow
    // later terms to leak back onto the shared tangential edge.

    lunar_bound_.Sync("Eedge");

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
