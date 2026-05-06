
#include "MercurySolver.h"

void MercurySolver::Build_E_explicit_edge_()
{
    AddIdealEdgeEMF_();

    // AddResistiveEdgeEMF_(); // Add magnetic diffusion in solid (and optionally fluid)
    AddPoleResistiveEdgeEMF_FromJcell_();

    // AddAmbipolarEdgeEMF_();

    // 后续CT只会用到inner的电场，不会使用nghost区域，因此只需对Pole处理即可
    // bound_.add_Edge_pole_boundary("E_xi"); // pole边界处理
    // bound_.add_Edge_pole_boundary("E_eta");
}

void MercurySolver::Calc_J_Edge()
{
    //  ComputeJ_AtEdges_Inner_();
    for (int iblk = 0; iblk < fld_->num_blocks(); ++iblk)
    {
        auto &Bxi = fld_->field(fid_.fid_B.xi, iblk);
        auto &Beta = fld_->field(fid_.fid_B.eta, iblk);
        auto &Bzeta = fld_->field(fid_.fid_B.zeta, iblk);

        auto &Jxi = fld_->field(fid_.fid_J.xi, iblk);
        auto &Jeta = fld_->field(fid_.fid_J.eta, iblk);
        auto &Jzeta = fld_->field(fid_.fid_J.zeta, iblk);

        auto &beta_xi = fld_->field(fid_.Face_beta.xi, iblk); // Hodge *: 2-form face -> 1-form face
        auto &beta_eta = fld_->field(fid_.Face_beta.eta, iblk);
        auto &beta_zeta = fld_->field(fid_.Face_beta.zeta, iblk);

        auto &alpha_xi = fld_->field(fid_.Edge_alpha.xi, iblk); // Hodge *: 2-form edge -> 1-form edge
        auto &alpha_eta = fld_->field(fid_.Edge_alpha.eta, iblk);
        auto &alpha_zeta = fld_->field(fid_.Edge_alpha.zeta, iblk);

        // compute J (edge 1-form) from face B (2-form)
        // multiper 用 +1.0 J =curl B。
        CTOperators::CurlAdjFaceToEdge(iblk,
                                       Bxi, Beta, Bzeta,
                                       beta_xi, beta_eta, beta_zeta,
                                       Jxi, Jeta, Jzeta,
                                       /*multiper=*/1.0);

        //  J_edge^(1-form) = (alpha_edge / mu0) * Jcirc_edge
        //     alpha = |e|/|S*|  (⋆1^{-1})
        {
            // Edge_xi
            Int3 lo = Jxi.inner_lo(), hi = Jxi.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        Jxi(i, j, k, 0) *= (alpha_xi(i, j, k, 0));
        }
        {
            // Edge_eta
            Int3 lo = Jeta.inner_lo(), hi = Jeta.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        Jeta(i, j, k, 0) *= (alpha_eta(i, j, k, 0));
        }
        {
            // Edge_zeta
            Int3 lo = Jzeta.inner_lo(), hi = Jzeta.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        Jzeta(i, j, k, 0) *= (alpha_zeta(i, j, k, 0));
        }
    }

    mercury_bound_.Sync("Jedge"); // ApplyBC_EdgeJ_();
}
