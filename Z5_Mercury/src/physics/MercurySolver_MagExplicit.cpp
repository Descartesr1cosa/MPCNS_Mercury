
#include "MercurySolver.h"

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

        auto &Hodge_star_2form_to_1form_face_xi = fld_->field(fid_.Hodge_star_2form_to_1form_face.xi, iblk);
        auto &Hodge_star_2form_to_1form_face_eta = fld_->field(fid_.Hodge_star_2form_to_1form_face.eta, iblk);
        auto &Hodge_star_2form_to_1form_face_zeta = fld_->field(fid_.Hodge_star_2form_to_1form_face.zeta, iblk);

        auto &Hodge_star_inverse_2form_to_1form_edge_xi = fld_->field(fid_.Hodge_star_inverse_2form_to_1form_edge.xi, iblk);
        auto &Hodge_star_inverse_2form_to_1form_edge_eta = fld_->field(fid_.Hodge_star_inverse_2form_to_1form_edge.eta, iblk);
        auto &Hodge_star_inverse_2form_to_1form_edge_zeta = fld_->field(fid_.Hodge_star_inverse_2form_to_1form_edge.zeta, iblk);

        // compute J (edge 1-form) from face B (2-form)
        // multiper 用 +1.0 J =curl B。
        CTOperators::CurlAdjFaceToEdge(iblk,
                                       Bxi, Beta, Bzeta,
                                       Hodge_star_2form_to_1form_face_xi, Hodge_star_2form_to_1form_face_eta, Hodge_star_2form_to_1form_face_zeta,
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
                        Jxi(i, j, k, 0) *= (Hodge_star_inverse_2form_to_1form_edge_xi(i, j, k, 0));
        }
        {
            // Edge_eta
            Int3 lo = Jeta.inner_lo(), hi = Jeta.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        Jeta(i, j, k, 0) *= (Hodge_star_inverse_2form_to_1form_edge_eta(i, j, k, 0));
        }
        {
            // Edge_zeta
            Int3 lo = Jzeta.inner_lo(), hi = Jzeta.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        Jzeta(i, j, k, 0) *= (Hodge_star_inverse_2form_to_1form_edge_zeta(i, j, k, 0));
        }
    }

    mercury_bound_.Sync("Jedge"); // ApplyBC_EdgeJ_();
}
