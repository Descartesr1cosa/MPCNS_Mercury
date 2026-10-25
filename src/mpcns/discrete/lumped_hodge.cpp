#include "mpcns/discrete/lumped_hodge.hpp"

namespace mpcns::discrete
{

LumpedHodge::LumpedHodge(LumpedHodgeCoefficients coefficients)
    : coefficients_(coefficients)
{
}

void LumpedHodge::apply_star_2(const state::PrimalFaceFormView &primal_face,
                               state::DualEdgeFormView dual_edge) const
{
    {
        const Int3 lo = dual_edge.xi.inner_lo();
        const Int3 hi = dual_edge.xi.inner_hi();
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                    dual_edge.xi(i, j, k, dual_edge.component) =
                        coefficients_.star_primal_2_to_dual_1.xi(
                            i, j, k, coefficients_.star_primal_2_to_dual_1.component) *
                        primal_face.xi(i, j, k, primal_face.component);
    }

    {
        const Int3 lo = dual_edge.eta.inner_lo();
        const Int3 hi = dual_edge.eta.inner_hi();
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                    dual_edge.eta(i, j, k, dual_edge.component) =
                        coefficients_.star_primal_2_to_dual_1.eta(
                            i, j, k, coefficients_.star_primal_2_to_dual_1.component) *
                        primal_face.eta(i, j, k, primal_face.component);
    }

    {
        const Int3 lo = dual_edge.zeta.inner_lo();
        const Int3 hi = dual_edge.zeta.inner_hi();
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                    dual_edge.zeta(i, j, k, dual_edge.component) =
                        coefficients_.star_primal_2_to_dual_1.zeta(
                            i, j, k, coefficients_.star_primal_2_to_dual_1.component) *
                        primal_face.zeta(i, j, k, primal_face.component);
    }
}

void LumpedHodge::apply_inverse_star_1(const state::DualFaceFormView &dual_face,
                                       state::PrimalEdgeFormView primal_edge) const
{
    {
        const Int3 lo = primal_edge.xi.inner_lo();
        const Int3 hi = primal_edge.xi.inner_hi();
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                    primal_edge.xi(i, j, k, primal_edge.component) =
                        coefficients_.inverse_star_primal_1.xi(
                            i, j, k, coefficients_.inverse_star_primal_1.component) *
                        dual_face.xi(i, j, k, dual_face.component);
    }

    {
        const Int3 lo = primal_edge.eta.inner_lo();
        const Int3 hi = primal_edge.eta.inner_hi();
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                    primal_edge.eta(i, j, k, primal_edge.component) =
                        coefficients_.inverse_star_primal_1.eta(
                            i, j, k, coefficients_.inverse_star_primal_1.component) *
                        dual_face.eta(i, j, k, dual_face.component);
    }

    {
        const Int3 lo = primal_edge.zeta.inner_lo();
        const Int3 hi = primal_edge.zeta.inner_hi();
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                    primal_edge.zeta(i, j, k, primal_edge.component) =
                        coefficients_.inverse_star_primal_1.zeta(
                            i, j, k, coefficients_.inverse_star_primal_1.component) *
                        dual_face.zeta(i, j, k, dual_face.component);
    }
}

} // namespace mpcns::discrete
