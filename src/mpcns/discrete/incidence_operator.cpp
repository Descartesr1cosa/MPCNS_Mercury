#include "mpcns/discrete/incidence_operator.hpp"

namespace mpcns::discrete
{

namespace
{

template <bool kSubtract>
inline void write_scaled(double &destination, double scale, double value)
{
    if constexpr (kSubtract)
        destination -= scale * value;
    else
        destination = scale * value;
}

template <bool kSubtract>
void apply_d1_impl(const state::PrimalEdgeFormView &edge,
                   state::PrimalFaceFormView face,
                   double scale)
{
    {
        const Int3 lo = face.xi.inner_lo();
        const Int3 hi = face.xi.inner_hi();
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                    write_scaled<kSubtract>(
                        face.xi(i, j, k, face.component), scale,
                        detail::d1_xi_at(edge, i, j, k));
    }

    {
        const Int3 lo = face.eta.inner_lo();
        const Int3 hi = face.eta.inner_hi();
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                    write_scaled<kSubtract>(
                        face.eta(i, j, k, face.component), scale,
                        detail::d1_eta_at(edge, i, j, k));
    }

    {
        const Int3 lo = face.zeta.inner_lo();
        const Int3 hi = face.zeta.inner_hi();
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                    write_scaled<kSubtract>(
                        face.zeta(i, j, k, face.component), scale,
                        detail::d1_zeta_at(edge, i, j, k));
    }
}

} // namespace

void apply_d1(const state::PrimalEdgeFormView &edge,
              state::PrimalFaceFormView face,
              double scale)
{
    apply_d1_impl<false>(edge, face, scale);
}

void subtract_d1(const state::PrimalEdgeFormView &edge,
                 state::PrimalFaceFormView face,
                 double scale)
{
    apply_d1_impl<true>(edge, face, scale);
}

void apply_d1_transpose(const state::DualEdgeFormView &dual_edge,
                        state::DualFaceFormView dual_face,
                        double scale)
{
    {
        const Int3 lo = dual_face.xi.inner_lo();
        const Int3 hi = dual_face.xi.inner_hi();
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                {
                    const double value =
                        (dual_edge.eta(i, j, k - 1, dual_edge.component) -
                         dual_edge.eta(i, j, k, dual_edge.component)) +
                        (dual_edge.zeta(i, j, k, dual_edge.component) -
                         dual_edge.zeta(i, j - 1, k, dual_edge.component));
                    dual_face.xi(i, j, k, dual_face.component) = scale * value;
                }
    }

    {
        const Int3 lo = dual_face.eta.inner_lo();
        const Int3 hi = dual_face.eta.inner_hi();
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                {
                    const double value =
                        (dual_edge.xi(i, j, k, dual_edge.component) -
                         dual_edge.xi(i, j, k - 1, dual_edge.component)) +
                        (dual_edge.zeta(i - 1, j, k, dual_edge.component) -
                         dual_edge.zeta(i, j, k, dual_edge.component));
                    dual_face.eta(i, j, k, dual_face.component) = scale * value;
                }
    }

    {
        const Int3 lo = dual_face.zeta.inner_lo();
        const Int3 hi = dual_face.zeta.inner_hi();
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                {
                    const double value =
                        (dual_edge.xi(i, j - 1, k, dual_edge.component) -
                         dual_edge.xi(i, j, k, dual_edge.component)) +
                        (dual_edge.eta(i, j, k, dual_edge.component) -
                         dual_edge.eta(i - 1, j, k, dual_edge.component));
                    dual_face.zeta(i, j, k, dual_face.component) = scale * value;
                }
    }
}

} // namespace mpcns::discrete
