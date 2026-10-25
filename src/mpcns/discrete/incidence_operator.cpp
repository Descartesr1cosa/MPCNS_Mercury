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
void apply_d1_impl(const PrimalEdgeFormView &edge,
                   PrimalFaceFormView face,
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

void apply_d1(const PrimalEdgeFormView &edge,
              PrimalFaceFormView face,
              double scale)
{
    apply_d1_impl<false>(edge, face, scale);
}

void subtract_d1(const PrimalEdgeFormView &edge,
                 PrimalFaceFormView face,
                 double scale)
{
    apply_d1_impl<true>(edge, face, scale);
}

} // namespace mpcns::discrete
