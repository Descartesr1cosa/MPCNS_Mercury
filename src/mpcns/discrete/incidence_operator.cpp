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

void apply_d0(const state::PrimalZeroFormView &node,
              state::PrimalEdgeFormView edge)
{
    {
        const Int3 lo = edge.xi.inner_lo();
        const Int3 hi = edge.xi.inner_hi();
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                    edge.xi(i, j, k, edge.component) =
                        node.node(i + 1, j, k, node.component) -
                        node.node(i, j, k, node.component);
    }
    {
        const Int3 lo = edge.eta.inner_lo();
        const Int3 hi = edge.eta.inner_hi();
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                    edge.eta(i, j, k, edge.component) =
                        node.node(i, j + 1, k, node.component) -
                        node.node(i, j, k, node.component);
    }
    {
        const Int3 lo = edge.zeta.inner_lo();
        const Int3 hi = edge.zeta.inner_hi();
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                    edge.zeta(i, j, k, edge.component) =
                        node.node(i, j, k + 1, node.component) -
                        node.node(i, j, k, node.component);
    }
}

void apply_d1(const state::PrimalEdgeFormView &edge,
              state::PrimalFaceFormView face,
              double scale)
{
    apply_d1_impl<false>(edge, face, scale);
}

void apply_d2(const state::PrimalFaceFormView &face,
              state::PrimalThreeFormView cell)
{
    const Int3 lo = cell.cell.inner_lo();
    const Int3 hi = cell.cell.inner_hi();
    for (int i = lo.i; i < hi.i; ++i)
        for (int j = lo.j; j < hi.j; ++j)
            for (int k = lo.k; k < hi.k; ++k)
                cell.cell(i, j, k, cell.component) =
                    (face.xi(i + 1, j, k, face.component) -
                     face.xi(i, j, k, face.component)) +
                    (face.eta(i, j + 1, k, face.component) -
                     face.eta(i, j, k, face.component)) +
                    (face.zeta(i, j, k + 1, face.component) -
                     face.zeta(i, j, k, face.component));
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
                    dual_face.xi(i, j, k, dual_face.component) = scale *
                        ((dual_edge.eta(i, j, k - 1, dual_edge.component) -
                          dual_edge.eta(i, j, k, dual_edge.component)) +
                         (dual_edge.zeta(i, j, k, dual_edge.component) -
                          dual_edge.zeta(i, j - 1, k, dual_edge.component)));
    }
    {
        const Int3 lo = dual_face.eta.inner_lo();
        const Int3 hi = dual_face.eta.inner_hi();
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                    dual_face.eta(i, j, k, dual_face.component) = scale *
                        ((dual_edge.xi(i, j, k, dual_edge.component) -
                          dual_edge.xi(i, j, k - 1, dual_edge.component)) +
                         (dual_edge.zeta(i - 1, j, k, dual_edge.component) -
                          dual_edge.zeta(i, j, k, dual_edge.component)));
    }
    {
        const Int3 lo = dual_face.zeta.inner_lo();
        const Int3 hi = dual_face.zeta.inner_hi();
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                    dual_face.zeta(i, j, k, dual_face.component) = scale *
                        ((dual_edge.xi(i, j - 1, k, dual_edge.component) -
                          dual_edge.xi(i, j, k, dual_edge.component)) +
                         (dual_edge.eta(i, j, k, dual_edge.component) -
                          dual_edge.eta(i - 1, j, k, dual_edge.component)));
    }
}

} // namespace mpcns::discrete
