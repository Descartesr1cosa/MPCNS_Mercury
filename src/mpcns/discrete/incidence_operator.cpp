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

struct Bounds { int ilo, jlo, klo, ihi, jhi, khi; };

Bounds inner_bounds(const FieldBlock &view)
{
    const Int3 lo = view.inner_lo();
    const Int3 hi = view.inner_hi();
    return {lo.i, lo.j, lo.k, hi.i, hi.j, hi.k};
}

Bounds inner_bounds(const state::StateView &view)
{
    const auto &lo = view.inner_extent().lo;
    const auto &hi = view.inner_extent().hi;
    return {lo.i, lo.j, lo.k, hi.i, hi.j, hi.k};
}

template <bool kSubtract, class EdgeForm, class FaceForm>
void apply_d1_impl(const EdgeForm &edge, FaceForm face, double scale)
{
    {
        const Bounds b = inner_bounds(face.xi);
        for (int i = b.ilo; i < b.ihi; ++i)
            for (int j = b.jlo; j < b.jhi; ++j)
                for (int k = b.klo; k < b.khi; ++k)
                    write_scaled<kSubtract>(
                        face.xi(i, j, k, face.component), scale,
                        (edge.eta(i, j, k, edge.component) - edge.eta(i, j, k + 1, edge.component)) +
                        (edge.zeta(i, j + 1, k, edge.component) - edge.zeta(i, j, k, edge.component)));
    }
    {
        const Bounds b = inner_bounds(face.eta);
        for (int i = b.ilo; i < b.ihi; ++i)
            for (int j = b.jlo; j < b.jhi; ++j)
                for (int k = b.klo; k < b.khi; ++k)
                    write_scaled<kSubtract>(
                        face.eta(i, j, k, face.component), scale,
                        (edge.xi(i, j, k + 1, edge.component) - edge.xi(i, j, k, edge.component)) +
                        (edge.zeta(i, j, k, edge.component) - edge.zeta(i + 1, j, k, edge.component)));
    }
    {
        const Bounds b = inner_bounds(face.zeta);
        for (int i = b.ilo; i < b.ihi; ++i)
            for (int j = b.jlo; j < b.jhi; ++j)
                for (int k = b.klo; k < b.khi; ++k)
                    write_scaled<kSubtract>(
                        face.zeta(i, j, k, face.component), scale,
                        (edge.xi(i, j, k, edge.component) - edge.xi(i, j + 1, k, edge.component)) +
                        (edge.eta(i + 1, j, k, edge.component) - edge.eta(i, j, k, edge.component)));
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

void apply_d1(const state::PrimalEdgeFormStateView &edge,
              state::PrimalFaceFormStateView face,
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

template <class DualEdgeForm, class DualFaceForm>
void apply_d1_transpose_impl(const DualEdgeForm &dual_edge,
                             DualFaceForm dual_face,
                             double scale)
{
    {
        const Bounds b = inner_bounds(dual_face.xi);
        for (int i = b.ilo; i < b.ihi; ++i)
            for (int j = b.jlo; j < b.jhi; ++j)
                for (int k = b.klo; k < b.khi; ++k)
                    dual_face.xi(i, j, k, dual_face.component) = scale *
                        ((dual_edge.eta(i, j, k - 1, dual_edge.component) -
                          dual_edge.eta(i, j, k, dual_edge.component)) +
                         (dual_edge.zeta(i, j, k, dual_edge.component) -
                          dual_edge.zeta(i, j - 1, k, dual_edge.component)));
    }
    {
        const Bounds b = inner_bounds(dual_face.eta);
        for (int i = b.ilo; i < b.ihi; ++i)
            for (int j = b.jlo; j < b.jhi; ++j)
                for (int k = b.klo; k < b.khi; ++k)
                    dual_face.eta(i, j, k, dual_face.component) = scale *
                        ((dual_edge.xi(i, j, k, dual_edge.component) -
                          dual_edge.xi(i, j, k - 1, dual_edge.component)) +
                         (dual_edge.zeta(i - 1, j, k, dual_edge.component) -
                          dual_edge.zeta(i, j, k, dual_edge.component)));
    }
    {
        const Bounds b = inner_bounds(dual_face.zeta);
        for (int i = b.ilo; i < b.ihi; ++i)
            for (int j = b.jlo; j < b.jhi; ++j)
                for (int k = b.klo; k < b.khi; ++k)
                    dual_face.zeta(i, j, k, dual_face.component) = scale *
                        ((dual_edge.xi(i, j - 1, k, dual_edge.component) -
                          dual_edge.xi(i, j, k, dual_edge.component)) +
                         (dual_edge.eta(i, j, k, dual_edge.component) -
                          dual_edge.eta(i - 1, j, k, dual_edge.component)));
    }
}

void apply_d1_transpose(const state::DualEdgeFormView &dual_edge,
                        state::DualFaceFormView dual_face,
                        double scale)
{
    apply_d1_transpose_impl(dual_edge, dual_face, scale);
}

void apply_d1_transpose(const state::DualEdgeFormStateView &dual_edge,
                        state::DualFaceFormView dual_face,
                        double scale)
{
    apply_d1_transpose_impl(dual_edge, dual_face, scale);
}

void apply_d1_transpose(const state::DualEdgeFormStateView &dual_edge,
                        state::DualFaceFormStateView dual_face,
                        double scale)
{
    apply_d1_transpose_impl(dual_edge, dual_face, scale);
}

} // namespace mpcns::discrete
