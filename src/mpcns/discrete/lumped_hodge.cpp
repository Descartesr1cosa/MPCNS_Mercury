#include "mpcns/discrete/lumped_hodge.hpp"

#include <algorithm>
#include <stdexcept>

namespace mpcns::discrete
{
namespace
{

struct Bounds
{
    int ilo, jlo, klo, ihi, jhi, khi;
};

Bounds full_bounds(const FieldBlock &view)
{
    const Int3 lo = view.get_lo();
    const Int3 hi = view.get_hi();
    return {lo.i, lo.j, lo.k, hi.i, hi.j, hi.k};
}

Bounds full_bounds(const state::StateView &view)
{
    const auto &lo = view.extent().lo;
    const auto &hi = view.extent().hi;
    return {lo.i, lo.j, lo.k, hi.i, hi.j, hi.k};
}

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

template <class Source, class Target>
void apply_snapshot_star_2_component(const Source &source, Target &&target,
                                     const geometry::ConstFieldView &coefficient,
                                     int source_component, int target_component)
{
    // Geometry coefficients carry one ghost layer fewer than the state fields
    // (see LegacyGeometryFieldRegistration::ghost_layers), so clip to what the
    // coefficient view actually covers.  The dual-edge consumer only reads one
    // layer beyond the inner range.
    Bounds b = full_bounds(target);
    const auto &ce = coefficient.extent();
    b.ilo = std::max(b.ilo, ce.lo.i);
    b.jlo = std::max(b.jlo, ce.lo.j);
    b.klo = std::max(b.klo, ce.lo.k);
    b.ihi = std::min(b.ihi, ce.hi.i);
    b.jhi = std::min(b.jhi, ce.hi.j);
    b.khi = std::min(b.khi, ce.hi.k);
    for (int i = b.ilo; i < b.ihi; ++i)
        for (int j = b.jlo; j < b.jhi; ++j)
            for (int k = b.klo; k < b.khi; ++k)
            {
                const core::LogicalIndex index{i, j, k};
                target(i, j, k, target_component) =
                    coefficient.value(index) * source(i, j, k, source_component);
            }
}

template <class PrimalFace, class DualEdge>
void apply_snapshot_star_2(const PrimalFace &primal_face, DualEdge dual_edge,
                           const geometry::LumpedHodgeGeometry &coefficients)
{
    apply_snapshot_star_2_component(primal_face.xi, dual_edge.xi,
                                    coefficients.star_primal_2_to_dual_1.xi,
                                    primal_face.component, dual_edge.component);
    apply_snapshot_star_2_component(primal_face.eta, dual_edge.eta,
                                    coefficients.star_primal_2_to_dual_1.eta,
                                    primal_face.component, dual_edge.component);
    apply_snapshot_star_2_component(primal_face.zeta, dual_edge.zeta,
                                    coefficients.star_primal_2_to_dual_1.zeta,
                                    primal_face.component, dual_edge.component);
}

template <class Source, class Target>
void apply_snapshot_inverse_star_1_component(
    const Source &source, Target &&target,
    const geometry::ConstFieldView &coefficient,
    int source_component, int target_component)
{
    const Bounds b = inner_bounds(target);
    for (int i = b.ilo; i < b.ihi; ++i)
        for (int j = b.jlo; j < b.jhi; ++j)
            for (int k = b.klo; k < b.khi; ++k)
            {
                const core::LogicalIndex index{i, j, k};
                target(i, j, k, target_component) =
                    source(i, j, k, source_component) * coefficient.value(index);
            }
}

template <class DualFace, class PrimalEdge>
void apply_snapshot_inverse_star_1(
    const DualFace &dual_face, PrimalEdge primal_edge,
    const geometry::LumpedHodgeGeometry &coefficients)
{
    apply_snapshot_inverse_star_1_component(
        dual_face.xi, primal_edge.xi, coefficients.inverse_star_primal_1.xi,
        dual_face.component, primal_edge.component);
    apply_snapshot_inverse_star_1_component(
        dual_face.eta, primal_edge.eta, coefficients.inverse_star_primal_1.eta,
        dual_face.component, primal_edge.component);
    apply_snapshot_inverse_star_1_component(
        dual_face.zeta, primal_edge.zeta, coefficients.inverse_star_primal_1.zeta,
        dual_face.component, primal_edge.component);
}

} // namespace

void HodgeOperator::apply_inverse_star_1(
    const state::DualFaceFormStateView &,
    state::PrimalEdgeFormStateView) const
{
    throw std::logic_error("HodgeOperator does not support StateView inverse star-1");
}

LumpedHodge::LumpedHodge(LumpedHodgeCoefficients coefficients)
    : coefficients_(coefficients)
{
}

void LumpedHodge::apply_star_2(const state::PrimalFaceFormView &primal_face,
                               state::DualEdgeFormView dual_edge) const
{
    {
        const Int3 lo = dual_edge.xi.get_lo();
        const Int3 hi = dual_edge.xi.get_hi();
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                    dual_edge.xi(i, j, k, dual_edge.component) =
                        coefficients_.star_primal_2_to_dual_1.xi(
                            i, j, k, coefficients_.star_primal_2_to_dual_1.component) *
                        primal_face.xi(i, j, k, primal_face.component);
    }

    {
        const Int3 lo = dual_edge.eta.get_lo();
        const Int3 hi = dual_edge.eta.get_hi();
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                    dual_edge.eta(i, j, k, dual_edge.component) =
                        coefficients_.star_primal_2_to_dual_1.eta(
                            i, j, k, coefficients_.star_primal_2_to_dual_1.component) *
                        primal_face.eta(i, j, k, primal_face.component);
    }

    {
        const Int3 lo = dual_edge.zeta.get_lo();
        const Int3 hi = dual_edge.zeta.get_hi();
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
                        dual_face.xi(i, j, k, dual_face.component) *
                        coefficients_.inverse_star_primal_1.xi(
                            i, j, k, coefficients_.inverse_star_primal_1.component);
    }

    {
        const Int3 lo = primal_edge.eta.inner_lo();
        const Int3 hi = primal_edge.eta.inner_hi();
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                    primal_edge.eta(i, j, k, primal_edge.component) =
                        dual_face.eta(i, j, k, dual_face.component) *
                        coefficients_.inverse_star_primal_1.eta(
                            i, j, k, coefficients_.inverse_star_primal_1.component);
    }

    {
        const Int3 lo = primal_edge.zeta.inner_lo();
        const Int3 hi = primal_edge.zeta.inner_hi();
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                    primal_edge.zeta(i, j, k, primal_edge.component) =
                        dual_face.zeta(i, j, k, dual_face.component) *
                        coefficients_.inverse_star_primal_1.zeta(
                            i, j, k, coefficients_.inverse_star_primal_1.component);
    }
}

SnapshotLumpedHodge::SnapshotLumpedHodge(
    const geometry::LumpedHodgeGeometry &coefficients)
    : coefficients_(coefficients)
{
    if (!coefficients_.is_complete())
        throw std::logic_error("SnapshotLumpedHodge requires complete geometry coefficients");
}

void SnapshotLumpedHodge::apply_star_2(
    const state::PrimalFaceFormView &primal_face,
    state::DualEdgeFormView dual_edge) const
{
    apply_snapshot_star_2(primal_face, dual_edge, coefficients_);
}

void SnapshotLumpedHodge::apply_star_2(
    const state::PrimalFaceFormStateView &primal_face,
    state::DualEdgeFormStateView dual_edge) const
{
    apply_snapshot_star_2(primal_face, dual_edge, coefficients_);
}

void SnapshotLumpedHodge::apply_inverse_star_1(
    const state::DualFaceFormView &dual_face,
    state::PrimalEdgeFormView primal_edge) const
{
    apply_snapshot_inverse_star_1(dual_face, primal_edge, coefficients_);
}

void SnapshotLumpedHodge::apply_inverse_star_1(
    const state::DualFaceFormStateView &dual_face,
    state::PrimalEdgeFormStateView primal_edge) const
{
    apply_snapshot_inverse_star_1(dual_face, primal_edge, coefficients_);
}

} // namespace mpcns::discrete
