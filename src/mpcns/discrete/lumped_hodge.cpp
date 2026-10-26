#include "mpcns/discrete/lumped_hodge.hpp"

#include <stdexcept>

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
    const auto apply = [](const FieldBlock &source, FieldBlock &target,
                          const geometry::ConstFieldView &coefficient,
                          int source_component, int target_component)
    {
        const Int3 lo = target.get_lo();
        const Int3 hi = target.get_hi();
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                {
                    const core::LogicalIndex index{i, j, k};
                    target(i, j, k, target_component) =
                        coefficient.value(index) * source(i, j, k, source_component);
                }
    };
    apply(primal_face.xi, dual_edge.xi,
          coefficients_.star_primal_2_to_dual_1.xi,
          primal_face.component, dual_edge.component);
    apply(primal_face.eta, dual_edge.eta,
          coefficients_.star_primal_2_to_dual_1.eta,
          primal_face.component, dual_edge.component);
    apply(primal_face.zeta, dual_edge.zeta,
          coefficients_.star_primal_2_to_dual_1.zeta,
          primal_face.component, dual_edge.component);
}

void SnapshotLumpedHodge::apply_inverse_star_1(
    const state::DualFaceFormView &dual_face,
    state::PrimalEdgeFormView primal_edge) const
{
    const auto apply = [](const FieldBlock &source, FieldBlock &target,
                          const geometry::ConstFieldView &coefficient,
                          int source_component, int target_component)
    {
        const Int3 lo = target.inner_lo();
        const Int3 hi = target.inner_hi();
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                {
                    const core::LogicalIndex index{i, j, k};
                    target(i, j, k, target_component) =
                        source(i, j, k, source_component) * coefficient.value(index);
                }
    };
    apply(dual_face.xi, primal_edge.xi,
          coefficients_.inverse_star_primal_1.xi,
          dual_face.component, primal_edge.component);
    apply(dual_face.eta, primal_edge.eta,
          coefficients_.inverse_star_primal_1.eta,
          dual_face.component, primal_edge.component);
    apply(dual_face.zeta, primal_edge.zeta,
          coefficients_.inverse_star_primal_1.zeta,
          dual_face.component, primal_edge.component);
}

} // namespace mpcns::discrete
