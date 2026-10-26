#pragma once

#include "mpcns/geometry/geometry_snapshot.hpp"
#include "mpcns/state/cochain_view.hpp"
#include "mpcns/compat/legacy_cochain_view.hpp"

namespace mpcns::discrete
{

class HodgeOperator
{
public:
    virtual ~HodgeOperator() = default;

    virtual void apply_star_2(const state::PrimalFaceFormView &primal_face,
                              state::DualEdgeFormView dual_edge) const = 0;
    virtual void apply_inverse_star_1(const state::DualFaceFormView &dual_face,
                                      state::PrimalEdgeFormView primal_edge) const = 0;
};

struct LumpedHodgeCoefficients
{
    // |dual edge| / |primal face|, stored at primal faces.
    state::PrimalFaceScalarView star_primal_2_to_dual_1;
    // |primal edge| / |dual face|, stored at primal edges.
    state::PrimalEdgeScalarView inverse_star_primal_1;
};

// Production diagonal Hodge.  It exposes the same mathematical operations
// that a future consistent/non-diagonal Hodge must implement.
class LumpedHodge final : public HodgeOperator
{
public:
    explicit LumpedHodge(LumpedHodgeCoefficients coefficients);

    void apply_star_2(const state::PrimalFaceFormView &primal_face,
                      state::DualEdgeFormView dual_edge) const override;
    void apply_inverse_star_1(const state::DualFaceFormView &dual_face,
                              state::PrimalEdgeFormView primal_edge) const override;

private:
    LumpedHodgeCoefficients coefficients_;
};

// Snapshot-backed production Hodge.  It reads immutable geometry through the
// typed GeometrySnapshot boundary, so operator clients need not look up
// legacy Field catalogue IDs.  Cochains remain transitional FieldBlock views.
class SnapshotLumpedHodge final : public HodgeOperator
{
public:
    explicit SnapshotLumpedHodge(
        const geometry::LumpedHodgeGeometry &coefficients);

    void apply_star_2(const state::PrimalFaceFormView &primal_face,
                      state::DualEdgeFormView dual_edge) const override;
    void apply_inverse_star_1(const state::DualFaceFormView &dual_face,
                              state::PrimalEdgeFormView primal_edge) const override;

private:
    const geometry::LumpedHodgeGeometry &coefficients_;
};

} // namespace mpcns::discrete
