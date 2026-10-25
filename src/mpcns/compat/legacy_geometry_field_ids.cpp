#include "mpcns/compat/legacy_geometry_field_ids.hpp"
#include "mpcns/compat/legacy_geometry_field_names.hpp"

#include "3_field/Field.h"

namespace mpcns::compat
{

namespace names = legacy_geometry_names;

namespace
{

DirectionalFieldIds directional_field_ids(const Field &fields,
                                          const char *xi_name,
                                          const char *eta_name,
                                          const char *zeta_name)
{
    return DirectionalFieldIds{
        fields.field_id(xi_name),
        fields.field_id(eta_name),
        fields.field_id(zeta_name)};
}

} // namespace

LegacyGeometryFieldIds bind_legacy_geometry_field_ids(const Field &fields)
{
    LegacyGeometryFieldIds ids;
    ids.primal_cell_volume = fields.field_id(names::kPrimalCellVolume);
    ids.primal_face_area_vector =
        directional_field_ids(fields, names::kPrimalFaceAreaVectorXi,
                              names::kPrimalFaceAreaVectorEta,
                              names::kPrimalFaceAreaVectorZeta);
    ids.primal_edge_vector =
        directional_field_ids(fields, names::kPrimalEdgeVectorXi,
                              names::kPrimalEdgeVectorEta,
                              names::kPrimalEdgeVectorZeta);
    ids.cell_covariant_basis =
        directional_field_ids(fields, names::kCellCovariantBasisXi,
                              names::kCellCovariantBasisEta,
                              names::kCellCovariantBasisZeta);

    ids.primal_face_area =
        directional_field_ids(fields, names::kPrimalFaceAreaXi,
                              names::kPrimalFaceAreaEta,
                              names::kPrimalFaceAreaZeta);
    ids.dual_edge_length_at_face =
        directional_field_ids(fields, names::kDualEdgeLengthAtPrimalFaceXi,
                              names::kDualEdgeLengthAtPrimalFaceEta,
                              names::kDualEdgeLengthAtPrimalFaceZeta);
    ids.primal_edge_length =
        directional_field_ids(fields, names::kPrimalEdgeLengthXi,
                              names::kPrimalEdgeLengthEta,
                              names::kPrimalEdgeLengthZeta);
    ids.dual_face_area_vector_at_edge =
        directional_field_ids(fields, names::kDualFaceAreaVectorAtPrimalEdgeXi,
                              names::kDualFaceAreaVectorAtPrimalEdgeEta,
                              names::kDualFaceAreaVectorAtPrimalEdgeZeta);
    ids.dual_face_area_at_edge =
        directional_field_ids(fields, names::kDualFaceAreaAtPrimalEdgeXi,
                              names::kDualFaceAreaAtPrimalEdgeEta,
                              names::kDualFaceAreaAtPrimalEdgeZeta);

    ids.cell_reconstruction_gradient_pseudoinverse =
        fields.field_id(names::kCellGradientPseudoinverse);
    ids.edge_reconstruction_gradient_pseudoinverse =
        directional_field_ids(fields, names::kEdgeGradientPseudoinverseXi,
                              names::kEdgeGradientPseudoinverseEta,
                              names::kEdgeGradientPseudoinverseZeta);
    ids.edge_reconstruction_area_pseudoinverse =
        directional_field_ids(fields, names::kEdgeAreaPseudoinverseXi,
                              names::kEdgeAreaPseudoinverseEta,
                              names::kEdgeAreaPseudoinverseZeta);

    ids.lumped_star_primal_0_to_dual_3 =
        fields.field_id(names::kLumpedStarPrimal0ToDual3);
    ids.lumped_star_primal_3_to_dual_0 =
        fields.field_id(names::kLumpedStarPrimal3ToDual0);
    ids.lumped_star_primal_1_to_dual_2 = directional_field_ids(
        fields,
        names::kLumpedStarPrimal1ToDual2Xi,
        names::kLumpedStarPrimal1ToDual2Eta,
        names::kLumpedStarPrimal1ToDual2Zeta);
    ids.lumped_star_primal_2_to_dual_1 = directional_field_ids(
        fields,
        names::kLumpedStarPrimal2ToDual1Xi,
        names::kLumpedStarPrimal2ToDual1Eta,
        names::kLumpedStarPrimal2ToDual1Zeta);
    ids.lumped_inverse_star_primal_1 = directional_field_ids(
        fields,
        names::kLumpedInverseStarPrimal1Xi,
        names::kLumpedInverseStarPrimal1Eta,
        names::kLumpedInverseStarPrimal1Zeta);
    return ids;
}

discrete::LumpedHodgeCoefficients bind_legacy_lumped_hodge_coefficients(
    Field &fields,
    int block)
{
    const LegacyGeometryFieldIds ids = bind_legacy_geometry_field_ids(fields);
    return discrete::LumpedHodgeCoefficients{
        state::PrimalFaceScalarView{
            fields.field(ids.lumped_star_primal_2_to_dual_1.xi, block),
            fields.field(ids.lumped_star_primal_2_to_dual_1.eta, block),
            fields.field(ids.lumped_star_primal_2_to_dual_1.zeta, block)},
        state::PrimalEdgeScalarView{
            fields.field(ids.lumped_inverse_star_primal_1.xi, block),
            fields.field(ids.lumped_inverse_star_primal_1.eta, block),
            fields.field(ids.lumped_inverse_star_primal_1.zeta, block)}};
}

} // namespace mpcns::compat
