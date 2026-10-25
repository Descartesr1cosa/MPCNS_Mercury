#pragma once

// The legacy Field catalog remains a temporary backend.  Keep every geometry
// catalog spelling here so registration and compatibility binding cannot drift
// independently while the storage migration is in progress.
namespace mpcns::compat::legacy_geometry_names
{

constexpr const char kPrimalCellVolume[] = "Jac";

constexpr const char kPrimalFaceAreaVectorXi[] = "JDxi";
constexpr const char kPrimalFaceAreaVectorEta[] = "JDet";
constexpr const char kPrimalFaceAreaVectorZeta[] = "JDze";

constexpr const char kCellCovariantBasisXi[] = "a_xi";
constexpr const char kCellCovariantBasisEta[] = "a_eta";
constexpr const char kCellCovariantBasisZeta[] = "a_zeta";

constexpr const char kPrimalEdgeVectorXi[] = "dr_xi";
constexpr const char kPrimalEdgeVectorEta[] = "dr_eta";
constexpr const char kPrimalEdgeVectorZeta[] = "dr_zeta";

constexpr const char kCellGradientPseudoinverse[] = "pinvGT_cell";
constexpr const char kEdgeGradientPseudoinverseXi[] = "pinvGT_xi";
constexpr const char kEdgeGradientPseudoinverseEta[] = "pinvGT_eta";
constexpr const char kEdgeGradientPseudoinverseZeta[] = "pinvGT_zeta";
constexpr const char kEdgeAreaPseudoinverseXi[] = "pinvAT_xi";
constexpr const char kEdgeAreaPseudoinverseEta[] = "pinvAT_eta";
constexpr const char kEdgeAreaPseudoinverseZeta[] = "pinvAT_zeta";

constexpr const char kPrimalFaceAreaXi[] = "Area_xi";
constexpr const char kPrimalFaceAreaEta[] = "Area_eta";
constexpr const char kPrimalFaceAreaZeta[] = "Area_zeta";
constexpr const char kDualEdgeLengthAtPrimalFaceXi[] = "dlstar_xi";
constexpr const char kDualEdgeLengthAtPrimalFaceEta[] = "dlstar_eta";
constexpr const char kDualEdgeLengthAtPrimalFaceZeta[] = "dlstar_zeta";

constexpr const char kPrimalEdgeLengthXi[] = "dl_xi";
constexpr const char kPrimalEdgeLengthEta[] = "dl_eta";
constexpr const char kPrimalEdgeLengthZeta[] = "dl_zeta";
constexpr const char kDualFaceAreaVectorAtPrimalEdgeXi[] = "Sstar_xi";
constexpr const char kDualFaceAreaVectorAtPrimalEdgeEta[] = "Sstar_eta";
constexpr const char kDualFaceAreaVectorAtPrimalEdgeZeta[] = "Sstar_zeta";
constexpr const char kDualFaceAreaAtPrimalEdgeXi[] = "Astar_xi";
constexpr const char kDualFaceAreaAtPrimalEdgeEta[] = "Astar_eta";
constexpr const char kDualFaceAreaAtPrimalEdgeZeta[] = "Astar_zeta";

constexpr const char kLumpedStarPrimal0ToDual3[] = "Hodge_star_0form_to_3form_cell_lumped";
constexpr const char kLumpedStarPrimal3ToDual0[] = "Hodge_star_3form_to_0form_cell_lumped";
constexpr const char kLumpedStarPrimal1ToDual2Xi[] = "Hodge_star_1form_to_2form_edge_xi_lumped";
constexpr const char kLumpedStarPrimal1ToDual2Eta[] = "Hodge_star_1form_to_2form_edge_eta_lumped";
constexpr const char kLumpedStarPrimal1ToDual2Zeta[] = "Hodge_star_1form_to_2form_edge_zeta_lumped";
constexpr const char kLumpedStarPrimal2ToDual1Xi[] = "Hodge_star_2form_to_1form_face_xi_lumped";
constexpr const char kLumpedStarPrimal2ToDual1Eta[] = "Hodge_star_2form_to_1form_face_eta_lumped";
constexpr const char kLumpedStarPrimal2ToDual1Zeta[] = "Hodge_star_2form_to_1form_face_zeta_lumped";
constexpr const char kLumpedInverseStarPrimal1Xi[] = "Hodge_star_inverse_2form_to_1form_edge_xi_lumped";
constexpr const char kLumpedInverseStarPrimal1Eta[] = "Hodge_star_inverse_2form_to_1form_edge_eta_lumped";
constexpr const char kLumpedInverseStarPrimal1Zeta[] = "Hodge_star_inverse_2form_to_1form_edge_zeta_lumped";

} // namespace mpcns::compat::legacy_geometry_names
