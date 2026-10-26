#include "mpcns/compat/legacy_geometry_field_registration.hpp"

#include "mpcns/compat/legacy_geometry_field_names.hpp"

#include "3_field/Field.h"

namespace mpcns::compat
{
namespace names = legacy_geometry_names;

namespace
{

using GhostWidth = LegacyGeometryGhostWidth;

const std::array<LegacyGeometryFieldRegistration,
                 kLegacyGeometryFieldRegistrationCount>
    kLegacyGeometryFieldRegistrations = {{
        {names::kPrimalCellVolume, StaggerLocation::Cell, 1, GhostWidth::GeometryMinusOne},
        {names::kPrimalFaceAreaVectorXi, StaggerLocation::FaceXi, 3, GhostWidth::Geometry},
        {names::kPrimalFaceAreaVectorEta, StaggerLocation::FaceEt, 3, GhostWidth::Geometry},
        {names::kPrimalFaceAreaVectorZeta, StaggerLocation::FaceZe, 3, GhostWidth::Geometry},
        {names::kCellCovariantBasisXi, StaggerLocation::Cell, 3, GhostWidth::Zero},
        {names::kCellCovariantBasisEta, StaggerLocation::Cell, 3, GhostWidth::Zero},
        {names::kCellCovariantBasisZeta, StaggerLocation::Cell, 3, GhostWidth::Zero},
        {names::kPrimalEdgeVectorXi, StaggerLocation::EdgeXi, 3, GhostWidth::GeometryMinusOne},
        {names::kPrimalEdgeVectorEta, StaggerLocation::EdgeEt, 3, GhostWidth::GeometryMinusOne},
        {names::kPrimalEdgeVectorZeta, StaggerLocation::EdgeZe, 3, GhostWidth::GeometryMinusOne},
        {names::kEdgeGradientPseudoinverseXi, StaggerLocation::EdgeXi, 9, GhostWidth::Zero},
        {names::kEdgeAreaPseudoinverseXi, StaggerLocation::EdgeXi, 9, GhostWidth::Zero},
        {names::kEdgeGradientPseudoinverseEta, StaggerLocation::EdgeEt, 9, GhostWidth::Zero},
        {names::kEdgeAreaPseudoinverseEta, StaggerLocation::EdgeEt, 9, GhostWidth::Zero},
        {names::kEdgeGradientPseudoinverseZeta, StaggerLocation::EdgeZe, 9, GhostWidth::Zero},
        {names::kEdgeAreaPseudoinverseZeta, StaggerLocation::EdgeZe, 9, GhostWidth::Zero},
        {names::kCellGradientPseudoinverse, StaggerLocation::Cell, 9, GhostWidth::Zero},
        {names::kPrimalFaceAreaXi, StaggerLocation::FaceXi, 1, GhostWidth::GeometryMinusOne},
        {names::kPrimalFaceAreaEta, StaggerLocation::FaceEt, 1, GhostWidth::GeometryMinusOne},
        {names::kPrimalFaceAreaZeta, StaggerLocation::FaceZe, 1, GhostWidth::GeometryMinusOne},
        {names::kDualEdgeLengthAtPrimalFaceXi, StaggerLocation::FaceXi, 1, GhostWidth::GeometryMinusOne},
        {names::kDualEdgeLengthAtPrimalFaceEta, StaggerLocation::FaceEt, 1, GhostWidth::GeometryMinusOne},
        {names::kDualEdgeLengthAtPrimalFaceZeta, StaggerLocation::FaceZe, 1, GhostWidth::GeometryMinusOne},
        {names::kPrimalEdgeLengthXi, StaggerLocation::EdgeXi, 1, GhostWidth::GeometryMinusOne},
        {names::kPrimalEdgeLengthEta, StaggerLocation::EdgeEt, 1, GhostWidth::GeometryMinusOne},
        {names::kPrimalEdgeLengthZeta, StaggerLocation::EdgeZe, 1, GhostWidth::GeometryMinusOne},
        {names::kDualFaceAreaVectorAtPrimalEdgeXi, StaggerLocation::EdgeXi, 3, GhostWidth::GeometryMinusOne},
        {names::kDualFaceAreaVectorAtPrimalEdgeEta, StaggerLocation::EdgeEt, 3, GhostWidth::GeometryMinusOne},
        {names::kDualFaceAreaVectorAtPrimalEdgeZeta, StaggerLocation::EdgeZe, 3, GhostWidth::GeometryMinusOne},
        {names::kDualFaceAreaAtPrimalEdgeXi, StaggerLocation::EdgeXi, 1, GhostWidth::GeometryMinusOne},
        {names::kDualFaceAreaAtPrimalEdgeEta, StaggerLocation::EdgeEt, 1, GhostWidth::GeometryMinusOne},
        {names::kDualFaceAreaAtPrimalEdgeZeta, StaggerLocation::EdgeZe, 1, GhostWidth::GeometryMinusOne},
        {names::kLumpedStarPrimal0ToDual3, StaggerLocation::Cell, 1, GhostWidth::GeometryMinusOne},
        {names::kLumpedStarPrimal3ToDual0, StaggerLocation::Cell, 1, GhostWidth::GeometryMinusOne},
        {names::kLumpedStarPrimal1ToDual2Xi, StaggerLocation::EdgeXi, 1, GhostWidth::GeometryMinusOne},
        {names::kLumpedStarPrimal1ToDual2Eta, StaggerLocation::EdgeEt, 1, GhostWidth::GeometryMinusOne},
        {names::kLumpedStarPrimal1ToDual2Zeta, StaggerLocation::EdgeZe, 1, GhostWidth::GeometryMinusOne},
        {names::kLumpedStarPrimal2ToDual1Xi, StaggerLocation::FaceXi, 1, GhostWidth::GeometryMinusOne},
        {names::kLumpedStarPrimal2ToDual1Eta, StaggerLocation::FaceEt, 1, GhostWidth::GeometryMinusOne},
        {names::kLumpedStarPrimal2ToDual1Zeta, StaggerLocation::FaceZe, 1, GhostWidth::GeometryMinusOne},
        {names::kLumpedInverseStarPrimal1Xi, StaggerLocation::EdgeXi, 1, GhostWidth::GeometryMinusOne},
        {names::kLumpedInverseStarPrimal1Eta, StaggerLocation::EdgeEt, 1, GhostWidth::GeometryMinusOne},
        {names::kLumpedInverseStarPrimal1Zeta, StaggerLocation::EdgeZe, 1, GhostWidth::GeometryMinusOne},
    }};

} // namespace

int LegacyGeometryFieldRegistration::ghost_layers(int geometry_ghost) const
{
    switch (ghost_width)
    {
    case LegacyGeometryGhostWidth::Geometry:
        return geometry_ghost;
    case LegacyGeometryGhostWidth::GeometryMinusOne:
        return geometry_ghost - 1;
    case LegacyGeometryGhostWidth::Zero:
        return 0;
    }
    return 0;
}

const std::array<LegacyGeometryFieldRegistration,
                 kLegacyGeometryFieldRegistrationCount> &
legacy_geometry_field_registrations()
{
    return kLegacyGeometryFieldRegistrations;
}

void register_legacy_geometry_fields(Field &fields, int geometry_ghost)
{
    for (const LegacyGeometryFieldRegistration &registration :
         legacy_geometry_field_registrations())
    {
        fields.register_field(FieldDescriptor{
            registration.legacy_name,
            registration.location,
            registration.component_count,
            registration.ghost_layers(geometry_ghost)});
    }
}

} // namespace mpcns::compat
