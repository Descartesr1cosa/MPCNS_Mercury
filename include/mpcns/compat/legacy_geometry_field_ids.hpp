#pragma once

#include "mpcns/discrete/lumped_hodge.hpp"

class Field;

namespace mpcns::compat
{

// Canonical semantic identities for geometry data that still reside in the
// legacy Field catalog.  The old catalog strings occur only in the adapter
// implementation; applications bind the IDs through this type.
struct DirectionalFieldIds
{
    int xi = -1;
    int eta = -1;
    int zeta = -1;
};

struct LegacyGeometryFieldIds
{
    int primal_cell_volume = -1;
    DirectionalFieldIds primal_face_area_vector;
    DirectionalFieldIds primal_edge_vector;
    DirectionalFieldIds cell_covariant_basis;

    DirectionalFieldIds primal_face_area;
    DirectionalFieldIds dual_edge_length_at_face;
    DirectionalFieldIds primal_edge_length;
    DirectionalFieldIds dual_face_area_vector_at_edge;
    DirectionalFieldIds dual_face_area_at_edge;

    int cell_reconstruction_gradient_pseudoinverse = -1;
    DirectionalFieldIds edge_reconstruction_gradient_pseudoinverse;
    DirectionalFieldIds edge_reconstruction_area_pseudoinverse;

    int lumped_star_primal_0_to_dual_3 = -1;
    int lumped_star_primal_3_to_dual_0 = -1;
    DirectionalFieldIds lumped_star_primal_1_to_dual_2;
    DirectionalFieldIds lumped_star_primal_2_to_dual_1;
    // |primal edge| / |dual face|, i.e. the inverse of the primal 1-form
    // Hodge action, stored at primal edges.
    DirectionalFieldIds lumped_inverse_star_primal_1;
};

LegacyGeometryFieldIds bind_legacy_geometry_field_ids(const Field &fields);

// Borrow the production diagonal Hodge coefficients for one local block.
// The returned views retain Field-owned storage; this is a bridge, not a copy
// or a new geometry calculation.
discrete::LumpedHodgeCoefficients bind_legacy_lumped_hodge_coefficients(
    Field &fields,
    int block);

} // namespace mpcns::compat
