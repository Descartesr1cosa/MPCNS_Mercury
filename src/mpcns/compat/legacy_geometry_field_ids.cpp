#include "mpcns/compat/legacy_geometry_field_ids.hpp"

#include "3_field/Field.h"

namespace mpcns::compat
{

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
    ids.primal_cell_volume = fields.field_id("Jac");
    ids.primal_face_area_vector =
        directional_field_ids(fields, "JDxi", "JDet", "JDze");
    ids.primal_edge_vector =
        directional_field_ids(fields, "dr_xi", "dr_eta", "dr_zeta");
    ids.cell_covariant_basis =
        directional_field_ids(fields, "a_xi", "a_eta", "a_zeta");

    ids.primal_face_area =
        directional_field_ids(fields, "Area_xi", "Area_eta", "Area_zeta");
    ids.dual_edge_length_at_face =
        directional_field_ids(fields, "dlstar_xi", "dlstar_eta", "dlstar_zeta");
    ids.primal_edge_length =
        directional_field_ids(fields, "dl_xi", "dl_eta", "dl_zeta");
    ids.dual_face_area_vector_at_edge =
        directional_field_ids(fields, "Sstar_xi", "Sstar_eta", "Sstar_zeta");
    ids.dual_face_area_at_edge =
        directional_field_ids(fields, "Astar_xi", "Astar_eta", "Astar_zeta");

    ids.cell_reconstruction_gradient_pseudoinverse =
        fields.field_id("pinvGT_cell");
    ids.edge_reconstruction_gradient_pseudoinverse =
        directional_field_ids(fields, "pinvGT_xi", "pinvGT_eta", "pinvGT_zeta");
    ids.edge_reconstruction_area_pseudoinverse =
        directional_field_ids(fields, "pinvAT_xi", "pinvAT_eta", "pinvAT_zeta");

    ids.lumped_star_primal_0_to_dual_3 =
        fields.field_id("Hodge_star_0form_to_3form_cell_lumped");
    ids.lumped_star_primal_3_to_dual_0 =
        fields.field_id("Hodge_star_3form_to_0form_cell_lumped");
    ids.lumped_star_primal_1_to_dual_2 = directional_field_ids(
        fields,
        "Hodge_star_1form_to_2form_edge_xi_lumped",
        "Hodge_star_1form_to_2form_edge_eta_lumped",
        "Hodge_star_1form_to_2form_edge_zeta_lumped");
    ids.lumped_star_primal_2_to_dual_1 = directional_field_ids(
        fields,
        "Hodge_star_2form_to_1form_face_xi_lumped",
        "Hodge_star_2form_to_1form_face_eta_lumped",
        "Hodge_star_2form_to_1form_face_zeta_lumped");
    ids.lumped_inverse_star_primal_1 = directional_field_ids(
        fields,
        "Hodge_star_inverse_2form_to_1form_edge_xi_lumped",
        "Hodge_star_inverse_2form_to_1form_edge_eta_lumped",
        "Hodge_star_inverse_2form_to_1form_edge_zeta_lumped");
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
