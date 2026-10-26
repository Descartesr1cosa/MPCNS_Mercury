#include "mpcns/compat/legacy_geometry_build_pipeline.hpp"

#include "mpcns/compat/legacy_geometry_field_registration.hpp"
#include "mpcns/compat/legacy_geometry_snapshot.hpp"
#include "mpcns/compat/legacy_lumped_metric_publisher.hpp"
#include "mpcns/compat/legacy_primal_metric_publisher.hpp"
#include "mpcns/compat/legacy_reconstruction_metric_publisher.hpp"

#include <utility>

namespace mpcns::compat
{

void LegacyGeometryBuildPipeline::register_fields(
    Field &fields,
    int geometry_ghost) const
{
    register_legacy_geometry_fields(fields, geometry_ghost);
}

void LegacyGeometryBuildPipeline::publish(
    const LegacyGeometryBuildRequest &request) const
{
    // This order is an explicit compatibility contract. New geometry
    // backends replace this producer, not consumers or the legacy catalog.
    populate_legacy_primal_face_and_cell_metrics(request.fields, request.grid);
    populate_legacy_lumped_face_and_edge_metrics(request.fields, request.grid);
    populate_legacy_reconstruction_cache(request.fields, request.grid);
    populate_legacy_primal_edge_vectors(request.fields, request.grid);
    populate_legacy_diagonal_hodge_closure(request.fields, request.grid);
}

geometry::GeometryBuildResult LegacyGeometryBuildPipeline::finalize(
    const LegacyGeometryBuildRequest &request) const
{
    const geometry::GeometrySnapshot borrowed =
        borrow_legacy_geometry_snapshot(request.grid, request.fields);

    geometry::GeometryBuildInput input;
    input.blocks = borrowed.blocks();
    return geometry::GeometryBuilder{}.finalize(std::move(input));
}

geometry::GeometryBuildResult LegacyGeometryBuildPipeline::build(
    const LegacyGeometryBuildRequest &request,
    int geometry_ghost) const
{
    register_fields(request.fields, geometry_ghost);
    publish(request);
    return finalize(request);
}

} // namespace mpcns::compat
