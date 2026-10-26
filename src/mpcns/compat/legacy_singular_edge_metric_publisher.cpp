#include "mpcns/compat/legacy_singular_edge_metric_publisher.hpp"

#include "3_field/Field.h"
#include "7_metric/SingularEdgeRegistry.h"

#include <string>

namespace mpcns::compat
{
namespace
{

bool contains_inner(const FieldBlock &field, int i, int j, int k)
{
    const Int3 lo = field.inner_lo();
    const Int3 hi = field.inner_hi();
    return i >= lo.i && i < hi.i &&
           j >= lo.j && j < hi.j &&
           k >= lo.k && k < hi.k;
}

const char *axis_suffix(TOPO::EntityAxis axis)
{
    return axis == TOPO::EntityAxis::Xi
               ? "xi"
               : (axis == TOPO::EntityAxis::Eta ? "eta" : "zeta");
}

} // namespace

void publish_legacy_singular_edge_hodge(
    Field &fields,
    const METRIC::SingularPhysicalEdge &edge,
    int local_rank)
{
    // The vector-valued Sstar remains untouched until a dual-polygon vector
    // reconstruction is requested by an operator.
    for (const METRIC::SingularEdgeAlias &alias : edge.aliases)
    {
        const TOPO::EntityKey &entity = alias.edge;
        if (entity.rank != local_rank ||
            entity.block < 0 ||
            entity.block >= fields.num_blocks())
            continue;

        FieldBlock &dual_area =
            fields.field(std::string("Astar_") + axis_suffix(entity.axis),
                         entity.block);
        FieldBlock &inverse_hodge =
            fields.field(
                std::string("Hodge_star_inverse_2form_to_1form_edge_") +
                    axis_suffix(entity.axis) + "_lumped",
                entity.block);
        if (contains_inner(dual_area, entity.i, entity.j, entity.k))
            dual_area(entity.i, entity.j, entity.k, 0) = edge.dual_area;
        if (contains_inner(inverse_hodge, entity.i, entity.j, entity.k))
            inverse_hodge(entity.i, entity.j, entity.k, 0) =
                edge.inverse_hodge;
    }
}

} // namespace mpcns::compat
