#include "mpcns/compat/legacy_singular_edge_collector.hpp"

#include "1_grid/1_MPCNS_Grid.h"
#include "3_field/Field.h"

#include <cmath>
#include <set>
#include <string>
#include <vector>

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

std::vector<TOPO::EntityKey> incident_cells(const TOPO::EntityKey &edge)
{
    std::vector<TOPO::EntityKey> result;
    const int axis = static_cast<int>(edge.axis);
    for (int first_offset : {-1, 0})
        for (int second_offset : {-1, 0})
        {
            int i = edge.i;
            int j = edge.j;
            int k = edge.k;
            if (axis == 0)
            {
                j += first_offset;
                k += second_offset;
            }
            if (axis == 1)
            {
                i += first_offset;
                k += second_offset;
            }
            if (axis == 2)
            {
                i += first_offset;
                j += second_offset;
            }
            result.push_back(TOPO::make_cell(edge.rank, edge.block, i, j, k));
        }
    return result;
}

std::vector<TOPO::EntityKey> incident_faces(const TOPO::EntityKey &edge)
{
    using TOPO::EntityAxis;
    std::vector<TOPO::EntityKey> result;
    if (edge.axis == EntityAxis::Xi)
    {
        result.push_back(TOPO::make_face(
            edge.rank, edge.block, edge.i, edge.j, edge.k, EntityAxis::Eta));
        result.push_back(TOPO::make_face(
            edge.rank, edge.block, edge.i, edge.j, edge.k - 1, EntityAxis::Eta));
        result.push_back(TOPO::make_face(
            edge.rank, edge.block, edge.i, edge.j, edge.k, EntityAxis::Zeta));
        result.push_back(TOPO::make_face(
            edge.rank, edge.block, edge.i, edge.j - 1, edge.k, EntityAxis::Zeta));
    }
    else if (edge.axis == EntityAxis::Eta)
    {
        result.push_back(TOPO::make_face(
            edge.rank, edge.block, edge.i, edge.j, edge.k, EntityAxis::Xi));
        result.push_back(TOPO::make_face(
            edge.rank, edge.block, edge.i, edge.j, edge.k - 1, EntityAxis::Xi));
        result.push_back(TOPO::make_face(
            edge.rank, edge.block, edge.i, edge.j, edge.k, EntityAxis::Zeta));
        result.push_back(TOPO::make_face(
            edge.rank, edge.block, edge.i - 1, edge.j, edge.k, EntityAxis::Zeta));
    }
    else
    {
        result.push_back(TOPO::make_face(
            edge.rank, edge.block, edge.i, edge.j, edge.k, EntityAxis::Xi));
        result.push_back(TOPO::make_face(
            edge.rank, edge.block, edge.i, edge.j - 1, edge.k, EntityAxis::Xi));
        result.push_back(TOPO::make_face(
            edge.rank, edge.block, edge.i, edge.j, edge.k, EntityAxis::Eta));
        result.push_back(TOPO::make_face(
            edge.rank, edge.block, edge.i - 1, edge.j, edge.k, EntityAxis::Eta));
    }
    return result;
}

const char *axis_suffix(TOPO::EntityAxis axis)
{
    return axis == TOPO::EntityAxis::Xi
               ? "xi"
               : (axis == TOPO::EntityAxis::Eta ? "eta" : "zeta");
}

} // namespace

LegacySingularEdgeLocalMeasures collect_legacy_singular_edge_measures(
    const TOPO::Topology &topology,
    const LegacySingularEdgeDataSource &source,
    int local_rank,
    const std::vector<mesh::SingularEdgeAlias> &aliases)
{
    LegacySingularEdgeLocalMeasures result;
    std::set<TOPO::EntityKey> seen_cells;
    std::set<TOPO::EntityKey> seen_faces;

    for (const mesh::SingularEdgeAlias &alias : aliases)
    {
        const TOPO::EntityKey &edge = alias.edge;
        if (edge.rank != local_rank ||
            edge.block < 0 ||
            edge.block >= source.fields.num_blocks())
            continue;

        FieldBlock &edge_length = source.fields.field(
            std::string("dl_") + axis_suffix(edge.axis), edge.block);
        if (contains_inner(edge_length, edge.i, edge.j, edge.k))
        {
            const double length = edge_length(edge.i, edge.j, edge.k, 0);
            if (std::isfinite(length) && length > 0.0)
            {
                result.primal_length_sum += length;
                ++result.primal_length_count;
            }

            FieldBlock &edge_vector = source.fields.field(
                std::string("dr_") + axis_suffix(edge.axis), edge.block);
            for (int component = 0; component < 3; ++component)
                result.tangent_sum[component] +=
                    alias.orientation *
                    edge_vector(edge.i, edge.j, edge.k, component);
        }

        FieldBlock &cell_volume = source.fields.field("Jac", edge.block);
        for (const TOPO::EntityKey &cell : incident_cells(edge))
            if (seen_cells.insert(cell).second &&
                contains_inner(cell_volume, cell.i, cell.j, cell.k))
            {
                const double measure =
                    std::abs(cell_volume(cell.i, cell.j, cell.k, 0));
                if (std::isfinite(measure) && measure > 1.e-30)
                {
                    result.incident_cells.push_back(
                        {cell, measure, 0.0, edge, alias.orientation});
                    result.cell_measure_sum += measure;
                    ++result.valid_cell_count;

                    Block &block = source.grid.grids(cell.block);
                    result.incident_cell_centers.push_back(
                        {block.dual_x(cell.i + 1, cell.j + 1, cell.k + 1),
                         block.dual_y(cell.i + 1, cell.j + 1, cell.k + 1),
                         block.dual_z(cell.i + 1, cell.j + 1, cell.k + 1)});
                }
            }

        for (const TOPO::EntityKey &face : incident_faces(edge))
            if (seen_faces.insert(face).second)
            {
                // Count and assemble each quotient (physical) face once.
                if (topology.faces.local_to_qkey.find(face) ==
                    topology.faces.local_to_qkey.end())
                    continue;
                if (!topology.is_owner(face))
                    continue;

                FieldBlock &face_area = source.fields.field(
                    std::string("Area_") + axis_suffix(face.axis), face.block);
                if (!contains_inner(face_area, face.i, face.j, face.k))
                    continue;

                const double measure = face_area(face.i, face.j, face.k, 0);
                if (std::isfinite(measure) && measure > 1.e-30)
                {
                    result.incident_faces.push_back(
                        {face, measure, 0.0, edge, alias.orientation});
                    result.face_measure_sum += measure;
                    ++result.valid_face_count;
                }
            }
    }
    return result;
}

} // namespace mpcns::compat
