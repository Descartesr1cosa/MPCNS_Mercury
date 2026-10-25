#include "mpcns/compat/legacy_geometry_snapshot.hpp"

#include "mpcns/compat/legacy_geometry_field_ids.hpp"

#include "1_grid/1_MPCNS_Grid.h"
#include "3_field/Field.h"

#include <stdexcept>
#include <utility>
#include <vector>

namespace mpcns::compat
{

namespace
{

double read_field_block(const void *storage,
                        core::LogicalIndex index,
                        int component)
{
    const auto *field = static_cast<const FieldBlock *>(storage);
    return (*field)(index.i, index.j, index.k, component);
}

double read_coordinate(const void *storage,
                       core::LogicalIndex index,
                       int /* component */)
{
    const auto *coordinate = static_cast<const double3D *>(storage);
    return (*coordinate)(index.i, index.j, index.k);
}

core::LogicalBox field_extent(const FieldBlock &field)
{
    const Int3 lo = field.get_lo();
    const Int3 hi = field.get_hi();
    return core::LogicalBox{{lo.i, lo.j, lo.k}, {hi.i, hi.j, hi.k}};
}

core::LogicalBox coordinate_extent(const double3D &coordinate)
{
    const int ghost = coordinate.Getghostmesh();
    return core::LogicalBox{
        {-ghost, -ghost, -ghost},
        {coordinate.Getsize1() - ghost,
         coordinate.Getsize2() - ghost,
         coordinate.Getsize3() - ghost}};
}

geometry::ConstFieldView borrow_field(const FieldBlock &field,
                                      const char *semantic_name)
{
    if (!field.is_allocated())
        throw std::logic_error(std::string("cannot build GeometrySnapshot from inactive field: ") + semantic_name);
    return geometry::ConstFieldView::borrow(
        &field,
        read_field_block,
        field_extent(field),
        field.descriptor().ncomp);
}

geometry::ConstFieldView borrow_coordinate(const double3D &coordinate)
{
    return geometry::ConstFieldView::borrow(
        &coordinate,
        read_coordinate,
        coordinate_extent(coordinate),
        1);
}

geometry::DirectionalFieldViews borrow_directional(Field &fields,
                                                    DirectionalFieldIds ids,
                                                    int block,
                                                    const char *semantic_name)
{
    if (ids.xi < 0 || ids.eta < 0 || ids.zeta < 0)
        throw std::logic_error(std::string("missing legacy geometry field IDs for ") + semantic_name);

    return geometry::DirectionalFieldViews{
        borrow_field(fields.field(ids.xi, block), semantic_name),
        borrow_field(fields.field(ids.eta, block), semantic_name),
        borrow_field(fields.field(ids.zeta, block), semantic_name)};
}

} // namespace

geometry::GeometrySnapshot borrow_legacy_geometry_snapshot(Grid &grid,
                                                            Field &fields)
{
    if (grid.nblock != fields.num_blocks())
        throw std::logic_error("Grid and Field disagree on the number of local geometry blocks");

    const LegacyGeometryFieldIds ids = bind_legacy_geometry_field_ids(fields);
    if (ids.primal_cell_volume < 0)
        throw std::logic_error("missing legacy geometry field ID for primal_cell_volume");

    std::vector<geometry::GeometryBlock> blocks;
    blocks.reserve(grid.nblock);
    for (int block = 0; block < grid.nblock; ++block)
    {
        Block &legacy_block = grid.grids(block);
        blocks.push_back(geometry::GeometryBlock{
            legacy_block.block_name,
            legacy_block.dimension,
            geometry::NodeCoordinates{
                borrow_coordinate(legacy_block.x),
                borrow_coordinate(legacy_block.y),
                borrow_coordinate(legacy_block.z)},
            geometry::PrimalGeometry{
                borrow_field(fields.field(ids.primal_cell_volume, block), "primal_cell_volume"),
                borrow_directional(fields, ids.primal_face_area_vector, block, "primal_face_area_vector"),
                borrow_directional(fields, ids.primal_edge_vector, block, "primal_edge_vector"),
                borrow_directional(fields, ids.primal_face_area, block, "primal_face_area"),
                borrow_directional(fields, ids.primal_edge_length, block, "primal_edge_length")},
            geometry::DualGeometry{
                borrow_directional(fields, ids.dual_edge_length_at_face, block, "dual_edge_length_at_primal_face"),
                borrow_directional(fields, ids.dual_face_area_vector_at_edge, block, "dual_face_area_vector_at_primal_edge"),
                borrow_directional(fields, ids.dual_face_area_at_edge, block, "dual_face_area_at_primal_edge")},
            geometry::ReconstructionCache{
                borrow_directional(fields, ids.cell_covariant_basis, block, "cell_covariant_basis"),
                borrow_field(fields.field(ids.cell_reconstruction_gradient_pseudoinverse, block),
                             "cell_reconstruction_gradient_pseudoinverse"),
                borrow_directional(fields, ids.edge_reconstruction_gradient_pseudoinverse, block,
                                   "edge_reconstruction_gradient_pseudoinverse"),
                borrow_directional(fields, ids.edge_reconstruction_area_pseudoinverse, block,
                                   "edge_reconstruction_area_pseudoinverse")}});
    }
    return geometry::GeometrySnapshot{std::move(blocks)};
}

} // namespace mpcns::compat
