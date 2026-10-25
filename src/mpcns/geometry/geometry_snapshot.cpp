#include "mpcns/geometry/geometry_snapshot.hpp"

#include <stdexcept>
#include <utility>

namespace mpcns::geometry
{

ConstFieldView ConstFieldView::borrow(const void *storage,
                                      ReadAccessor read,
                                      core::LogicalBox extent,
                                      int component_count)
{
    ConstFieldView view;
    view.storage_ = storage;
    view.read_ = read;
    view.extent_ = extent;
    view.component_count_ = component_count;
    return view;
}

bool ConstFieldView::is_bound() const
{
    return storage_ != nullptr && read_ != nullptr && component_count_ > 0;
}

double ConstFieldView::value(core::LogicalIndex index, int component) const
{
    if (!is_bound())
        throw std::logic_error("attempted to read an unbound geometry view");
    if (!extent_.contains(index))
        throw std::out_of_range("geometry index lies outside the view extent");
    if (component < 0 || component >= component_count_)
        throw std::out_of_range("geometry component lies outside the view extent");
    return read_(storage_, index, component);
}

const ConstFieldView &DirectionalFieldViews::at(core::Axis axis) const
{
    switch (axis)
    {
    case core::Axis::Xi:
        return xi;
    case core::Axis::Eta:
        return eta;
    case core::Axis::Zeta:
        return zeta;
    }
    throw std::logic_error("unrecognised structured-grid axis");
}

CartesianPoint NodeCoordinates::position(core::LogicalIndex index) const
{
    return CartesianPoint{x.value(index), y.value(index), z.value(index)};
}

bool NodeCoordinates::is_complete() const
{
    return x.is_bound() && y.is_bound() && z.is_bound();
}

namespace
{

bool directional_views_are_complete(const DirectionalFieldViews &views)
{
    return views.xi.is_bound() && views.eta.is_bound() && views.zeta.is_bound();
}

} // namespace

bool PrimalGeometry::is_complete() const
{
    return cell_volume.is_bound() &&
           directional_views_are_complete(face_area_vector) &&
           directional_views_are_complete(edge_vector) &&
           directional_views_are_complete(face_area) &&
           directional_views_are_complete(edge_length);
}

bool DualGeometry::is_complete() const
{
    return directional_views_are_complete(edge_length_at_primal_face) &&
           directional_views_are_complete(face_area_vector_at_primal_edge) &&
           directional_views_are_complete(face_area_at_primal_edge);
}

bool GeometryBlock::is_complete() const
{
    return spatial_dimension >= 2 && spatial_dimension <= 3 &&
           coordinates.is_complete() && primal.is_complete() && dual.is_complete();
}

GeometrySnapshot::GeometrySnapshot(std::vector<GeometryBlock> blocks)
    : blocks_(std::move(blocks))
{
}

const GeometryBlock &GeometrySnapshot::block(std::size_t index) const
{
    if (index >= blocks_.size())
        throw std::out_of_range("geometry block index lies outside the snapshot");
    return blocks_[index];
}

bool GeometrySnapshot::is_complete() const
{
    if (blocks_.empty())
        return false;
    for (const GeometryBlock &block : blocks_)
    {
        if (!block.is_complete())
            return false;
    }
    return true;
}

} // namespace mpcns::geometry
