#pragma once

#include "mpcns/core/axis.hpp"
#include "mpcns/core/logical_index.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace mpcns::geometry
{

// A read-only, non-owning view.  The geometry module deliberately knows
// nothing about Grid, Field, MPI, or a concrete array container.  compat/
// supplies the type-erased read accessor during the migration.
class ConstFieldView
{
public:
    using ReadAccessor = double (*)(const void *storage,
                                    core::LogicalIndex index,
                                    int component);

    ConstFieldView() = default;

    static ConstFieldView borrow(const void *storage,
                                 ReadAccessor read,
                                 core::LogicalBox extent,
                                 int component_count);

    bool is_bound() const;
    const core::LogicalBox &extent() const { return extent_; }
    int component_count() const { return component_count_; }
    double value(core::LogicalIndex index, int component = 0) const;

private:
    const void *storage_ = nullptr;
    ReadAccessor read_ = nullptr;
    core::LogicalBox extent_;
    int component_count_ = 0;
};

struct DirectionalFieldViews
{
    ConstFieldView xi;
    ConstFieldView eta;
    ConstFieldView zeta;

    const ConstFieldView &at(core::Axis axis) const;
};

struct CartesianPoint
{
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct NodeCoordinates
{
    ConstFieldView x;
    ConstFieldView y;
    ConstFieldView z;

    CartesianPoint position(core::LogicalIndex index) const;
    bool is_complete() const;
};

// Metrics of the primal complex.  Vector-valued entries use components
// (x, y, z); scalar entries use component zero.
struct PrimalGeometry
{
    ConstFieldView cell_volume;
    DirectionalFieldViews face_area_vector;
    DirectionalFieldViews edge_vector;
    DirectionalFieldViews face_area;
    DirectionalFieldViews edge_length;

    bool is_complete() const;
};

// Metrics of the dual complex represented on complementary primal locations.
struct DualGeometry
{
    DirectionalFieldViews edge_length_at_primal_face;
    DirectionalFieldViews face_area_vector_at_primal_edge;
    DirectionalFieldViews face_area_at_primal_edge;

    bool is_complete() const;
};

// Diagonal production Hodge coefficients.  They are geometry data and are
// intentionally kept separate from the mutable cochains on which a Hodge
// operation acts.
struct LumpedHodgeGeometry
{
    // |dual edge| / |primal face|, stored at primal faces.
    DirectionalFieldViews star_primal_2_to_dual_1;
    // |primal edge| / |dual face|, stored at primal edges.
    DirectionalFieldViews inverse_star_primal_1;

    bool is_complete() const;
};

// Cached local reconstruction maps.  These are geometric linear-algebra
// data, not Hodge coefficients and not evolving simulation state.
struct ReconstructionCache
{
    DirectionalFieldViews cell_covariant_basis;
    ConstFieldView cell_gradient_pseudoinverse;
    DirectionalFieldViews edge_gradient_pseudoinverse;
    DirectionalFieldViews edge_area_pseudoinverse;

    bool is_complete() const;
};

struct GeometryBlock
{
    std::string name;
    int spatial_dimension = 3;
    NodeCoordinates coordinates;
    PrimalGeometry primal;
    DualGeometry dual;
    LumpedHodgeGeometry lumped_hodge;
    ReconstructionCache reconstruction;

    bool is_complete() const;
};

// Immutable-in-use view of all local geometry blocks.  It does not own the
// storage it views; the compatibility builder documents the required lifetime.
// A future owned geometry backend will retain this public API.
class GeometrySnapshot
{
public:
    GeometrySnapshot() = default;
    explicit GeometrySnapshot(std::vector<GeometryBlock> blocks);

    std::size_t block_count() const { return blocks_.size(); }
    const GeometryBlock &block(std::size_t index) const;
    const std::vector<GeometryBlock> &blocks() const { return blocks_; }
    bool is_complete() const;

private:
    std::vector<GeometryBlock> blocks_;
};

} // namespace mpcns::geometry
