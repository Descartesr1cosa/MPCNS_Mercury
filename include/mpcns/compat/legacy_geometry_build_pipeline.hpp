#pragma once

#include "mpcns/geometry/geometry_builder.hpp"

class Field;
class Grid;

namespace mpcns::compat
{

// All legacy dependencies are confined to this request. The geometry module
// receives only typed GeometryBlock views through GeometryBuildInput.
struct LegacyGeometryBuildRequest
{
    Field &fields;
    Grid &grid;
};

// Coordinates the established publisher order, then hands a validated,
// immutable-in-use snapshot to geometry consumers. It owns no MPI calls and
// retains no Field or Grid state after a call returns.
class LegacyGeometryBuildPipeline
{
public:
    void register_fields(Field &fields, int geometry_ghost) const;
    void publish(const LegacyGeometryBuildRequest &request) const;
    geometry::GeometryBuildResult finalize(
        const LegacyGeometryBuildRequest &request) const;

    geometry::GeometryBuildResult build(
        const LegacyGeometryBuildRequest &request,
        int geometry_ghost) const;
};

} // namespace mpcns::compat
