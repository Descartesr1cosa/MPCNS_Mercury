#pragma once

#include "mpcns/geometry/geometry_snapshot.hpp"

#include <cstddef>
#include <vector>

namespace mpcns::geometry
{

// Typed hand-off between a geometry producer and the immutable geometry
// backend. The input contains only geometry views: it deliberately carries
// no Field, Grid, MPI, topology, or solver state.
struct GeometryBuildInput
{
    std::vector<GeometryBlock> blocks;
};

// Result of structural geometry finalization. The snapshot still borrows its
// storage during this migration; a future owned backend can preserve this
// result type while changing only the producer.
struct GeometryBuildResult
{
    GeometrySnapshot snapshot;
    std::size_t completed_block_count = 0;
};

// The builder owns the validation boundary between a geometry producer and
// all geometry consumers. It performs no numerical stencil calculation.
class GeometryBuilder
{
public:
    GeometryBuildResult finalize(GeometryBuildInput input) const;
};

} // namespace mpcns::geometry
