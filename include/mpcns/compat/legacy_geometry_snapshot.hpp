#pragma once

#include "mpcns/geometry/geometry_snapshot.hpp"

class Field;
class Grid;

namespace mpcns::compat
{

// Create a read-only geometry snapshot over the exact existing Grid and
// Field-backed metric data.  This performs no MPI, allocation, interpolation,
// or numerical geometry calculation.  The caller must keep grid and fields
// alive and must not mutate their geometry storage while using the snapshot.
geometry::GeometrySnapshot borrow_legacy_geometry_snapshot(Grid &grid,
                                                            Field &fields);

} // namespace mpcns::compat
