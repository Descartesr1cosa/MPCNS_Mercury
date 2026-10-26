#pragma once

class Field;
class Grid;

namespace mpcns::compat
{

// Publishes historical primal face-area vectors and oriented cell volumes.
// Call before every publisher that consumes JDxi, JDet, JDze, or Jac.
void populate_legacy_primal_face_and_cell_metrics(Field &fields, Grid &grid);

// Publishes the three historical primal edge-vector fields from node coordinates.
// Kept separate so its legacy execution position remains explicit.
void populate_legacy_primal_edge_vectors(Field &fields, Grid &grid);

} // namespace mpcns::compat
