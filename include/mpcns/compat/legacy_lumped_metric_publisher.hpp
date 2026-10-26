#pragma once

class Field;
class Grid;

namespace mpcns::compat
{

// Publishes the legacy Field-backed face/edge measures and diagonal Hodge
// coefficients from existing Grid coordinates. This explicit compatibility
// boundary owns all Field and Grid access; pure geometry and Hodge policies
// remain in their respective bottom-layer modules.
void populate_legacy_lumped_face_and_edge_metrics(Field &fields, Grid &grid);

} // namespace mpcns::compat
