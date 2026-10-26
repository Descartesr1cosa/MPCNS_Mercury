#pragma once

class Field;
class Grid;

namespace mpcns::compat
{

// Publishes the historical covariant-basis and cell/edge pseudo-inverse caches.
// Call after primal face vectors are available; this adapter owns no MPI or Hodge fields.
void populate_legacy_reconstruction_cache(Field &fields, Grid &grid);
// Publishes diagonal star 0<->3 and forward star 1->2 from current metric fields.
// Call after the inverse edge Hodge values have been published.
void populate_legacy_diagonal_hodge_closure(Field &fields, Grid &grid);

} // namespace mpcns::compat
