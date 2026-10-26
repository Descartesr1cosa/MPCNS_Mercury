#pragma once

class Field;

namespace METRIC
{
struct SingularPhysicalEdge;
}

namespace mpcns::compat
{

// Legacy Field-catalog publication for a singular-edge geometry closure.
// Topology, MPI reduction, and geometry closure remain outside this adapter.
void publish_legacy_singular_edge_hodge(
    Field &fields,
    const METRIC::SingularPhysicalEdge &edge,
    int local_rank);

} // namespace mpcns::compat
