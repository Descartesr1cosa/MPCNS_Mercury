#pragma once

#include "mpcns/parallel/oriented_vertex_plan.hpp"

class Halo;
struct FieldHaloRequest;

namespace mpcns::compat
{
class LegacyOrientedVertexPlanAdapter
{
public:
    // Converts a runtime snapshot captured during the legacy executor.
    static parallel::OrientedVertexTransferPlan build(
        const mesh::MeshPartition &partition,
        const Halo &halo,
        const FieldHaloRequest &request);

    // Converts independently prepared initialization-time metadata only.
    static parallel::OrientedVertexTransferPlan build_prepared(
        const mesh::MeshPartition &partition,
        const Halo &halo,
        const FieldHaloRequest &request);
};
} // namespace mpcns::compat
