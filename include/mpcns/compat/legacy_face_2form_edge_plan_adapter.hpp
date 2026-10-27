#pragma once

#include "mpcns/parallel/face_2form_edge_plan.hpp"

class Halo;
struct FieldHaloRequest;

namespace mpcns::compat
{
class LegacyFace2FormEdgePlanAdapter
{
public:
    // Converts a runtime snapshot captured during the legacy executor.
    static parallel::Face2FormEdgeTransferPlan build(
        const mesh::MeshPartition &partition,
        const Halo &halo,
        const FieldHaloRequest &request);

    // Converts independently prepared initialization-time metadata only.
    static parallel::Face2FormEdgeTransferPlan build_prepared(
        const mesh::MeshPartition &partition,
        const Halo &halo,
        const FieldHaloRequest &request);
};
} // namespace mpcns::compat
