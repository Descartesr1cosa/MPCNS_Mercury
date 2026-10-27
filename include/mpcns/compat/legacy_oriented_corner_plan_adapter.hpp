#pragma once

#include "mpcns/parallel/oriented_corner_plan.hpp"

class Halo;
struct FieldHaloRequest;

namespace mpcns::compat
{
class LegacyOrientedCornerPlanAdapter
{
public:
    // Converts a runtime snapshot captured during the legacy executor.
    static parallel::OrientedCornerTransferPlan build(
        const mesh::MeshPartition &partition,
        const Halo &halo,
        const FieldHaloRequest &request);

    // Converts the independently prepared initialization-time metadata.
    // It never triggers prepare, MPI, or a field synchronization.
    static parallel::OrientedCornerTransferPlan build_prepared(
        const mesh::MeshPartition &partition,
        const Halo &halo,
        const FieldHaloRequest &request);
};
} // namespace mpcns::compat
