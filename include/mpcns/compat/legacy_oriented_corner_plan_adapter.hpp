#pragma once

#include "mpcns/parallel/oriented_corner_plan.hpp"

class Halo;
struct FieldHaloRequest;

namespace mpcns::compat
{
class LegacyOrientedCornerPlanAdapter
{
public:
    static parallel::OrientedCornerTransferPlan build(
        const mesh::MeshPartition &partition,
        const Halo &halo,
        const FieldHaloRequest &request);
};
} // namespace mpcns::compat
