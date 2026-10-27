#pragma once

#include "mpcns/parallel/oriented_vertex_plan.hpp"

class Halo;
struct FieldHaloRequest;

namespace mpcns::compat
{
class LegacyOrientedVertexPlanAdapter
{
public:
    static parallel::OrientedVertexTransferPlan build(
        const mesh::MeshPartition &partition,
        const Halo &halo,
        const FieldHaloRequest &request);
};
} // namespace mpcns::compat
