#pragma once

#include "mpcns/parallel/face_2form_edge_plan.hpp"

class Halo;
struct FieldHaloRequest;

namespace mpcns::compat
{
class LegacyFace2FormEdgePlanAdapter
{
public:
    static parallel::Face2FormEdgeTransferPlan build(
        const mesh::MeshPartition &partition,
        const Halo &halo,
        const FieldHaloRequest &request);
};
} // namespace mpcns::compat
