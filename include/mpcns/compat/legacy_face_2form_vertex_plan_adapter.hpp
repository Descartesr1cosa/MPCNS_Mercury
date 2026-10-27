#pragma once

#include "mpcns/parallel/face_2form_vertex_plan.hpp"

class Halo;
struct FieldHaloRequest;

namespace mpcns::compat
{
class LegacyFace2FormVertexPlanAdapter
{
public:
    static parallel::Face2FormVertexTransferPlan build(
        const mesh::MeshPartition &partition,
        const Halo &halo,
        const FieldHaloRequest &request);
    static parallel::Face2FormVertexTransferPlan build_prepared(
        const mesh::MeshPartition &partition,
        const Halo &halo,
        const FieldHaloRequest &request);

};
} // namespace mpcns::compat
