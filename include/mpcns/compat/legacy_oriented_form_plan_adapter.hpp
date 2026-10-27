#pragma once

#include "mpcns/parallel/oriented_form_plan.hpp"

#include <string>

class Halo;
struct FieldHaloRequest;

namespace mpcns::compat
{
class LegacyOrientedFormPlanAdapter
{
public:
    static parallel::OrientedFormTransferPlan build(
        const mesh::MeshPartition &partition,
        const Halo &halo,
        const FieldHaloRequest &request);

private:
    static parallel::OrientedFormTransferPlan build_group(
        const mesh::MeshPartition &partition,
        const Halo &halo,
        const std::string &sync_group);
};
} // namespace mpcns::compat
