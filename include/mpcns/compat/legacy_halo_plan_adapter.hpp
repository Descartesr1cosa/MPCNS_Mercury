#pragma once

#include "mpcns/parallel/exchange_plan.hpp"

class Halo;
struct FieldHaloRequest;

namespace mpcns::compat
{
class LegacyHaloPlanAdapter
{
public:
    static parallel::ExchangePlan build(const mesh::MeshPartition &partition,
                                        const Halo &halo,
                                        const FieldHaloRequest &request);
};
} // namespace mpcns::compat
