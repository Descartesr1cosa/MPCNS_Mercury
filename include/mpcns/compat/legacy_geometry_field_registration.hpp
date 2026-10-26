#pragma once

#include "0_basic/StaggerLocation.h"

#include <array>
#include <cstddef>

class Field;

namespace mpcns::compat
{

// This is the publishing contract for the Field-backed geometry cache.  It is
// intentionally independent of numerical metric formulas: a future geometry
// builder can publish the same legacy catalog without owning its storage.
enum class LegacyGeometryGhostWidth
{
    Geometry,
    GeometryMinusOne,
    Zero
};

struct LegacyGeometryFieldRegistration
{
    const char *legacy_name = nullptr;
    StaggerLocation location = StaggerLocation::Cell;
    int component_count = 0;
    LegacyGeometryGhostWidth ghost_width = LegacyGeometryGhostWidth::Zero;

    int ghost_layers(int geometry_ghost) const;
};

constexpr std::size_t kLegacyGeometryFieldRegistrationCount = 43;

const std::array<LegacyGeometryFieldRegistration,
                 kLegacyGeometryFieldRegistrationCount> &
legacy_geometry_field_registrations();

// Register the legacy Field catalog in its historical order.  This is a
// temporary publishing adapter; it does not calculate a geometric quantity.
void register_legacy_geometry_fields(Field &fields, int geometry_ghost);

} // namespace mpcns::compat
