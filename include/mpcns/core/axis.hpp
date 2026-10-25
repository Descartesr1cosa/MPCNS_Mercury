#pragma once

#include <cstdint>

namespace mpcns::core
{

// Logical directions of the structured primal complex.  New APIs use this
// type instead of integer direction codes or legacy Et/Ze spellings.
enum class Axis : std::uint8_t
{
    Xi,
    Eta,
    Zeta
};

} // namespace mpcns::core
