#pragma once

#include "Z0_Diagnostics.h"

#include <iosfwd>

class Field;
class Grid;
class Param;

namespace Z0
{
    TestResult test_location_output_smoke(Param &par,
                                          Grid &grid,
                                          Field &fields,
                                          int my_rank,
                                          std::ostream &os);
}
