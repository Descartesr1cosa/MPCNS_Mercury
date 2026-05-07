#pragma once

class Field;
class Grid;
class Param;

namespace Z0_NULL
{
    void write_tecplot_output(Param &par,
                              Grid &grid,
                              Field &fields,
                              int step,
                              double time);
}
