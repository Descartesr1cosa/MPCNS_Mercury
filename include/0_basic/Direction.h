#pragma once

namespace DIR
{
    using Code = int;

    bool is_valid(Code d);
    int axis(Code d);
    int axis1(Code d);
    int sign(Code d);
    Code make(int axis0, int sign);
    Code opposite(Code d);
    bool same_axis(Code a, Code b);
    bool distinct_axes(Code a, Code b);
    bool distinct_axes(Code a, Code b, Code c);
    const char *name(Code d);
}
