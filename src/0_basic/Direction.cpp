#include "0_basic/Direction.h"

#include "0_basic/Error.h"

namespace
{
    int abs_int(int v)
    {
        return v < 0 ? -v : v;
    }

    void assert_valid_dir(DIR::Code d, const char *where)
    {
        if (!DIR::is_valid(d))
            ERROR::Abort(where);
    }
}

namespace DIR
{
    bool is_valid(Code d)
    {
        const int a = abs_int(d);
        return a >= 1 && a <= 3;
    }

    int axis(Code d)
    {
        assert_valid_dir(d, "DIR::axis: invalid direction code");
        return abs_int(d) - 1;
    }

    int axis1(Code d)
    {
        assert_valid_dir(d, "DIR::axis1: invalid direction code");
        return abs_int(d);
    }

    int sign(Code d)
    {
        assert_valid_dir(d, "DIR::sign: invalid direction code");
        return d < 0 ? -1 : +1;
    }

    Code make(int axis0, int sgn)
    {
        if (axis0 < 0 || axis0 > 2)
            ERROR::Abort("DIR::make: invalid axis");
        if (sgn != -1 && sgn != +1)
            ERROR::Abort("DIR::make: invalid sign");
        return sgn * (axis0 + 1);
    }

    Code opposite(Code d)
    {
        assert_valid_dir(d, "DIR::opposite: invalid direction code");
        return -d;
    }

    bool same_axis(Code a, Code b)
    {
        return axis(a) == axis(b);
    }

    bool distinct_axes(Code a, Code b)
    {
        return !same_axis(a, b);
    }

    bool distinct_axes(Code a, Code b, Code c)
    {
        const int ax_a = axis(a);
        const int ax_b = axis(b);
        const int ax_c = axis(c);
        return ax_a != ax_b && ax_a != ax_c && ax_b != ax_c;
    }

    const char *name(Code d)
    {
        switch (d)
        {
        case -1:
            return "XMinus";
        case +1:
            return "XPlus";
        case -2:
            return "YMinus";
        case +2:
            return "YPlus";
        case -3:
            return "ZMinus";
        case +3:
            return "ZPlus";
        default:
            return "Invalid";
        }
    }
}
