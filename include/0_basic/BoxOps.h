#pragma once

#include "0_basic/TYPES.h"

#include <string>

namespace BOX
{
    Int3 size(const Box3 &b);
    bool empty(const Box3 &b);
    bool nonnegative(const Box3 &b);
    int volume(const Box3 &b);
    Box3 intersect(const Box3 &a, const Box3 &b);
    bool contains(const Box3 &outer, const Box3 &inner);
    bool contains_point(const Box3 &b, const Int3 &p);
    Box3 shifted(const Box3 &b, const Int3 &s);
    void assert_nonnegative(const Box3 &b, const char *where);
    void assert_nonempty(const Box3 &b, const char *where);
    void assert_inside(const Box3 &inner, const Box3 &outer, const char *where);
    std::string to_string(const Box3 &b);
}
