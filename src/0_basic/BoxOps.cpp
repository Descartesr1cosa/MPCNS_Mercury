#include "0_basic/BoxOps.h"

#include "0_basic/Error.h"

namespace
{
    int min_int(int a, int b)
    {
        return a < b ? a : b;
    }

    int max_int(int a, int b)
    {
        return a > b ? a : b;
    }

    std::string where_or_default(const char *where, const char *fallback)
    {
        return where ? where : fallback;
    }
}

namespace BOX
{
    Int3 size(const Box3 &b)
    {
        return {b.hi.i - b.lo.i,
                b.hi.j - b.lo.j,
                b.hi.k - b.lo.k};
    }

    bool empty(const Box3 &b)
    {
        const Int3 s = size(b);
        return s.i <= 0 || s.j <= 0 || s.k <= 0;
    }

    bool nonnegative(const Box3 &b)
    {
        const Int3 s = size(b);
        return s.i >= 0 && s.j >= 0 && s.k >= 0;
    }

    int volume(const Box3 &b)
    {
        if (empty(b))
            return 0;

        const Int3 s = size(b);
        return s.i * s.j * s.k;
    }

    Box3 intersect(const Box3 &a, const Box3 &b)
    {
        return Box3{{max_int(a.lo.i, b.lo.i),
                     max_int(a.lo.j, b.lo.j),
                     max_int(a.lo.k, b.lo.k)},
                    {min_int(a.hi.i, b.hi.i),
                     min_int(a.hi.j, b.hi.j),
                     min_int(a.hi.k, b.hi.k)}};
    }

    bool contains(const Box3 &outer, const Box3 &inner)
    {
        return inner.lo.i >= outer.lo.i &&
               inner.lo.j >= outer.lo.j &&
               inner.lo.k >= outer.lo.k &&
               inner.hi.i <= outer.hi.i &&
               inner.hi.j <= outer.hi.j &&
               inner.hi.k <= outer.hi.k;
    }

    bool contains_point(const Box3 &b, const Int3 &p)
    {
        return p.i >= b.lo.i && p.i < b.hi.i &&
               p.j >= b.lo.j && p.j < b.hi.j &&
               p.k >= b.lo.k && p.k < b.hi.k;
    }

    Box3 shifted(const Box3 &b, const Int3 &s)
    {
        return Box3{{b.lo.i + s.i,
                     b.lo.j + s.j,
                     b.lo.k + s.k},
                    {b.hi.i + s.i,
                     b.hi.j + s.j,
                     b.hi.k + s.k}};
    }

    void assert_nonnegative(const Box3 &b, const char *where)
    {
        if (!nonnegative(b))
            ERROR::Abort(where_or_default(where, "BOX::assert_nonnegative"));
    }

    void assert_nonempty(const Box3 &b, const char *where)
    {
        if (empty(b))
            ERROR::Abort(where_or_default(where, "BOX::assert_nonempty"));
    }

    void assert_inside(const Box3 &inner, const Box3 &outer, const char *where)
    {
        if (!contains(outer, inner))
            ERROR::Abort(where_or_default(where, "BOX::assert_inside"));
    }

    std::string to_string(const Box3 &b)
    {
        return "[(" + std::to_string(b.lo.i) + "," +
               std::to_string(b.lo.j) + "," +
               std::to_string(b.lo.k) + "),(" +
               std::to_string(b.hi.i) + "," +
               std::to_string(b.hi.j) + "," +
               std::to_string(b.hi.k) + "))";
    }
}
