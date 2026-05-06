#include "2_topology/2_MPCNS_Topology.h"

void TOPO::node_box_from_subsup(const int sub[3], const int sup[3], Box3 &box)
{
    for (int d = 0; d < 3; ++d)
    {
        int s = std::abs(sub[d]);
        int t = std::abs(sup[d]);
        int a = std::min(s, t);
        int b = std::max(s, t);

        if (a == b)
        {
            // 单层：厚度 1
            if (d == 0)
            {
                box.lo.i = a;
                box.hi.i = a + 1;
            }
            if (d == 1)
            {
                box.lo.j = a;
                box.hi.j = a + 1;
            }
            if (d == 2)
            {
                box.lo.k = a;
                box.hi.k = a + 1;
            }
        }
        else
        {
            // 区间 [a,b] → [a,b+1)
            if (d == 0)
            {
                box.lo.i = a;
                box.hi.i = b + 1;
            }
            if (d == 1)
            {
                box.lo.j = a;
                box.hi.j = b + 1;
            }
            if (d == 2)
            {
                box.lo.k = a;
                box.hi.k = b + 1;
            }
        }
    }
}