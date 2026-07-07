#include "Z0_Tests.h"

#include "Z0_Boundary.h"
#include "Z0_TestCommon.h"
#include "8_dec/DecOps.h"

#include <cmath>
#include <iostream>
#include <sstream>

namespace
{
    void fill_phi(Field &field)
    {
        for (int ib = 0; ib < field.num_blocks(); ++ib)
        {
            FieldBlock &phi = field.field("phi", ib);
            if (!phi.is_allocated())
                continue;
            const Int3 lo = phi.get_lo();
            const Int3 hi = phi.get_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        phi(i, j, k, 0) =
                            0.125 * i * i - 0.25 * j * k + 0.5 * k + 3.0;
        }
    }

    void fill_edge_form(Field &field)
    {
        for (int ib = 0; ib < field.num_blocks(); ++ib)
        {
            auto fill = [&](const char *name, double bias)
            {
                FieldBlock &fb = field.field(name, ib);
                if (!fb.is_allocated())
                    return;
                const Int3 lo = fb.get_lo();
                const Int3 hi = fb.get_hi();
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                            fb(i, j, k, 0) =
                                bias + 0.375 * i - 0.5 * j + 0.625 * k +
                                0.0625 * i * j - 0.03125 * j * k;
            };
            fill("E_xi", 1.0);
            fill("E_eta", -2.0);
            fill("E_zeta", 0.75);
        }
    }

    double max_abs_owned(Field &field, const std::string &name)
    {
        double max_abs = 0.0;
        const int fid = field.field_id(name);
        const FieldDescriptor &desc = field.descriptor(fid);
        for (int ib = 0; ib < field.num_blocks(); ++ib)
        {
            FieldBlock &fb = field.field(fid, ib);
            if (!fb.is_allocated())
                continue;
            const Box3 owned = Z0_TEST::owned_box(field, ib, desc.location);
            for (int i = owned.lo.i; i < owned.hi.i; ++i)
                for (int j = owned.lo.j; j < owned.hi.j; ++j)
                    for (int k = owned.lo.k; k < owned.hi.k; ++k)
                        for (int m = 0; m < desc.ncomp; ++m)
                            max_abs = std::max(max_abs, std::abs(fb(i, j, k, m)));
        }
        return Z0_TEST::global_max(max_abs);
    }
}

namespace Z0_TEST
{
    bool RunDecChainTests(Field &field, Halo &halo, Z0_Boundary &boundary, Param &param)
    {
        (void)halo;
        (void)param;
        if (rank() == 0)
            std::cout << "Z0 DEC chain tests\n";

        bool passed = true;
        const DEC::EdgeFormNames e{"E_xi", "E_eta", "E_zeta"};
        const DEC::FaceFormNames b{"B_xi", "B_eta", "B_zeta"};

        fill_all_with_nan(field);
        fill_phi(field);
        DEC::d0_node_to_edge(field, "phi", e);
        DEC::d1_edge_to_face(field, e, b);
        const double d1d0 = std::max(max_abs_owned(field, "B_xi"),
                                     std::max(max_abs_owned(field, "B_eta"),
                                              max_abs_owned(field, "B_zeta")));
        passed &= print_result("D1 * D0 == 0", d1d0 < 1.0e-12, max_detail(d1d0));

        fill_all_with_nan(field);
        fill_edge_form(field);
        DEC::d1_edge_to_face(field, e, b);
        DEC::d2_face_to_cell(field, b, "divB");
        const double d2d1 = max_abs_owned(field, "divB");
        passed &= print_result("D2 * D1 == 0", d2d1 < 1.0e-12, max_detail(d2d1));

        return passed;
    }
}
