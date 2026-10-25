#include "Z0_Tests.h"

#include "Z0_Boundary.h"
#include "Z0_TestCommon.h"
#include "8_dec/DecOps.h"
#include "mpcns/discrete/incidence_operator.hpp"

#include <algorithm>
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

    // Regression oracle for the production CT curl stencil as it existed
    // before the d1 implementation was centralized.  Keep this formula in
    // the test only: numerical code must use mpcns::discrete::apply_d1.
    double max_difference_from_legacy_ct_d1(Field &field,
                                            const DEC::EdgeFormNames &edge_names,
                                            const DEC::FaceFormNames &face_names)
    {
        double max_difference = 0.0;
        for (int ib = 0; ib < field.num_blocks(); ++ib)
        {
            FieldBlock &edge_xi = field.field(edge_names.xi, ib);
            FieldBlock &edge_eta = field.field(edge_names.eta, ib);
            FieldBlock &edge_zeta = field.field(edge_names.zeta, ib);
            FieldBlock &face_xi = field.field(face_names.xi, ib);
            FieldBlock &face_eta = field.field(face_names.eta, ib);
            FieldBlock &face_zeta = field.field(face_names.zeta, ib);

            if (!edge_xi.is_allocated() || !edge_eta.is_allocated() ||
                !edge_zeta.is_allocated() || !face_xi.is_allocated() ||
                !face_eta.is_allocated() || !face_zeta.is_allocated())
                continue;

            const auto update_max = [&](const FieldBlock &face, int i, int j, int k, double expected)
            {
                max_difference = std::max(
                    max_difference,
                    std::abs(face(i, j, k, 0) - expected));
            };

            {
                const Int3 lo = face_xi.inner_lo();
                const Int3 hi = face_xi.inner_hi();
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                            update_max(
                                face_xi, i, j, k,
                                (edge_eta(i, j, k, 0) - edge_eta(i, j, k + 1, 0)) +
                                    (edge_zeta(i, j + 1, k, 0) - edge_zeta(i, j, k, 0)));
            }

            {
                const Int3 lo = face_eta.inner_lo();
                const Int3 hi = face_eta.inner_hi();
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                            update_max(
                                face_eta, i, j, k,
                                (edge_xi(i, j, k + 1, 0) - edge_xi(i, j, k, 0)) +
                                    (edge_zeta(i, j, k, 0) - edge_zeta(i + 1, j, k, 0)));
            }

            {
                const Int3 lo = face_zeta.inner_lo();
                const Int3 hi = face_zeta.inner_hi();
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                            update_max(
                                face_zeta, i, j, k,
                                (edge_xi(i, j, k, 0) - edge_xi(i, j + 1, k, 0)) +
                                    (edge_eta(i + 1, j, k, 0) - edge_eta(i, j, k, 0)));
            }
        }
        return Z0_TEST::global_max(max_difference);
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
        for (int ib = 0; ib < field.num_blocks(); ++ib)
        {
            auto &edge_xi = field.field(e.xi, ib);
            auto &edge_eta = field.field(e.eta, ib);
            auto &edge_zeta = field.field(e.zeta, ib);
            auto &face_xi = field.field(b.xi, ib);
            auto &face_eta = field.field(b.eta, ib);
            auto &face_zeta = field.field(b.zeta, ib);

            if (!edge_xi.is_allocated() || !edge_eta.is_allocated() ||
                !edge_zeta.is_allocated() || !face_xi.is_allocated() ||
                !face_eta.is_allocated() || !face_zeta.is_allocated())
                continue;

            const mpcns::state::PrimalEdgeFormView edge{
                edge_xi, edge_eta, edge_zeta};
            const mpcns::state::PrimalFaceFormView face{
                face_xi, face_eta, face_zeta};
            mpcns::discrete::apply_d1(edge, face);
        }
        const double d1_legacy_difference =
            max_difference_from_legacy_ct_d1(field, e, b);
        passed &= print_result("MPCNS d1 matches legacy CT stencil exactly",
                               d1_legacy_difference == 0.0,
                               max_detail(d1_legacy_difference));
        DEC::d2_face_to_cell(field, b, "divB");
        const double d2d1 = max_abs_owned(field, "divB");
        passed &= print_result("MPCNS D2 * D1 == 0", d2d1 < 1.0e-12, max_detail(d2d1));

        return passed;
    }
}
