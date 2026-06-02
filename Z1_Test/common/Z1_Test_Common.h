#pragma once

#include "0_basic/1_MPCNS_Parameter.h"
#include "0_basic/MPI_WRAPPER.h"
#include "1_grid/1_MPCNS_Grid.h"
#include "2_topology/TopologyBuilder.h"
#include "3_field/Field.h"
#include "3_field/FieldDescriptor.h"

#include <cmath>
#include <iostream>
#include <string>

namespace Z1_TEST
{
    struct CaseContext
    {
        int myid = 0;
        Param param;
        Grid grid;
        TOPO::Topology topology;
    };

    inline CaseContext load_case_and_topology()
    {
        CaseContext ctx;
        PARALLEL::mpi_rank(&ctx.myid);
        ctx.param.ReadParam(ctx.myid);
        ctx.grid.Grid_Preprocess(&ctx.param);
        ctx.topology = TOPO::build_topology(ctx.grid, ctx.myid, ctx.param.GetInt("dimension"));
        return ctx;
    }

    inline FieldDescriptor descriptor(const char *name,
                                      StaggerLocation location,
                                      FieldValueKind kind,
                                      int ncomp,
                                      int nghost)
    {
        FieldDescriptor d;
        d.name = name;
        d.location = location;
        d.value_kind = kind;
        d.ncomp = ncomp;
        d.nghost = nghost;
        return d;
    }

    inline void register_dec_test_fields(Field &field, int nghost)
    {
        field.register_field(descriptor("E_xi", StaggerLocation::EdgeXi, FieldValueKind::EdgeCovariant1Form, 1, nghost));
        field.register_field(descriptor("E_eta", StaggerLocation::EdgeEt, FieldValueKind::EdgeCovariant1Form, 1, nghost));
        field.register_field(descriptor("E_zeta", StaggerLocation::EdgeZe, FieldValueKind::EdgeCovariant1Form, 1, nghost));
        field.register_field(descriptor("B_xi", StaggerLocation::FaceXi, FieldValueKind::FaceContravariant2Form, 1, nghost));
        field.register_field(descriptor("B_eta", StaggerLocation::FaceEt, FieldValueKind::FaceContravariant2Form, 1, nghost));
        field.register_field(descriptor("B_zeta", StaggerLocation::FaceZe, FieldValueKind::FaceContravariant2Form, 1, nghost));
        field.register_field(descriptor("divB", StaggerLocation::Cell, FieldValueKind::Scalar, 1, nghost));
    }

    inline void fill_edge_test_form(Field &field)
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
                            fb(i, j, k, 0) = bias + 0.125 * i - 0.25 * j + 0.5 * k;
            };
            fill("E_xi", 1.0);
            fill("E_eta", -2.0);
            fill("E_zeta", 0.75);
        }
    }

    inline bool print_result(int myid, const char *status, const std::string &name, const std::string &detail = "")
    {
        if (myid == 0)
        {
            std::cout << "  [" << status << "] " << name;
            if (!detail.empty())
                std::cout << " " << detail;
            std::cout << "\n";
        }
        return std::string(status) == "PASS" || std::string(status) == "SKIP";
    }

    inline bool print_pass_fail(int myid, const std::string &name, bool passed, const std::string &detail = "")
    {
        return print_result(myid, passed ? "PASS" : "FAIL", name, detail);
    }

    inline bool finite_positive_field(Field &field, const char *name, bool require_positive)
    {
        for (int ib = 0; ib < field.num_blocks(); ++ib)
        {
            FieldBlock &fb = field.field(name, ib);
            if (!fb.is_allocated())
                continue;
            const Int3 lo = fb.get_lo();
            const Int3 hi = fb.get_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const double v = fb(i, j, k, 0);
                        if (!std::isfinite(v))
                            return false;
                        if (require_positive && !(v > 0.0))
                            return false;
                    }
        }
        return true;
    }
}
