#pragma once

#include "0_basic/1_MPCNS_Parameter.h"
#include "0_basic/LayoutTraits.h"
#include "0_basic/MPI_WRAPPER.h"
#include "3_field/Field.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace Z0_TEST
{
    inline int rank()
    {
        int myid = 0;
        PARALLEL::mpi_rank(&myid);
        return myid;
    }

    inline bool print_result(const std::string &name, bool passed, const std::string &detail)
    {
        if (rank() == 0)
        {
            std::cout << "  [" << (passed ? "PASS" : "FAIL") << "] " << name;
            if (!detail.empty())
                std::cout << " " << detail;
            std::cout << "\n";
        }
        return passed;
    }

    inline long long global_sum(long long local)
    {
        double in = static_cast<double>(local);
        double out = 0.0;
        PARALLEL::mpi_sum(&in, &out, 1);
        return static_cast<long long>(out);
    }

    inline double global_max(double local)
    {
        double out = 0.0;
        PARALLEL::mpi_max(&local, &out, 1);
        return out;
    }

    inline Box3 owned_box(Field &field, int ib, StaggerLocation loc)
    {
        const Block &blk = field.grd->grids(ib);
        return LAYOUT::owned_box_from_cells({blk.mx, blk.my, blk.mz}, loc);
    }

    inline bool in_box(const Box3 &b, int i, int j, int k)
    {
        return i >= b.lo.i && i < b.hi.i &&
               j >= b.lo.j && j < b.hi.j &&
               k >= b.lo.k && k < b.hi.k;
    }

    inline std::vector<std::string> registered_test_fields()
    {
        return {"phi", "E_xi", "E_eta", "E_zeta",
                "B_xi", "B_eta", "B_zeta",
                "divB", "U", "Bcell"};
    }

    inline void fill_all_with_nan(Field &field)
    {
        const double qnan = std::numeric_limits<double>::quiet_NaN();
        for (const FieldDescriptor &desc : field.descriptors())
        {
            const int fid = field.field_id(desc.name);
            for (int ib = 0; ib < field.num_blocks(); ++ib)
            {
                FieldBlock &fb = field.field(fid, ib);
                if (!fb.is_allocated())
                    continue;
                const Int3 lo = fb.get_lo();
                const Int3 hi = fb.get_hi();
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                            for (int m = 0; m < desc.ncomp; ++m)
                                fb(i, j, k, m) = qnan;
            }
        }
    }

    inline double unique_code(int myid, int ib, int i, int j, int k, int m)
    {
        return 1.0e9 * static_cast<double>(myid + 1) +
               1.0e6 * static_cast<double>(ib + 1) +
               1.0e3 * static_cast<double>(m + 1) +
               100.0 * static_cast<double>(i) +
               10.0 * static_cast<double>(j) +
               static_cast<double>(k);
    }

    inline void fill_owned_unique(Field &field, const std::string &name)
    {
        const int fid = field.field_id(name);
        const FieldDescriptor &desc = field.descriptor(fid);
        const int myid = rank();
        for (int ib = 0; ib < field.num_blocks(); ++ib)
        {
            FieldBlock &fb = field.field(fid, ib);
            if (!fb.is_allocated())
                continue;
            const Box3 owned = owned_box(field, ib, desc.location);
            for (int i = owned.lo.i; i < owned.hi.i; ++i)
                for (int j = owned.lo.j; j < owned.hi.j; ++j)
                    for (int k = owned.lo.k; k < owned.hi.k; ++k)
                        for (int m = 0; m < desc.ncomp; ++m)
                            fb(i, j, k, m) = unique_code(myid, ib, i, j, k, m);
        }
    }

    inline long long count_nonfinite_owned(Field &field, const std::string &name)
    {
        long long count = 0;
        const int fid = field.field_id(name);
        const FieldDescriptor &desc = field.descriptor(fid);
        for (int ib = 0; ib < field.num_blocks(); ++ib)
        {
            FieldBlock &fb = field.field(fid, ib);
            if (!fb.is_allocated())
                continue;
            const Box3 owned = owned_box(field, ib, desc.location);
            for (int i = owned.lo.i; i < owned.hi.i; ++i)
                for (int j = owned.lo.j; j < owned.hi.j; ++j)
                    for (int k = owned.lo.k; k < owned.hi.k; ++k)
                        for (int m = 0; m < desc.ncomp; ++m)
                            if (!std::isfinite(fb(i, j, k, m)))
                                ++count;
        }
        return count;
    }

    inline void map_index(const TOPO::IndexTransform &tr,
                          int i, int j, int k,
                          int &io, int &jo, int &ko)
    {
        const int src[3] = {i, j, k};
        int dst[3] = {0, 0, 0};
        const int offset[3] = {tr.offset.i, tr.offset.j, tr.offset.k};
        for (int axis = 0; axis < 3; ++axis)
            dst[tr.perm[axis]] = tr.sign[axis] * src[axis] + offset[axis];
        io = dst[0];
        jo = dst[1];
        ko = dst[2];
    }

    inline std::string max_detail(double value)
    {
        return "max=" + std::to_string(value);
    }
}
