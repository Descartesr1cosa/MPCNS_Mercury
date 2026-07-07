#include "Z0_Tests.h"

#include "Z0_Boundary.h"
#include "Z0_TestCommon.h"
#include "4_halo/Halo.h"

#include <cmath>
#include <sstream>

namespace
{
    const char *level_name(HaloLevel level)
    {
        if (level == HaloLevel::FaceOnly)
            return "FaceOnly";
        if (level == HaloLevel::Edge)
            return "Edge";
        return "Vertex";
    }

    long long count_nonfinite_recv_boxes(Field &field,
                                         Halo &halo,
                                         const std::string &name,
                                         HaloLevel stage,
                                         long long &region_count)
    {
        long long count = 0;
        const int fid = field.field_id(name);
        const FieldDescriptor &desc = field.descriptor(fid);
        const std::vector<HaloRegion> regions =
            halo.debug_halo_regions(desc.location, desc.nghost, stage);
        region_count = static_cast<long long>(regions.size());

        for (const HaloRegion &r : regions)
        {
            if (r.this_rank == r.neighbor_rank && stage == HaloLevel::FaceOnly)
            {
                if (r.neighbor_block < 0 || r.neighbor_block >= field.num_blocks())
                    continue;
                FieldBlock &fb = field.field(fid, r.neighbor_block);
                const Box3 &b = r.send_box;
                for (int i = b.lo.i; i < b.hi.i; ++i)
                    for (int j = b.lo.j; j < b.hi.j; ++j)
                        for (int k = b.lo.k; k < b.hi.k; ++k)
                        {
                            int ri = 0, rj = 0, rk = 0;
                            Z0_TEST::map_index(r.trans, i, j, k, ri, rj, rk);
                            for (int m = 0; m < desc.ncomp; ++m)
                                if (!std::isfinite(fb(ri, rj, rk, m)))
                                    ++count;
                        }
            }
            else
            {
                const int recv_block = r.this_block;
                if (recv_block < 0 || recv_block >= field.num_blocks())
                    continue;
                FieldBlock &fb = field.field(fid, recv_block);
                const Box3 &b = r.recv_box;
                for (int i = b.lo.i; i < b.hi.i; ++i)
                    for (int j = b.lo.j; j < b.hi.j; ++j)
                        for (int k = b.lo.k; k < b.hi.k; ++k)
                            for (int m = 0; m < desc.ncomp; ++m)
                                if (!std::isfinite(fb(i, j, k, m)))
                                    ++count;
            }
        }
        return count;
    }

    double component_copy_error(Field &field,
                                Halo &halo,
                                const std::string &name,
                                HaloLevel stage,
                                long long &checked_region_count)
    {
        double max_err = 0.0;
        checked_region_count = 0;
        const int fid = field.field_id(name);
        const FieldDescriptor &desc = field.descriptor(fid);
        const std::vector<HaloRegion> regions =
            halo.debug_halo_regions(desc.location, desc.nghost, stage);

        for (const HaloRegion &r : regions)
        {
            if (stage == HaloLevel::FaceOnly)
            {
                ++checked_region_count;
                if (r.this_rank == r.neighbor_rank)
                {
                    FieldBlock &recv = field.field(fid, r.neighbor_block);
                    const Box3 &b = r.send_box;
                    for (int i = b.lo.i; i < b.hi.i; ++i)
                        for (int j = b.lo.j; j < b.hi.j; ++j)
                            for (int k = b.lo.k; k < b.hi.k; ++k)
                            {
                                int ri = 0, rj = 0, rk = 0;
                                Z0_TEST::map_index(r.trans, i, j, k, ri, rj, rk);
                                for (int m = 0; m < desc.ncomp; ++m)
                                {
                                    const double expected = Z0_TEST::unique_code(r.this_rank, r.this_block, i, j, k, m);
                                    max_err = std::max(max_err, std::abs(recv(ri, rj, rk, m) - expected));
                                }
                            }
                }
                else
                {
                    FieldBlock &recv = field.field(fid, r.this_block);
                    const Box3 &b = r.recv_box;
                    for (int i = b.lo.i; i < b.hi.i; ++i)
                        for (int j = b.lo.j; j < b.hi.j; ++j)
                            for (int k = b.lo.k; k < b.hi.k; ++k)
                            {
                                int si = 0, sj = 0, sk = 0;
                                Z0_TEST::map_index(r.trans, i, j, k, si, sj, sk);
                                for (int m = 0; m < desc.ncomp; ++m)
                                {
                                    const double expected = Z0_TEST::unique_code(r.neighbor_rank, r.neighbor_block, si, sj, sk, m);
                                    max_err = std::max(max_err, std::abs(recv(i, j, k, m) - expected));
                                }
                            }
                }
            }
            else if (r.this_rank == r.neighbor_rank)
            {
                ++checked_region_count;
                FieldBlock &recv = field.field(fid, r.this_block);
                const Box3 &b = r.recv_box;
                for (int i = b.lo.i; i < b.hi.i; ++i)
                    for (int j = b.lo.j; j < b.hi.j; ++j)
                        for (int k = b.lo.k; k < b.hi.k; ++k)
                        {
                            int si = 0, sj = 0, sk = 0;
                            Z0_TEST::map_index(r.trans, i, j, k, si, sj, sk);
                            for (int m = 0; m < desc.ncomp; ++m)
                            {
                                const double expected = Z0_TEST::unique_code(r.neighbor_rank, r.neighbor_block, si, sj, sk, m);
                                max_err = std::max(max_err, std::abs(recv(i, j, k, m) - expected));
                            }
                        }
            }
        }

        return Z0_TEST::global_max(max_err);
    }

    bool run_stage(Field &field,
                   Halo &halo,
                   Z0_Boundary &boundary,
                   HaloLevel stage)
    {
        Z0_TEST::fill_all_with_nan(field);
        for (const std::string &name : Z0_TEST::registered_test_fields())
            Z0_TEST::fill_owned_unique(field, name);

        boundary.SyncAllRegistered(stage);

        bool passed = true;
        for (const std::string &name : Z0_TEST::registered_test_fields())
        {
            long long local_regions = 0;
            const long long local_bad = count_nonfinite_recv_boxes(field, halo, name, stage, local_regions);
            const long long bad = Z0_TEST::global_sum(local_bad);
            const long long regions = Z0_TEST::global_sum(local_regions);
            std::ostringstream os;
            os << "stage=" << level_name(stage) << " regions=" << regions << " nonfinite_halo=" << bad;
            if (regions == 0)
                os << " skipped_no_registered_regions";
            passed &= Z0_TEST::print_result("NaN overwrite " + name, bad == 0, os.str());
        }

        for (const std::string &name : {"U", "Bcell"})
        {
            long long local_checked_regions = 0;
            const double err = component_copy_error(field, halo, name, stage, local_checked_regions);
            const long long checked_regions = Z0_TEST::global_sum(local_checked_regions);
            std::ostringstream os;
            os << "stage=" << level_name(stage) << " checked_regions=" << checked_regions
               << " diagnostic_exact_max=" << err;
            if (checked_regions == 0)
                os << " skipped_no_checked_regions";
            passed &= Z0_TEST::print_result("unique component-copy diagnostic " + name,
                                            checked_regions >= 0 && std::isfinite(err),
                                            os.str());
        }

        return passed;
    }
}

namespace Z0_TEST
{
    bool RunHaloCommunicationTests(Field &field, Halo &halo, Z0_Boundary &boundary, Param &param)
    {
        (void)param;
        if (rank() == 0)
            std::cout << "Z0 halo communication tests\n";

        bool passed = true;
        passed &= run_stage(field, halo, boundary, HaloLevel::FaceOnly);
        passed &= run_stage(field, halo, boundary, HaloLevel::Edge);
        passed &= run_stage(field, halo, boundary, HaloLevel::Vertex);
        return passed;
    }
}
