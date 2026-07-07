#include "Z0_Tests.h"

#include "Z0_TestCommon.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <sstream>

namespace
{
    std::string io_file_name(int myid)
    {
        std::ostringstream os;
        os << "z0_physical_io_rank" << myid << ".txt";
        return os.str();
    }

    void initialize_io_fields(Field &field)
    {
        Z0_TEST::fill_all_with_nan(field);
        for (const std::string &name : Z0_TEST::registered_test_fields())
            Z0_TEST::fill_owned_unique(field, name);
    }

    bool write_owned_fields(Field &field, const std::string &path)
    {
        std::ofstream out(path.c_str());
        if (!out)
            return false;

        out.precision(17);
        for (const std::string &name : Z0_TEST::registered_test_fields())
        {
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
                                out << name << ' ' << ib << ' ' << i << ' ' << j << ' '
                                    << k << ' ' << m << ' ' << fb(i, j, k, m) << '\n';
            }
        }

        return static_cast<bool>(out);
    }

    double read_and_check(Field &field, const std::string &path, long long &records)
    {
        std::ifstream in(path.c_str());
        if (!in)
            return std::numeric_limits<double>::infinity();

        std::string name;
        int ib = 0, i = 0, j = 0, k = 0, m = 0;
        double value = 0.0;
        double max_err = 0.0;
        records = 0;
        while (in >> name >> ib >> i >> j >> k >> m >> value)
        {
            const double got = field.field(name, ib)(i, j, k, m);
            max_err = std::max(max_err, std::abs(got - value));
            ++records;
        }

        return max_err;
    }
}

namespace Z0_TEST
{
    bool RunPhysicalIoTests(Field &field, Param &param)
    {
        (void)param;
        if (rank() == 0)
            std::cout << "Z0 physical field IO tests\n";

        initialize_io_fields(field);

        const std::string path = io_file_name(rank());
        const bool wrote = write_owned_fields(field, path);
        bool passed = print_result("write owned physical fields", wrote, path);

        long long local_records = 0;
        const double local_err = read_and_check(field, path, local_records);
        const double err = global_max(local_err);
        const long long records = global_sum(local_records);

        std::ostringstream os;
        os << "records=" << records << " max_roundtrip_error=" << err;
        passed &= print_result("read owned physical fields", records > 0 && err == 0.0, os.str());

        std::remove(path.c_str());
        return passed;
    }
}
