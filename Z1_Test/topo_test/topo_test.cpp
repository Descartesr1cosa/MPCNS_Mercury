#include "0_basic/1_MPCNS_Parameter.h"
#include "0_basic/MPI_WRAPPER.h"
#include "1_grid/1_MPCNS_Grid.h"
#include "2_topology/GlobalIncidence.h"
#include "2_topology/LocalIncidence.h"
#include "2_topology/TopologyBuilder.h"

#include <iostream>
#include <stdexcept>
#include <string>

namespace
{
    struct TestSummary
    {
        long long checked = 0;
        long long failed = 0;
    };

    TestSummary check_local_boundary_boundary_face(const TOPO::Topology &topology)
    {
        TestSummary result;
        for (const auto &entry : topology.face2key)
        {
            ++result.checked;
            if (!TOPO::check_boundary_boundary_face_zero(entry.first))
                ++result.failed;
        }
        return result;
    }

    TestSummary check_local_boundary_boundary_cell(const TOPO::Topology &topology)
    {
        TestSummary result;
        for (const auto &entry : topology.cell_to_id)
        {
            ++result.checked;
            if (!TOPO::check_boundary_boundary_cell_zero(entry.first))
                ++result.failed;
        }
        return result;
    }

    bool print_test_result(const std::string &name, bool passed, int myid)
    {
        if (myid == 0)
            std::cout << "  [" << (passed ? "PASS" : "FAIL") << "] " << name << "\n";
        return passed;
    }

    bool print_counted_test_result(const std::string &name, const TestSummary &summary, int myid)
    {
        const bool passed = summary.failed == 0;
        if (myid == 0)
        {
            std::cout << "  [" << (passed ? "PASS" : "FAIL") << "] " << name
                      << " checked=" << summary.checked
                      << " failed=" << summary.failed << "\n";
        }
        return passed;
    }
}

int main(int argc, char **argv)
{
    PARALLEL::mpi_initial(argc, argv);

    int exit_code = 0;
    try
    {
        int myid = 0;
        PARALLEL::mpi_rank(&myid);

        Param param;
        param.ReadParam(myid);

        Grid grid;
        grid.Grid_Preprocess(&param);

        const int dimension = param.GetInt("dimension");
        TOPO::Topology topology = TOPO::build_topology(grid, myid, dimension);
        TOPO::GlobalIncidence incidence(topology);

        if (myid == 0)
            std::cout << "Topology tests\n";

        bool passed = true;
        passed &= print_counted_test_result(
            "local boundary(boundary(face)) == 0",
            check_local_boundary_boundary_face(topology),
            myid);
        passed &= print_counted_test_result(
            "local boundary(boundary(cell)) == 0",
            check_local_boundary_boundary_cell(topology),
            myid);
        passed &= print_test_result("global D1 * D0 == 0",
                                    incidence.check_d1_d0_zero(),
                                    myid);
        passed &= print_test_result("global D2 * D1 == 0",
                                    incidence.check_d2_d1_zero(),
                                    myid);

        exit_code = passed ? 0 : 1;
    }
    catch (const std::exception &ex)
    {
        int myid = 0;
        PARALLEL::mpi_rank(&myid);
        if (myid == 0)
            std::cerr << "topo_test failed with exception: " << ex.what() << "\n";
        exit_code = 1;
    }

    PARALLEL::mpi_finalize();
    return exit_code;
}
