#include "Z0_Control.h"
#include "Z0_Boundary.h"
#include "Z0_FieldCatalog.h"
#include "Z0_Solver.h"
#include "Z0_Tests.h"

#include "0_basic/1_MPCNS_Parameter.h"
#include "0_basic/MPI_WRAPPER.h"
#include "1_grid/1_MPCNS_Grid.h"
#include "2_topology/TopologyBuilder.h"
#include "3_field/Field.h"
#include "4_halo/Halo.h"

#include <iostream>

int main(int argc, char **argv)
{
    PARALLEL::mpi_initial(argc, argv);

    int exit_code = 0;
    {
        int myid = 0;
        PARALLEL::mpi_rank(&myid);
        Z0::PrepareCaseWorkdirIfNeeded(myid);

        Param param;
        param.ReadParam(myid);

        Grid grid;
        grid.Grid_Preprocess(&param);

        const int dimension = param.GetInt("dimension");
        TOPO::Topology topology = TOPO::build_topology(grid, myid, dimension);

        const int nghost = param.GetInt("ngg");
        Field field(&grid, &param, nghost);
        Z0::RegisterFields(field, param, nghost);
        Z0::RegisterCouplingChannels(field, topology, param, nghost, dimension);

        Halo halo(&field, &topology);
        halo.set_topology_equiv(&topology);
        Z0::RegisterHaloFields(field, halo);

        Z0_Boundary boundary(&grid, &field, &halo, &topology, &topology, &param);
        Z0_Solver solver(&grid, &field, &halo, &topology, &topology, &boundary, &param);
        solver.Advance();

        bool tests_passed = true;
        tests_passed &= Z0_TEST::RunTopologyTests(topology, myid);
        tests_passed &= Z0_TEST::RunPhysicalIoTests(field, param);
        tests_passed &= Z0_TEST::RunHaloCommunicationTests(field, halo, boundary, param);
        tests_passed &= Z0_TEST::RunDecChainTests(field, halo, boundary, param);

        if (myid == 0)
            std::cout << "Z0_NULL finished " << (tests_passed ? "normally" : "with test failures") << ".\n" << std::flush;

        if (!tests_passed)
            exit_code = 1;
    }

    PARALLEL::mpi_finalize();
    return exit_code;
}
