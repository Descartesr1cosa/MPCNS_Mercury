// Z1_NULL is a full empty physics template.
// It demonstrates the standard structure for adding a new solver module.
// It does not solve any PDE.
// For low-level framework tests, use Z0_CoreDebug.
// For real Mercury MHD/Hall-MHD, use Z4_Mercury.

#include "Z1_Boundary.h"
#include "Z1_Control.h"
#include "Z1_FieldCatalog.h"
#include "Z1_Solver.h"

#include "0_basic/1_MPCNS_Parameter.h"
#include "0_basic/MPI_WRAPPER.h"
#include "1_grid/1_MPCNS_Grid.h"
#include "2_topology/TopologyBuilder.h"
#include "2_topology/TopologyEquiv.h"
#include "3_field/Field.h"
#include "4_halo/Halo.h"

int main(int argc, char **argv)
{
    PARALLEL::mpi_initial(argc, argv);

    {
        int myid = 0;
        PARALLEL::mpi_rank(&myid);
        Z1::PrepareCaseWorkdirIfNeeded(myid);

        Param param;
        param.ReadParam(myid);

        Grid grid;
        grid.Grid_Preprocess(&param);

        const int dimension = param.GetInt("dimension");
        const int nghost = param.GetInt("ngg");

        TOPO::Topology topology = TOPO::build_topology(grid, myid, dimension);

        TOPO::TopologyEquiv topology_equiv;
        TOPO::build_topology_equiv(topology, grid, myid, dimension, topology_equiv);

        Field field(&grid, &param, nghost);
        Z1::RegisterFields(field, param, nghost);
        Z1::RegisterCouplingChannels(field, topology, param, nghost, dimension);

        Halo halo(&field, &topology);
        halo.set_topology_equiv(&topology_equiv);

        Z1::RegisterHaloFields(field, halo);
        halo.build_registered_patterns();

        Z1_Boundary boundary(&grid, &field, &halo, &topology, &topology_equiv, &param);

        Z1_Solver solver(&grid, &field, &halo, &topology, &topology_equiv, &boundary, &param);
        solver.Initialize();
        solver.Advance();
    }

    PARALLEL::mpi_finalize();
    return 0;
}
