#include "Z1_Control.h"

#include "0_basic/1_MPCNS_Parameter.h"
#include "0_basic/MPI_WRAPPER.h"
#include "1_grid/1_MPCNS_Grid.h"
#include "2_topology/TopologyBuilder.h"

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
        TOPO::Topology topology = TOPO::build_topology(grid, myid, dimension);
    }

    PARALLEL::mpi_finalize();
    return 0;
}
