#include "0_basic/MPI_WRAPPER.h"
#include "1_grid/1_MPCNS_Grid.h"
#include "2_topology/TopologyBuilder.h"
#include "3_field/Field.h"
#include "4_halo/Halo.h"
#include "7_metric/SingularEdgeRegistry.h"
#include "MarsCase.h"
#include "MarsSolver.h"

#include <iostream>

int main(int argc, char **argv)
{
    PARALLEL::mpi_initial(argc, argv);
    int rank = 0;
    PARALLEL::mpi_rank(&rank);
    {
        MARS::PrepareCaseWorkdirIfNeeded(rank);
        Param param;
        param.ReadParam(rank);
        Grid grid;
        grid.Grid_Preprocess(&param);
        TOPO::Topology topology = TOPO::build_topology(grid, rank, param.GetInt("dimension"));
        const int ngh = param.GetInt("ngg");
        Field field(&grid, &param, ngh);
        MarsSolver::RegisterFields(&field, ngh);
        METRIC::SingularEdgeRegistry singular_edges;
        singular_edges.build(topology, field, grid, rank);
        singular_edges.validate_or_abort();
        if (rank == 0)
            std::cout << "[SingularEdgeRegistry] physical singular edges="
                      << singular_edges.size() << '\n';
        Halo halo(&field, &topology);
        halo.set_topology_equiv(&topology);
        MarsSolver::RegisterHaloFields(&field, &halo);
        MarsSolver solver(&grid, &topology, &field, &halo, &param, &singular_edges);
        solver.Advance();
        if (rank == 0)
            std::cout << "MarsZ5 finished normally\n";
    }
    PARALLEL::mpi_finalize();
    return 0;
}
