//==============================================================================
//-------------->>>Multi-Physics Coupling Numerical Simulation<<<---------------
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>M P C N S<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//==============================================================================

//==============================================================================
#include "00_Lunar_Const.h"
#include "1_grid/1_MPCNS_Grid.h"
#include "0_basic/MPI_WRAPPER.h"
#include "2_topology/TopologyBuilder.h"
#include "7_metric/SingularEdgeRegistry.h"
#include "2_topology/Topology.h"
#include "3_field/Field.h"
#include "4_halo/Halo.h"
#include "4_halo/HaloEdgeOwner.h"

#include "LunarSolver.h"
#include "LunarCase.h"
#include "LunarGridValidation.h"

//==============================================================================

//==============================================================================

//==============================================================================

int main(int arg, char **argv)
{
    //=============================================================================================
    // MPI initialization
        int myid;
    PARALLEL::mpi_initial(arg, argv);
    PARALLEL::mpi_rank(&myid);
    //=============================================================================================
    {
        //=============================================================================================
        //--------------------------------------------------------------------------
        // Read control parameters
        LUNAR::PrepareCaseWorkdirIfNeeded(myid);
        Param *par = new Param;
        par->ReadParam(myid);
        //--------------------------------------------------------------------------
        // Read Grid and Preprocess the Grid related info
        Grid *grd = new Grid;
        grd->Grid_Preprocess(par);
        // Z3_Lunar has a single physical region: every grid block is fluid.
        for (int ib = 0; ib < grd->nblock; ++ib)
            grd->grids(ib).block_name = "Fluid";
        LUNAR::ValidateCubicSphereGridOrAbort(*grd, par->GetInt("ngg"), myid);
        //--------------------------------------------------------------------------
        // Build topology
        TOPO::Topology topology = TOPO::build_topology(*grd, myid, par->GetInt("dimension"));
        //--------------------------------------------------------------------------
        int ngg = par->GetInt("ngg");
        // Build Field
        Field *fld = new Field(grd, par, ngg);
        LunarSolver::RegisterFields(fld, ngg);
        LunarSolver::RegisterCouplingChannels(fld, topology, par->GetInt("dimension"), ngg);
        METRIC::SingularEdgeRegistry singular_edges;
        singular_edges.build(topology, *fld, *grd, myid);
        singular_edges.validate_or_abort();
        if (myid == 0)
            std::cout << "[SingularEdgeRegistry] physical singular edges="
                      << singular_edges.size() << "\n";
        //--------------------------------------------------------------------------
        // Build Halo Communicator
        Halo *hal = new Halo(fld, &topology);
        hal->set_topology_equiv(&topology);
        LunarSolver::RegisterHaloFields(fld, hal);
        //--------------------------------------------------------------------------
        HALO_OWNER::EdgeOwnerSyncPattern edge_owner_pat;
        HALO_OWNER::build_edge_owner_sync_pattern(topology, edge_owner_pat);
        //=============================================================================================

        //=============================================================================================
        LunarSolver solver(grd, &topology, fld, hal, par,
                             &topology,
                             &edge_owner_pat,
                             &singular_edges);
        solver.Advance();
        if (myid == 0)
            std::cout << "Program is finished normally ! !  ^_^\n"
                      << std::flush;
        //=============================================================================================

        //--------------------------------------------------------------------------
        // Release the allocated memory, recommend to release with inverse order of building
        delete hal;
        delete fld;
        delete par;
        delete grd;
        //=============================================================================================
    }
    //=============================================================================================
    //--------------------------------------------------------------------------
    // MPI finalization

    PARALLEL::mpi_finalize();
    return 0;
}
