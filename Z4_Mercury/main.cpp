//==============================================================================
//-------------->>>Multi-Physics Coupling Numerical Simulation<<<---------------
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>M P C N S<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//==============================================================================

//==============================================================================
#include "00_Mercury_Const.h"
#include "1_grid/1_MPCNS_Grid.h"
#include "0_basic/MPI_WRAPPER.h"
#include "2_topology/TopologyBuilder.h"
#include "7_metric/SingularEdgeRegistry.h"
#include "2_topology/Topology.h"
#include "3_field/Field.h"
#include "4_halo/Halo.h"
#include "4_halo/HaloEdgeOwner.h"

#include "MercurySolver.h"
#include "MercuryCase.h"
#include "MercuryGridValidation.h"

#if HALL_IMPLICIT == 1
// #include "4_solver/ImplicitHall_Solver.h"
#endif
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
        PetscInitialize(&arg, &argv, NULL, NULL); // 再 PETSc 初始化
    //=============================================================================================
    {
        //=============================================================================================
        //--------------------------------------------------------------------------
        // Read control parameters
        MERCURY::PrepareCaseWorkdirIfNeeded(myid);
        Param *par = new Param;
        par->ReadParam(myid);
        //--------------------------------------------------------------------------
        // Read Grid and Preprocess the Grid related info
        Grid *grd = new Grid;
        grd->Grid_Preprocess(par);
        MERCURY::ValidateCubicSphereGridOrAbort(*grd, par->GetInt("ngg"), myid);
        //--------------------------------------------------------------------------
        // Build topology
        TOPO::Topology topology = TOPO::build_topology(*grd, myid, par->GetInt("dimension"));
        //--------------------------------------------------------------------------
        int ngg = par->GetInt("ngg");
        // Build Field
        Field *fld = new Field(grd, par, ngg);
        MercurySolver::RegisterFields(fld, ngg);
        MercurySolver::RegisterCouplingChannels(fld, topology, par->GetInt("dimension"), ngg);
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
        MercurySolver::RegisterHaloFields(fld, hal);
        //--------------------------------------------------------------------------
        HALO_OWNER::EdgeOwnerSyncPattern edge_owner_pat;
        HALO_OWNER::build_edge_owner_sync_pattern(topology, edge_owner_pat);
        //=============================================================================================

        //=============================================================================================
        MercurySolver solver(grd, &topology, fld, hal, par,
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

    PetscFinalize(); // 先 PETSc finalize
    PARALLEL::mpi_finalize();
    return 0;
}
