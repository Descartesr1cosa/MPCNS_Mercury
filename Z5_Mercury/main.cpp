//==============================================================================
//-------------->>>Multi-Physics Coupling Numerical Simulation<<<---------------
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>M P C N S<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//==============================================================================

//==============================================================================
#include "00_Mercury_Const.h"
#include "1_grid/1_MPCNS_Grid.h"
#include "0_basic/MPI_WRAPPER.h"
#include "2_topology/TopologyBuilder.h"
#include "2_topology/Topology.h"
#include "3_field/Field.h"
#include "4_halo/Halo.h"
#include "4_halo/HaloEdgeOwner.h"

#include "MercurySolver.h"

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
#if HALL_IMPLICIT == 1
    PetscInitialize(&arg, &argv, NULL, NULL); // 再 PETSc 初始化
#endif
    //=============================================================================================
    {
        //=============================================================================================
        //--------------------------------------------------------------------------
        // Read control parameters
        Param *par = new Param;
        par->ReadParam(myid);
        //--------------------------------------------------------------------------
        // Read Grid and Preprocess the Grid related info
        Grid *grd = new Grid;
        grd->Grid_Preprocess(par);
        //--------------------------------------------------------------------------
        // Build topology
        TOPO::Topology topology = TOPO::build_topology(*grd, myid, par->GetInt("dimension"));
        TOPO::build_topology_equivalence(topology, *grd, myid, par->GetInt("dimension"));
        //--------------------------------------------------------------------------
        int ngg = par->GetInt("ngg");
        // Build Field
        Field *fld = new Field(grd, par, ngg);
        MercurySolver::RegisterFields(fld, ngg);
        MercurySolver::RegisterCouplingChannels(fld, topology, par->GetInt("dimension"), ngg);
        //--------------------------------------------------------------------------
        // Build Halo Communicator
        Halo *hal = new Halo(fld, &topology);
        MercurySolver::RegisterHaloFields(fld, hal);
        //--------------------------------------------------------------------------
        // Build owner sync pattern once; both explicit and implicit edge fields use it.
        HALO_OWNER::EdgeOwnerSyncPattern edge_owner_pattern;
        HALO_OWNER::build_edge_owner_sync_pattern(topology, edge_owner_pattern);
        //=============================================================================================

        //=============================================================================================
        MercurySolver solver(grd, &topology, fld, hal, par,
                             &topology,
                             &edge_owner_pattern);
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

#if HALL_IMPLICIT == 1
    PetscFinalize(); // 先 PETSc finalize
#endif
    PARALLEL::mpi_finalize();
    return 0;
}
