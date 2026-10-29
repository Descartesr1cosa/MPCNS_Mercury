//==============================================================================
//-------------->>>Multi-Physics Coupling Numerical Simulation<<<---------------
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>M P C N S<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//==============================================================================

//==============================================================================
#include "00_Mercury_Const.h"
#include "0_basic/MPI_WRAPPER.h"
#include "7_metric/SingularEdgeRegistry.h"
#include "mpcns/compat/legacy_singular_edge_collector.hpp"
#include "4_halo/HaloEdgeOwner.h"
#include "9_runtime/SimulationContext.h"

#include "MercurySolver.h"
#include "MercuryCase.h"
#include "0_MercuryFieldCatalog.h"
#include "MercuryGridValidation.h"

#include <utility>

//==============================================================================

//==============================================================================

//==============================================================================

int main(int arg, char **argv)
{
    //=============================================================================================
    // MPI/PETSc remain process-level services.  SimulationContext borrows
    // them; it must be destroyed before either library is finalized.
    int myid;
    PARALLEL::mpi_initial(arg, argv);
    PARALLEL::mpi_rank(&myid);
    PetscInitialize(&arg, &argv, NULL, NULL);
    //=============================================================================================
    {
        //=============================================================================================
        //--------------------------------------------------------------------------
        // Read control parameters
        MERCURY::PrepareCaseWorkdirIfNeeded(myid);
        Param par;
        par.ReadParam(myid);

        // The runtime is the sole owner of the generic numerical substrate:
        // Grid -> Topology -> geometry/Field -> Coupling -> Halo.  It freezes
        // the common registry collectively before Mercury creates any
        // case-specific object, preventing rank-dependent field layouts.
        const int ngg = par.GetInt("ngg");
        MPCNS::SimulationConfig runtime_config;
        runtime_config.geometry_ghost = ngg;
        runtime_config.fields = MERCURY_FIELD::RuntimeFieldDescriptors(ngg);
        runtime_config.coupling_pairs = MERCURY_FIELD::CouplingPairs();
        runtime_config.enable_metric_fields = true;
        // MercurySolver owns its legacy-compatible IOModule and the special
        // MercuryBoundary scheduler.  Disable only the runtime's optional IO;
        // Grid, topology, field, metric, and Halo ownership stay centralized.
        runtime_config.io.enabled = false;
        runtime_config.physical_corner_mode = MPCNS::PhysicalCornerMode::SolverManaged;

        MPCNS::SimulationContext runtime(par, myid, std::move(runtime_config));
        runtime.initialize();
        MERCURY::ValidateCubicSphereGridOrAbort(runtime.grid(), ngg, myid);

        // The generic runtime deliberately has no Mercury physics.  Build the
        // variable-valence physical-edge registry only after its immutable
        // topology and metric fields exist.
        mpcns::compat::LegacySingularEdgeDataSource singular_edge_source{
            runtime.field(), runtime.grid()};
        METRIC::SingularEdgeRegistry singular_edges;
        singular_edges.build(runtime.topology(), singular_edge_source, myid);
        singular_edges.validate_or_abort();
        if (myid == 0)
            std::cout << "[SingularEdgeRegistry] physical singular edges="
                      << singular_edges.size() << "\n";
        //--------------------------------------------------------------------------
        HALO_OWNER::EdgeOwnerSyncPattern edge_owner_pat;
        HALO_OWNER::build_edge_owner_sync_pattern(runtime.topology(), edge_owner_pat);
        //=============================================================================================

        //=============================================================================================
        MercurySolver solver(runtime, &edge_owner_pat, &singular_edges);
        solver.Advance();
        if (myid == 0)
            std::cout << "Program is finished normally ! !  ^_^\n"
                      << std::flush;
        //=============================================================================================

        // `solver`, the singular-edge registry, and `runtime` are destroyed in
        // reverse order before PETSc/MPI finalization below.
    }
    //=============================================================================================
    //--------------------------------------------------------------------------
    // MPI finalization

    PetscFinalize(); // 先 PETSc finalize
    PARALLEL::mpi_finalize();
    return 0;
}
