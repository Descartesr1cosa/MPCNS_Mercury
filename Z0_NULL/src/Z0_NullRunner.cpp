#include "Z0_NullRunner.h"

#include "Z0_NullConfig.h"
#include "Z0_NullDiagnostics.h"
#include "Z0_NullFieldCatalog.h"
#include "Z0_NullIO.h"
#include "Z0_NullInitializer.h"
#include "0_basic/1_MPCNS_Parameter.h"
#include "0_basic/MPI_WRAPPER.h"
#include "1_grid/1_MPCNS_Grid.h"
#include "2_topology/2_MPCNS_Topology.h"
#include "2_topology/2_MPCNS_Topology_Equiv.h"
#include "3_field/2_MPCNS_Field.h"
#include "4_halo/1_MPCNS_Halo.h"

#include <iostream>
#include <memory>

namespace Z0_NULL
{
    int run(int argc, char **argv)
    {
        int myid = 0;
        PARALLEL::mpi_rank(&myid);

        const NullConfig cfg = parse_config(argc, argv);
        use_z4_case_workdir_if_needed(myid);

        if (myid == 0)
            print_banner();

        auto par = std::make_unique<Param>();
        par->ReadParam(myid);

        auto grd = std::make_unique<Grid>();
        grd->Grid_Preprocess(par.get());

        const int dimension = par->GetInt("dimension");
        const int nghost = par->GetInt("ngg");

        TOPO::Topology topology = TOPO::build_topology(*grd, myid, dimension);

        TOPO::TopologyEquiv topology_equiv;
        TOPO::build_topology_equiv(topology, *grd, myid, dimension, topology_equiv);

        auto fields = std::make_unique<Field>(grd.get(), par.get(), nghost);
        register_null_fields(*fields, nghost);
        initialize_null_fields(*fields);

        auto halo = std::make_unique<Halo>(fields.get(), &topology);
        halo->set_topology_equiv(&topology_equiv);
        halo->register_halo_fields(fields->halo_requests());
        halo->build_registered_patterns();

        if (cfg.sync_test)
        {
            initialize_halo_smoke_fields(*fields);
            halo->sync_registered();
            PARALLEL::mpi_barrier();

            if (myid == 0)
            {
                std::cout << "Deterministic halo sync smoke test completed.\n"
                          << std::flush;
            }
        }

        halo->sync_registered();
        PARALLEL::mpi_barrier();

        if (cfg.write_tecplot)
            write_tecplot_output(*par, *grd, *fields, cfg.output_step, cfg.output_time);

        if (myid == 0)
        {
            print_diagnostics(*fields, topology_equiv, dimension, nghost);
            dump_halo_registry_if_requested(*halo, cfg.dump_registry);

            std::cout << "Z0_NULL framework example finished.\n"
                      << std::flush;
        }

        return 0;
    }
}
