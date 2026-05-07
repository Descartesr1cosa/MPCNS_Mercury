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
        {
            print_banner();
            std::cout << "[Z0_NULL] mode = " << mode_name(cfg.mode) << "\n"
                      << std::flush;
        }

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

        auto halo = std::make_unique<Halo>(fields.get(), &topology);
        halo->set_topology_equiv(&topology_equiv);
        halo->register_halo_fields(fields->halo_requests());
        halo->build_registered_patterns();

        InitContext init_ctx;
        init_ctx.my_rank = myid;
        init_ctx.dimension = dimension;
        initialize_all_fields(*fields, init_ctx);

        bool pass = true;
        const bool do_sync = (cfg.mode == NullMode::Sync || cfg.mode == NullMode::All);
        const bool do_io = cfg.write_tecplot;

        if (do_sync)
        {
            pass = DIAG::run_basic_halo_validation(*fields, *halo, topology_equiv, myid, std::cout) && pass;
            PARALLEL::mpi_barrier();
        }

        if (myid == 0)
        {
            print_diagnostics(*fields, topology_equiv, dimension, nghost);
            DIAG::dump_topology_equiv_summary(topology_equiv, myid, std::cout);
            DIAG::dump_field_catalog(*fields, myid, std::cout);
            DIAG::dump_field_block_summary(*fields, myid, std::cout);
            halo->dump_sync_registry(std::cout);
        }

        if (do_io)
            write_null_tecplot(*par, *grd, *fields, cfg.output_step, cfg.output_time);

        double pass_local = pass ? 1.0 : 0.0;
        double pass_global = pass_local;
        PARALLEL::mpi_min(&pass_local, &pass_global, 1);
        pass = (pass_global > 0.5);

        if (myid == 0)
        {
            std::cout << "Z0_NULL result: " << (pass ? "PASS" : "FAIL") << "\n"
                      << std::flush;
        }

        return pass ? 0 : 1;
    }
}
