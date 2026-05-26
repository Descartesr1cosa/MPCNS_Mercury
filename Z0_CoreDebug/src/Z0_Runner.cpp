#include "Z0_Runner.h"

#include "Z0_Config.h"
#include "Z0_CouplingTests.h"
#include "Z0_Diagnostics.h"
#include "Z0_FieldCatalog.h"
#include "Z0_Initializer.h"
#include "Z0_OutputTests.h"
#include "Z0_SyncTests.h"
#include "0_basic/1_MPCNS_Parameter.h"
#include "0_basic/MPI_WRAPPER.h"
#include "1_grid/1_MPCNS_Grid.h"
#include "2_topology/TopologyBuilder.h"
#include "2_topology/Topology.h"
#include "3_field/Field.h"
#include "4_halo/Halo.h"

#include <iostream>
#include <memory>
#include <set>
#include <vector>

namespace
{
    void collect_coupling_pairs(const TOPO::Topology &topology,
                                std::vector<Field::PairKey> &pairs)
    {
        std::set<Field::PairKey> unique;
        auto add = [&](const auto &p)
        {
            if (!p.is_coupling)
                return;
            if (p.nb_block_name.empty() || p.this_block_name.empty())
                return;
            unique.insert(Field::PairKey{p.nb_block_name, p.this_block_name});
        };
        for (const auto &p : topology.inner_patches)
            add(p);
        for (const auto &p : topology.parallel_patches)
            add(p);
        for (const auto &p : topology.inner_edge_patches)
            add(p);
        for (const auto &p : topology.parallel_edge_patches)
            add(p);
        for (const auto &p : topology.inner_vertex_patches)
            add(p);
        for (const auto &p : topology.parallel_vertex_patches)
            add(p);
        pairs.assign(unique.begin(), unique.end());
    }

    void merge_result(Z0::TestResult &total, const Z0::TestResult &one)
    {
        total.pass = total.pass && one.pass;
        total.max_error = std::max(total.max_error, one.max_error);
    }
}

namespace Z0
{
    int run(int argc, char **argv)
    {
        int myid = 0;
        PARALLEL::mpi_rank(&myid);

        const Config cfg = parse_config(argc, argv);
        use_z4_case_workdir_if_needed(myid);

        if (myid == 0)
        {
            print_banner();
            std::cout << "[Z0_CoreDebug] mode = " << mode_name(cfg.mode) << "\n";
        }

        auto par = std::make_unique<Param>();
        par->ReadParam(myid);

        auto grd = std::make_unique<Grid>();
        grd->Grid_Preprocess(par.get());

        const int dimension = par->GetInt("dimension");
        const int nghost = par->GetInt("ngg");

        TOPO::Topology topology = TOPO::build_topology(*grd, myid, dimension);
        TOPO::build_topology_equivalence(topology, *grd, myid, dimension);

        auto fields = std::make_unique<Field>(grd.get(), par.get(), nghost);
        register_core_debug_fields(*fields, nghost);

        std::vector<Field::PairKey> coupling_pairs;
        collect_coupling_pairs(topology, coupling_pairs);
        fields->register_declared_coupling_channels(coupling_pairs);
        fields->build_coupling_buffers(topology, dimension);

        auto halo = std::make_unique<Halo>(fields.get(), &topology);
        halo->set_topology_equiv(&topology);
        halo->register_halo_fields(fields->halo_requests());
        halo->build_registered_patterns();

        InitContext init_ctx;
        init_ctx.my_rank = myid;
        init_ctx.dimension = dimension;
        initialize_all_fields(*fields, init_ctx);

        TestResult total;
        if (cfg.mode == Mode::Summary || cfg.mode == Mode::All)
        {
            if (myid == 0)
            {
                print_diagnostics(*fields, topology, topology, dimension, nghost, std::cout);
                halo->dump_sync_registry(std::cout);
            }
        }

        if (cfg.mode == Mode::Sync || cfg.mode == Mode::All)
        {
            merge_result(total, test_field_extents(*fields, *grd, dimension, myid, std::cout));
            merge_result(total, test_component_halo(*fields, *halo, topology, myid, std::cout));
            merge_result(total, test_edge_1form_triplet_halo(*fields, *halo, topology, myid, std::cout));
            merge_result(total, test_face_2form_triplet_halo(*fields, *halo, topology, myid, std::cout));
            merge_result(total, test_owner_alias_sync(*fields, *halo, topology, myid, std::cout));
            merge_result(total, test_sync_group_order(*fields, topology, myid, std::cout));
        }

        if (cfg.mode == Mode::Coupling || cfg.mode == Mode::All)
        {
            halo->sync_registered();
            merge_result(total, test_coupling_cell_scalar(*fields, *halo, topology, dimension, myid, std::cout));
            merge_result(total, test_coupling_edge_1form(*fields, *halo, topology, dimension, myid, std::cout));
            merge_result(total, test_coupling_face_2form(*fields, *halo, topology, dimension, myid, std::cout));
        }

        if (cfg.mode == Mode::Output || cfg.mode == Mode::All)
            merge_result(total, test_location_output_smoke(*par, *grd, *fields, myid, std::cout));

        double pass_local = total.pass ? 1.0 : 0.0;
        double pass_global = pass_local;
        double err_local = total.max_error;
        double err_global = err_local;
        PARALLEL::mpi_min(&pass_local, &pass_global, 1);
        PARALLEL::mpi_max(&err_local, &err_global, 1);

        const bool pass = pass_global > 0.5;
        if (myid == 0)
        {
            std::cout << "Z0_CoreDebug result: " << (pass ? "PASS" : "FAIL")
                      << " max_error=" << err_global << "\n";
        }
        return pass ? 0 : 1;
    }
}
