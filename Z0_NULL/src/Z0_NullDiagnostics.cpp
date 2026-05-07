#include "Z0_NullDiagnostics.h"

#include "0_basic/LayoutTraits.h"
#include "3_field/2_MPCNS_Field.h"
#include "4_halo/1_MPCNS_Halo.h"

#include <iostream>

namespace Z0_NULL
{
    void print_banner()
    {
        std::cout << "========== Z0_NULL Synchronization Framework Example ==========\n"
                  << "Purpose:\n"
                  << "  Build Grid / Topology / TopologyEquiv / Field / Halo.\n"
                  << "  Register representative fields.\n"
                  << "  Initialize representative physical fields.\n"
                  << "  Exercise synchronization framework mechanisms.\n"
                  << "  Write Tecplot output.\n"
                  << "  No physical solver is executed.\n"
                  << std::flush;
    }

    void print_diagnostics(const Field &fields,
                           const TOPO::TopologyEquiv &topology_equiv,
                           int dimension,
                           int nghost)
    {
        std::cout << "---------- Z0_NULL Diagnostics ----------\n"
                  << "dimension              = " << dimension << "\n"
                  << "ngg                    = " << nghost << "\n"
                  << "local blocks           = " << fields.num_blocks() << "\n"
                  << "registered fields      = " << fields.num_fields() << "\n"
                  << "halo requests          = " << fields.halo_requests().size() << "\n"
                  << "edge owner classes     = " << topology_equiv.edge_classes_general.size() << "\n"
                  << "face owner classes     = " << topology_equiv.face_classes.size() << "\n"
                  << "local edge owners      = " << topology_equiv.n_local_edge_owner << "\n"
                  << "global edge owners     = " << topology_equiv.n_global_edge_owner << "\n"
                  << "local face owners      = " << topology_equiv.n_local_face_owner << "\n"
                  << "global face owners     = " << topology_equiv.n_global_face_owner << "\n"
                  << std::flush;

        std::cout << "registered field names:\n";
        for (const auto &desc : fields.descriptors())
        {
            std::cout << "  - " << desc.name
                      << " loc=" << LAYOUT::location_name(desc.location)
                      << " kind=" << field_value_kind_name(desc.value_kind)
                      << " halo=" << (desc.sync.do_halo ? "true" : "false")
                      << "\n";
        }
        std::cout << std::flush;
    }

    void dump_halo_registry_if_requested(const Halo &halo, bool dump_registry)
    {
        if (dump_registry)
            halo.dump_sync_registry(std::cout);
    }
}
