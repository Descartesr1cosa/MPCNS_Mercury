#pragma once

#include "2_topology/2_MPCNS_Topology_Equiv.h"

#include <iosfwd>
#include <string>

class Field;
class Halo;

namespace Z0_NULL
{
    void print_banner();

    void print_diagnostics(const Field &fields,
                           const TOPO::TopologyEquiv &topology_equiv,
                           int dimension,
                           int nghost);

    void dump_halo_registry_if_requested(const Halo &halo, bool dump_registry);

    namespace DIAG
    {
        void dump_field_catalog(const Field &fields,
                                int my_rank,
                                std::ostream &os);

        void dump_field_block_summary(const Field &fields,
                                      int my_rank,
                                      std::ostream &os);

        void dump_topology_equiv_summary(const TOPO::TopologyEquiv &equiv,
                                         int my_rank,
                                         std::ostream &os);

        bool check_field_finite(const Field &fields,
                                const std::string &field_name,
                                int my_rank,
                                std::ostream &os);

        double check_edge_owner_alias_error(const Field &fields,
                                            const TOPO::TopologyEquiv &equiv,
                                            const std::string &field_name,
                                            int my_rank,
                                            std::ostream &os);

        double check_face_owner_alias_error(const Field &fields,
                                            const TOPO::TopologyEquiv &equiv,
                                            const std::string &field_name,
                                            int my_rank,
                                            std::ostream &os);

        bool run_basic_halo_validation(Field &fields,
                                       Halo &halo,
                                       const TOPO::TopologyEquiv &equiv,
                                       int my_rank,
                                       std::ostream &os);
    }
}
