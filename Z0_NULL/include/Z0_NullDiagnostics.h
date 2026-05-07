#pragma once

#include "2_topology/2_MPCNS_Topology_Equiv.h"

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
}
