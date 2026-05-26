#include "Z0_Diagnostics.h"

#include "0_basic/LayoutTraits.h"
#include "3_field/Field.h"
#include "4_halo/Halo.h"

#include <cmath>
#include <iostream>

namespace
{
    const char *halo_level_name(HaloLevel level)
    {
        switch (level)
        {
        case HaloLevel::FaceOnly:
            return "FaceOnly";
        case HaloLevel::Edge:
            return "Edge";
        case HaloLevel::Vertex:
            return "Vertex";
        }
        return "Unknown";
    }
}

namespace Z0
{
    void print_banner()
    {
        std::cout << "========== Z0_CoreDebug Framework Test ==========\n"
                  << "Field / Topology / Topology / Halo / Coupling / IO validation only.\n"
                  << "No Mercury PDE solver or implicit Hall path is executed.\n"
                  << "=================================================\n";
    }

    void dump_field_catalog(const Field &fields, int my_rank, std::ostream &os)
    {
        os << "[rank " << my_rank << "] FieldCatalog (" << fields.num_fields() << " fields)\n";
        for (int fid = 0; fid < fields.num_fields(); ++fid)
        {
            const FieldDescriptor &d = fields.descriptor(fid);
            os << "  fid=" << fid
               << " name=" << d.name
               << " loc=" << LAYOUT::location_name(d.location)
               << " kind=" << field_value_kind_name(d.value_kind)
               << " ncomp=" << d.ncomp
               << " nghost=" << d.nghost
               << " group=" << d.sync.group
               << " halo=" << (d.sync.do_halo ? "true" : "false")
               << " coupling=" << (d.sync.do_coupling ? "true" : "false")
               << " level=" << halo_level_name(d.sync.halo_level)
               << " owner=" << owner_sync_policy_name(d.sync.owner_sync)
               << " orientation=" << (d.sync.orientation_aware ? "true" : "false")
               << "\n";
        }
    }

    void dump_topology_equiv_summary(const TOPO::Topology &equiv, int my_rank, std::ostream &os)
    {
        os << "[rank " << my_rank << "] Topology summary"
           << " node_classes=" << equiv.classes(TOPO::EntityDim::Node).size()
           << " edge_classes=" << equiv.classes(TOPO::EntityDim::Edge).size()
           << " face_classes=" << equiv.classes(TOPO::EntityDim::Face).size()
           << "\n";
    }

    void print_diagnostics(const Field &fields,
                           const TOPO::Topology &topology,
                           const TOPO::Topology &equiv,
                           int dimension,
                           int nghost,
                           std::ostream &os)
    {
        os << "[Z0] dimension=" << dimension
           << " nghost=" << nghost
           << " inner_faces=" << topology.inner_patches.size()
           << " parallel_faces=" << topology.parallel_patches.size()
           << " physical_faces=" << topology.physical_patches.size()
           << "\n";
        dump_field_catalog(fields, 0, os);
        dump_topology_equiv_summary(equiv, 0, os);
    }

    void report_test(const std::string &name, const TestResult &result, std::ostream &os)
    {
        os << "[Z0][" << name << "] "
           << (result.pass ? "PASS" : "FAIL")
           << " max_error=" << result.max_error << "\n";
    }

    void report_failure(const char *test,
                        const Field &fields,
                        const FieldDescriptor &desc,
                        int rank,
                        int block,
                        int i,
                        int j,
                        int k,
                        int comp,
                        double expected,
                        double actual,
                        std::ostream &os,
                        const TOPO::IndexTransform *tr)
    {
        os << "[Z0][" << test << "] failure"
           << " field=" << desc.name
           << " rank=" << rank
           << " block=" << block
           << " i=" << i
           << " j=" << j
           << " k=" << k
           << " comp=" << comp
           << " expected=" << expected
           << " actual=" << actual
           << " abs_error=" << std::abs(actual - expected)
           << " location=" << LAYOUT::location_name(desc.location)
           << " value_kind=" << field_value_kind_name(desc.value_kind);

        if (tr)
        {
            os << " transform_perm=(" << tr->perm[0] << "," << tr->perm[1] << "," << tr->perm[2] << ")"
               << " transform_sign=(" << tr->sign[0] << "," << tr->sign[1] << "," << tr->sign[2] << ")"
               << " transform_offset=(" << tr->offset.i << "," << tr->offset.j << "," << tr->offset.k << ")";
        }
        os << "\n";
        (void)fields;
    }
}
