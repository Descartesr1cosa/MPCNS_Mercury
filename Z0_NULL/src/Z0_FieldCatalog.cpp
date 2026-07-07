#include "Z0_FieldCatalog.h"

#include "3_field/Field.h"
#include "4_halo/Halo.h"

#include <set>
#include <vector>

namespace
{
    FieldSyncContract sync_contract(const char *group,
                                    bool do_halo,
                                    bool do_physical,
                                    bool do_coupling,
                                    HaloLevel level = HaloLevel::Vertex,
                                    bool orientation_aware = false,
                                    OwnerSyncPolicy owner_sync = OwnerSyncPolicy::None)
    {
        FieldSyncContract sync;
        sync.group = group ? group : "";
        sync.do_halo = do_halo;
        sync.do_physical = do_physical;
        sync.do_coupling = do_coupling;
        sync.halo_level = level;
        sync.orientation_aware = orientation_aware;
        sync.owner_sync = owner_sync;
        return sync;
    }

    FieldDescriptor descriptor(const char *name,
                               StaggerLocation location,
                               FieldValueKind kind,
                               int ncomp,
                               int nghost,
                               FieldSyncContract sync)
    {
        FieldDescriptor d;
        d.name = name;
        d.location = location;
        d.value_kind = kind;
        d.ncomp = ncomp;
        d.nghost = nghost;
        d.physics = "";
        d.sync = sync;
        return d;
    }

    void collect_coupling_pairs(const TOPO::Topology &topology,
                                std::vector<Field::PairKey> &pairs)
    {
        std::set<Field::PairKey> unique;
        auto add = [&](const auto &p)
        {
            if (!p.is_coupling || p.nb_block_name.empty() || p.this_block_name.empty())
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
}

namespace Z0
{
    void RegisterFields(Field &field, Param &param, int nghost)
    {
        (void)param;
        const HaloLevel default_level = HaloLevel::FaceOnly;
        field.register_field(descriptor("phi", StaggerLocation::Node, FieldValueKind::Scalar,
                                        1, nghost, sync_contract("phi", false, false, false,
                                                                 default_level)));
        field.register_field(descriptor("E_xi", StaggerLocation::EdgeXi, FieldValueKind::EdgeCovariant1Form,
                                        1, nghost, sync_contract("Eedge", false, false, false,
                                                                 default_level)));
        field.register_field(descriptor("E_eta", StaggerLocation::EdgeEt, FieldValueKind::EdgeCovariant1Form,
                                        1, nghost, sync_contract("Eedge", false, false, false,
                                                                 default_level)));
        field.register_field(descriptor("E_zeta", StaggerLocation::EdgeZe, FieldValueKind::EdgeCovariant1Form,
                                        1, nghost, sync_contract("Eedge", false, false, false,
                                                                 default_level)));
        field.register_field(descriptor("B_xi", StaggerLocation::FaceXi, FieldValueKind::FaceContravariant2Form,
                                        1, nghost, sync_contract("Bface", false, false, false,
                                                                 default_level)));
        field.register_field(descriptor("B_eta", StaggerLocation::FaceEt, FieldValueKind::FaceContravariant2Form,
                                        1, nghost, sync_contract("Bface", false, false, false,
                                                                 default_level)));
        field.register_field(descriptor("B_zeta", StaggerLocation::FaceZe, FieldValueKind::FaceContravariant2Form,
                                        1, nghost, sync_contract("Bface", false, false, false,
                                                                 default_level)));
        field.register_field(descriptor("divB", StaggerLocation::Cell, FieldValueKind::Scalar,
                                        1, nghost, sync_contract("divB", true, false, false,
                                                                 default_level)));
        field.register_field(descriptor("U", StaggerLocation::Cell, FieldValueKind::ConservedVector,
                                        5, nghost, sync_contract("U", true, false, false,
                                                                 default_level)));
        field.register_field(descriptor("Bcell", StaggerLocation::Cell, FieldValueKind::CartesianVector,
                                        3, nghost, sync_contract("Bcell", true, false, false,
                                                                 default_level)));
    }

    void RegisterCouplingChannels(Field &field,
                                  const TOPO::Topology &topology,
                                  Param &param,
                                  int nghost,
                                  int dimension)
    {
        (void)param;
        (void)nghost;
        std::vector<Field::PairKey> pairs;
        collect_coupling_pairs(topology, pairs);
        if (pairs.empty())
            return;

        field.register_declared_coupling_channels(pairs);
        field.build_coupling_buffers(topology, dimension);
    }

    void RegisterHaloFields(Field &field, Halo &halo)
    {
        halo.register_halo_fields(field.halo_requests());
        halo.build_registered_patterns();
    }
}
