#include "Z0_FieldCatalog.h"

#include "3_field/Field.h"
#include "4_halo/Halo.h"

#include <cstdlib>
#include <set>
#include <string>
#include <vector>

namespace
{
    FieldSyncContract sync_contract(const char *group,
                                    bool do_halo,
                                    bool do_physical,
                                    bool do_coupling,
                                    HaloLevel level = HaloLevel::Corner3D,
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

    std::string env_string(const char *name, const char *fallback)
    {
        const char *value = std::getenv(name);
        return value ? std::string(value) : std::string(fallback);
    }

    HaloLevel env_halo_level()
    {
        const std::string level = env_string("Z0_HALO_LEVEL", "face");
        if (level == "edge" || level == "Edge" || level == "2d" || level == "2D")
            return HaloLevel::Corner2D;
        if (level == "vertex" || level == "Vertex" || level == "3d" || level == "3D")
            return HaloLevel::Corner3D;
        return HaloLevel::Corner1D;
    }

    bool halo_enabled_for(const std::string &mode, const char *category)
    {
        if (mode == "all")
            return true;
        if (mode == "none")
            return false;
        if (mode == "node")
            return std::string(category) == "node";
        if (mode == "forms")
            return std::string(category) == "forms";
        if (mode == "edgeforms")
            return std::string(category) == "edgeforms";
        if (mode == "faceforms")
            return std::string(category) == "faceforms";
        return std::string(category) == "cell";
    }

    bool env_owner_sync_enabled()
    {
        const std::string value = env_string("Z0_OWNER_SYNC", "0");
        return value == "1" || value == "true" || value == "TRUE" || value == "on" || value == "ON";
    }

}

namespace Z0
{
    void RegisterFields(Field &field, Param &param, int nghost)
    {
        (void)param;
        const HaloLevel default_level = env_halo_level();
        const std::string halo_mode = env_string("Z0_HALO_MODE", "cell");
        const bool halo_node = halo_enabled_for(halo_mode, "node");
        const bool halo_edge_forms = halo_enabled_for(halo_mode, "forms") || halo_enabled_for(halo_mode, "edgeforms");
        const bool halo_face_forms = halo_enabled_for(halo_mode, "forms") || halo_enabled_for(halo_mode, "faceforms");
        const bool halo_cell = halo_enabled_for(halo_mode, "cell");
        const bool owner_sync = env_owner_sync_enabled();
        const bool form_orientation = env_string("Z0_FORM_ORIENTATION", "1") != "0";
        field.register_field(descriptor("phi", StaggerLocation::Node, FieldValueKind::Scalar,
                                        1, nghost, sync_contract("phi", halo_node, false, false,
                                                                 default_level, false,
                                                                 owner_sync ? OwnerSyncPolicy::NodeOwner : OwnerSyncPolicy::None)));
        field.register_field(descriptor("E_xi", StaggerLocation::EdgeXi, FieldValueKind::EdgeCovariant1Form,
                                        1, nghost, sync_contract("Eedge", halo_edge_forms, false, false,
                                                                 default_level, form_orientation,
                                                                 owner_sync ? OwnerSyncPolicy::EdgeOwner : OwnerSyncPolicy::None)));
        field.register_field(descriptor("E_eta", StaggerLocation::EdgeEt, FieldValueKind::EdgeCovariant1Form,
                                        1, nghost, sync_contract("Eedge", halo_edge_forms, false, false,
                                                                 default_level, form_orientation,
                                                                 owner_sync ? OwnerSyncPolicy::EdgeOwner : OwnerSyncPolicy::None)));
        field.register_field(descriptor("E_zeta", StaggerLocation::EdgeZe, FieldValueKind::EdgeCovariant1Form,
                                        1, nghost, sync_contract("Eedge", halo_edge_forms, false, false,
                                                                 default_level, form_orientation,
                                                                 owner_sync ? OwnerSyncPolicy::EdgeOwner : OwnerSyncPolicy::None)));
        field.register_field(descriptor("B_xi", StaggerLocation::FaceXi, FieldValueKind::FaceContravariant2Form,
                                        1, nghost, sync_contract("Bface", halo_face_forms, false, false,
                                                                 default_level, form_orientation,
                                                                 owner_sync ? OwnerSyncPolicy::FaceOwner : OwnerSyncPolicy::None)));
        field.register_field(descriptor("B_eta", StaggerLocation::FaceEt, FieldValueKind::FaceContravariant2Form,
                                        1, nghost, sync_contract("Bface", halo_face_forms, false, false,
                                                                 default_level, form_orientation,
                                                                 owner_sync ? OwnerSyncPolicy::FaceOwner : OwnerSyncPolicy::None)));
        field.register_field(descriptor("B_zeta", StaggerLocation::FaceZe, FieldValueKind::FaceContravariant2Form,
                                        1, nghost, sync_contract("Bface", halo_face_forms, false, false,
                                                                 default_level, form_orientation,
                                                                 owner_sync ? OwnerSyncPolicy::FaceOwner : OwnerSyncPolicy::None)));
        field.register_field(descriptor("divB", StaggerLocation::Cell, FieldValueKind::Scalar,
                                        1, nghost, sync_contract("divB", halo_cell, false, false,
                                                                 default_level)));
        field.register_field(descriptor("U", StaggerLocation::Cell, FieldValueKind::ConservedVector,
                                        5, nghost, sync_contract("U", halo_cell, false, false,
                                                                 default_level)));
        field.register_field(descriptor("Bcell", StaggerLocation::Cell, FieldValueKind::CartesianVector,
                                        3, nghost, sync_contract("Bcell", halo_cell, false, false,
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
