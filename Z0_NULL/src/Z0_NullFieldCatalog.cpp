#include "Z0_NullFieldCatalog.h"

#include "3_field/2_MPCNS_Field.h"

#ifndef Z0_ENABLE_EDGE_OWNER
#define Z0_ENABLE_EDGE_OWNER 1
#endif

#ifndef Z0_ENABLE_FACE_OWNER
#define Z0_ENABLE_FACE_OWNER 1
#endif

namespace
{
    FieldSyncContract sync_contract(const char *group,
                                    bool do_halo,
                                    HaloLevel halo_level = HaloLevel::Vertex,
                                    bool orientation_aware = false,
                                    OwnerSyncPolicy owner_sync = OwnerSyncPolicy::None)
    {
        FieldSyncContract sync;
        sync.group = group ? group : "";
        sync.do_halo = do_halo;
        sync.halo_level = halo_level;
        sync.orientation_aware = orientation_aware;
        sync.owner_sync = owner_sync;
        return sync;
    }

    FieldDescriptor descriptor(const char *name,
                               StaggerLocation location,
                               FieldValueKind value_kind,
                               int ncomp,
                               int nghost,
                               FieldSyncContract sync)
    {
        FieldDescriptor desc;
        desc.name = name;
        desc.location = location;
        desc.value_kind = value_kind;
        desc.ncomp = ncomp;
        desc.nghost = nghost;
        desc.physics = "";
        desc.sync = sync;
        return desc;
    }
}

namespace Z0_NULL
{
    void register_null_fields(Field &fields, int nghost)
    {
#if Z0_ENABLE_EDGE_OWNER
        constexpr OwnerSyncPolicy edge_owner = OwnerSyncPolicy::EdgeOwner;
#else
        constexpr OwnerSyncPolicy edge_owner = OwnerSyncPolicy::None;
#endif

#if Z0_ENABLE_FACE_OWNER
        constexpr OwnerSyncPolicy face_owner = OwnerSyncPolicy::FaceOwner;
#else
        constexpr OwnerSyncPolicy face_owner = OwnerSyncPolicy::None;
#endif

        fields.register_field(descriptor("phi",
                                         StaggerLocation::Cell,
                                         FieldValueKind::Scalar,
                                         1,
                                         nghost,
                                         sync_contract("phi", true)));

        fields.register_field(descriptor("U",
                                         StaggerLocation::Cell,
                                         FieldValueKind::ConservedVector,
                                         5,
                                         nghost,
                                         sync_contract("U", true)));

        fields.register_field(descriptor("V_cart",
                                         StaggerLocation::Cell,
                                         FieldValueKind::CartesianVector,
                                         3,
                                         nghost,
                                         sync_contract("V_cart", true)));

        fields.register_field(descriptor("E_xi",
                                         StaggerLocation::EdgeXi,
                                         FieldValueKind::EdgeCovariant1Form,
                                         1,
                                         nghost,
                                         sync_contract("Eedge", true, HaloLevel::Vertex, true, edge_owner)));
        fields.register_field(descriptor("E_eta",
                                         StaggerLocation::EdgeEt,
                                         FieldValueKind::EdgeCovariant1Form,
                                         1,
                                         nghost,
                                         sync_contract("Eedge", true, HaloLevel::Vertex, true, edge_owner)));
        fields.register_field(descriptor("E_zeta",
                                         StaggerLocation::EdgeZe,
                                         FieldValueKind::EdgeCovariant1Form,
                                         1,
                                         nghost,
                                         sync_contract("Eedge", true, HaloLevel::Vertex, true, edge_owner)));

        fields.register_field(descriptor("B_xi",
                                         StaggerLocation::FaceXi,
                                         FieldValueKind::FaceContravariant2Form,
                                         1,
                                         nghost,
                                         sync_contract("Bface", true, HaloLevel::Vertex, true, face_owner)));
        fields.register_field(descriptor("B_eta",
                                         StaggerLocation::FaceEt,
                                         FieldValueKind::FaceContravariant2Form,
                                         1,
                                         nghost,
                                         sync_contract("Bface", true, HaloLevel::Vertex, true, face_owner)));
        fields.register_field(descriptor("B_zeta",
                                         StaggerLocation::FaceZe,
                                         FieldValueKind::FaceContravariant2Form,
                                         1,
                                         nghost,
                                         sync_contract("Bface", true, HaloLevel::Vertex, true, face_owner)));

        fields.register_field(descriptor("tmp",
                                         StaggerLocation::Cell,
                                         FieldValueKind::Temporary,
                                         1,
                                         nghost,
                                         sync_contract("tmp", false, HaloLevel::FaceOnly)));
    }
}
