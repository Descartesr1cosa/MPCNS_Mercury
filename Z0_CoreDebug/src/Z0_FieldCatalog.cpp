#include "Z0_FieldCatalog.h"

#include "3_field/2_MPCNS_Field.h"

namespace
{
    FieldSyncContract sync_contract(const char *group,
                                    bool do_halo,
                                    bool do_coupling,
                                    HaloLevel halo_level = HaloLevel::Vertex,
                                    bool orientation_aware = false,
                                    OwnerSyncPolicy owner_sync = OwnerSyncPolicy::None)
    {
        FieldSyncContract sync;
        sync.group = group ? group : "";
        sync.do_halo = do_halo;
        sync.do_coupling = do_coupling;
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

namespace Z0
{
    void register_core_debug_fields(Field &fields, int nghost)
    {
        fields.register_field(descriptor("phi_cell", StaggerLocation::Cell, FieldValueKind::Scalar,
                                         1, nghost, sync_contract("phi_cell", true, true)));
        fields.register_field(descriptor("U_cell", StaggerLocation::Cell, FieldValueKind::ConservedVector,
                                         5, nghost, sync_contract("U_cell", true, true)));
        fields.register_field(descriptor("V_cell", StaggerLocation::Cell, FieldValueKind::CartesianVector,
                                         3, nghost, sync_contract("V_cell", true, true)));
        fields.register_field(descriptor("psi_node", StaggerLocation::Node, FieldValueKind::Scalar,
                                         1, nghost, sync_contract("psi_node", true, true)));

        fields.register_field(descriptor("B_xi", StaggerLocation::FaceXi, FieldValueKind::FaceContravariant2Form,
                                         1, nghost, sync_contract("Bface", true, true, HaloLevel::Vertex, true, OwnerSyncPolicy::FaceOwner)));
        fields.register_field(descriptor("B_eta", StaggerLocation::FaceEt, FieldValueKind::FaceContravariant2Form,
                                         1, nghost, sync_contract("Bface", true, true, HaloLevel::Vertex, true, OwnerSyncPolicy::FaceOwner)));
        fields.register_field(descriptor("B_zeta", StaggerLocation::FaceZe, FieldValueKind::FaceContravariant2Form,
                                         1, nghost, sync_contract("Bface", true, true, HaloLevel::Vertex, true, OwnerSyncPolicy::FaceOwner)));

        fields.register_field(descriptor("E_xi", StaggerLocation::EdgeXi, FieldValueKind::EdgeCovariant1Form,
                                         1, nghost, sync_contract("Eedge", true, true, HaloLevel::Vertex, true, OwnerSyncPolicy::EdgeOwner)));
        fields.register_field(descriptor("E_eta", StaggerLocation::EdgeEt, FieldValueKind::EdgeCovariant1Form,
                                         1, nghost, sync_contract("Eedge", true, true, HaloLevel::Vertex, true, OwnerSyncPolicy::EdgeOwner)));
        fields.register_field(descriptor("E_zeta", StaggerLocation::EdgeZe, FieldValueKind::EdgeCovariant1Form,
                                         1, nghost, sync_contract("Eedge", true, true, HaloLevel::Vertex, true, OwnerSyncPolicy::EdgeOwner)));

        fields.register_field(descriptor("FaceXi_cart", StaggerLocation::FaceXi, FieldValueKind::CartesianVector,
                                         3, nghost, sync_contract("FaceXi_cart", true, true)));
        fields.register_field(descriptor("FaceEt_cart", StaggerLocation::FaceEt, FieldValueKind::CartesianVector,
                                         3, nghost, sync_contract("FaceEt_cart", true, true)));
        fields.register_field(descriptor("FaceZe_cart", StaggerLocation::FaceZe, FieldValueKind::CartesianVector,
                                         3, nghost, sync_contract("FaceZe_cart", true, true)));
        fields.register_field(descriptor("EdgeXi_cart", StaggerLocation::EdgeXi, FieldValueKind::CartesianVector,
                                         3, nghost, sync_contract("EdgeXi_cart", true, true)));
        fields.register_field(descriptor("EdgeEt_cart", StaggerLocation::EdgeEt, FieldValueKind::CartesianVector,
                                         3, nghost, sync_contract("EdgeEt_cart", true, true)));
        fields.register_field(descriptor("EdgeZe_cart", StaggerLocation::EdgeZe, FieldValueKind::CartesianVector,
                                         3, nghost, sync_contract("EdgeZe_cart", true, true)));

        fields.register_field(descriptor("tmp_cell", StaggerLocation::Cell, FieldValueKind::Temporary,
                                         1, nghost, sync_contract("tmp", false, false, HaloLevel::FaceOnly)));
    }
}
