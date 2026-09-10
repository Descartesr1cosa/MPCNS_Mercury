#pragma once

#include "MarsSpecies.h"
#include "3_field/Field.h"
#include "4_halo/Halo.h"

#include <string>
#include <vector>

namespace MARS
{
inline FieldSyncContract Sync(const char *group, bool physical, bool halo,
                              HaloLevel level = HaloLevel::Corner3D,
                              OwnerSyncPolicy owner = OwnerSyncPolicy::None,
                              bool oriented = false)
{
    FieldSyncContract c;
    c.group = group;
    c.do_physical = physical;
    c.do_halo = halo;
    c.halo_level = level;
    c.owner_sync = owner;
    c.orientation_aware = oriented;
    return c;
}

inline FieldDescriptor Descriptor(std::string name, StaggerLocation location,
                                  int ncomp, int nghost, std::string physics,
                                  FieldSyncContract sync = {},
                                  FieldValueKind kind = FieldValueKind::Scalar)
{
    FieldDescriptor d;
    d.name = std::move(name);
    d.location = location;
    d.ncomp = ncomp;
    d.nghost = nghost;
    d.physics = std::move(physics);
    d.sync = std::move(sync);
    d.value_kind = kind;
    return d;
}

inline void RegisterMarsFields(Field &field, int ngh)
{
    for (std::string_view species : SpeciesNames)
    {
        const std::string suffix(species);
        field.register_field(Descriptor("U_" + suffix, StaggerLocation::Cell, 5, ngh, "Fluid",
                                        Sync("Ucell", true, true, HaloLevel::Corner1D),
                                        FieldValueKind::ConservedVector));
        field.register_field(Descriptor("PV_" + suffix, StaggerLocation::Cell, 5, ngh, "Fluid"));
        field.register_field(Descriptor("RHS_" + suffix, StaggerLocation::Cell, 5, 0, "Fluid"));
    }

    field.register_field(Descriptor("neutral_density", StaggerLocation::Cell,
                                    static_cast<int>(NumNeutral), ngh, "Fluid"));
    field.register_field(Descriptor("photo_rate", StaggerLocation::Cell,
                                    static_cast<int>(NumPhoto), ngh, "Fluid"));
    field.register_field(Descriptor("U_plus", StaggerLocation::Cell, 3, ngh, "Fluid",
                                    Sync("Uplus", true, true)));
    field.register_field(Descriptor("B_cell", StaggerLocation::Cell, 3, ngh, "Fluid",
                                    Sync("Bcell", true, true), FieldValueKind::CartesianVector));
    field.register_field(Descriptor("Bind_cell", StaggerLocation::Cell, 3, ngh, "Fluid",
                                    Sync("Bcell", true, true), FieldValueKind::CartesianVector));
    field.register_field(Descriptor("J_cell", StaggerLocation::Cell, 3, ngh, "Fluid",
                                    Sync("Jcell", true, true), FieldValueKind::CartesianVector));
    field.register_field(Descriptor("divB", StaggerLocation::Cell, 1, 1, "Fluid"));
    field.register_field(Descriptor("divB_total", StaggerLocation::Cell, 1, 1, "Fluid"));

    const char *suffix[3] = {"xi", "eta", "zeta"};
    const StaggerLocation face_loc[3] = {StaggerLocation::FaceXi, StaggerLocation::FaceEt,
                                         StaggerLocation::FaceZe};
    const StaggerLocation edge_loc[3] = {StaggerLocation::EdgeXi, StaggerLocation::EdgeEt,
                                         StaggerLocation::EdgeZe};
    for (int d = 0; d < 3; ++d)
    {
        // Badd = prescribed IMF + crustal flux; Aadd is crustal potential only.
        // No physical handler may overwrite an active Badd face after curl(A).
        field.register_field(Descriptor(std::string("Badd_") + suffix[d], face_loc[d], 1, ngh, "Fluid",
                                        Sync("Badd", false, true, HaloLevel::Corner3D,
                                             OwnerSyncPolicy::FaceOwner, true),
                                        FieldValueKind::FaceContravariant2Form));
        field.register_field(Descriptor(std::string("Aadd_") + suffix[d], edge_loc[d], 1, ngh, "Fluid",
                                        Sync("Aadd", false, true, HaloLevel::Corner3D,
                                             OwnerSyncPolicy::EdgeOwner, true),
                                        FieldValueKind::EdgeCovariant1Form));
        field.register_field(Descriptor(std::string("B_") + suffix[d], face_loc[d], 1, ngh, "Fluid",
                                        Sync("Bface", true, true, HaloLevel::Corner3D,
                                             OwnerSyncPolicy::FaceOwner, true),
                                        FieldValueKind::FaceContravariant2Form));
        field.register_field(Descriptor(std::string("E_") + suffix[d], edge_loc[d], 1, ngh, "Fluid",
                                        Sync("Eedge", true, true, HaloLevel::Corner3D,
                                             OwnerSyncPolicy::EdgeOwner, true),
                                        FieldValueKind::EdgeCovariant1Form));
        // Mercury-compatible ideal-UCT work tensor.  Its three components are
        // block-local covariant EMF candidates, so it must not be exchanged
        // across rotated block interfaces.
        field.register_field(Descriptor(std::string("Eface_") + suffix[d], face_loc[d], 3, ngh,
                                        "Fluid"));
        // Keep Hall separate until its own oriented edge synchronization has
        // completed, matching Mercury's induction assembly order.
        field.register_field(Descriptor(std::string("Ehall_") + suffix[d], edge_loc[d], 1, ngh,
                                        "Fluid",
                                        Sync("Ehall", true, true, HaloLevel::Corner3D,
                                             OwnerSyncPolicy::EdgeOwner, true),
                                        FieldValueKind::EdgeCovariant1Form));
        field.register_field(Descriptor(std::string("RHS_B_") + suffix[d], face_loc[d], 1, 0, "Fluid"));
        field.register_field(Descriptor(std::string("F_") + suffix[d], face_loc[d], 5, 0, "Fluid"));
    }

    // Cell scratch used by Mercury's Hall Rusanov edge reconstruction.  During
    // a step it first caches (ne_true,ne_eff,Pe), then Hall overwrites it with
    // its Cartesian cell electric field before output.
    field.register_field(Descriptor("E_hall_cell", StaggerLocation::Cell, 3, ngh, "Fluid"));
    field.register_field(Descriptor("Hall_beta", StaggerLocation::Cell, 1, ngh, "Fluid"));
}

inline void RegisterMarsHalo(Field &field, Halo &halo)
{
    halo.register_halo_fields(field.halo_requests());
    halo.build_registered_patterns();
}
}
