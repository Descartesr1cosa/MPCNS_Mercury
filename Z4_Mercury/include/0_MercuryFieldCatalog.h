#pragma once

#include <cstddef>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "3_field/Field.h"
#include "4_halo/Halo.h"

namespace MERCURY_FIELD
{
constexpr int UseRuntimeGhost = -1;

struct FieldSpec
{
    const char *name;
    StaggerLocation location;
    int ncomp;
    int nghost; // UseRuntimeGhost means the runtime parameter ngg.
    const char *physics;
    FieldSyncContract sync;
    FieldValueKind value_kind = FieldValueKind::Scalar;

    FieldDescriptor descriptor(int runtime_nghost) const
    {
        FieldDescriptor desc;
        desc.name = name;
        desc.location = location;
        desc.ncomp = ncomp;
        desc.nghost = nghost == UseRuntimeGhost ? runtime_nghost : nghost;
        desc.physics = physics ? physics : "";
        desc.sync = sync;
        desc.value_kind = value_kind;
        return desc;
    }
};

struct SyncGroupSpec
{
    std::string name;
    std::vector<std::string> fields;

    bool do_coupling;
    bool do_physical;
    bool do_halo;
    HaloLevel halo_level;
};

inline FieldSyncContract SyncContract(const char *group,
                                      bool do_coupling,
                                      bool do_physical,
                                      bool do_halo,
                                      HaloLevel halo_level = HaloLevel::Corner3D,
                                      OwnerSyncPolicy owner_sync = OwnerSyncPolicy::None,
                                      bool orientation_aware = false)
{
    FieldSyncContract c;
    c.group = group ? group : "";
    c.do_coupling = do_coupling;
    c.do_physical = do_physical;
    c.do_halo = do_halo;
    c.halo_level = halo_level;
    c.owner_sync = owner_sync;
    c.orientation_aware = orientation_aware;
    return c;
}

inline std::vector<FieldSpec> FieldSpecs()
{
    // Mercury field single source of truth:
    // name/location/ncomp/nghost/physics decide allocation;
    // SyncContract decides runtime group, coupling channel, physical BC, and halo.
    return {
        // Conservative fluid state.
        {"U_H", StaggerLocation::Cell, 5, UseRuntimeGhost, "Fluid",
         SyncContract("Ucell", false, true, true, HaloLevel::Corner1D),
         FieldValueKind::ConservedVector},
        {"U_Na", StaggerLocation::Cell, 5, UseRuntimeGhost, "Fluid",
         SyncContract("Ucell", false, true, true, HaloLevel::Corner1D),
         FieldValueKind::ConservedVector},

        // Magnetic field, electric field, current, and derived cell fields.
        {"B_xi", StaggerLocation::FaceXi, 1, UseRuntimeGhost, "",
         SyncContract("Bface", true, true, true, HaloLevel::Corner3D, OwnerSyncPolicy::FaceOwner, true),
         FieldValueKind::FaceContravariant2Form},
        {"B_eta", StaggerLocation::FaceEt, 1, UseRuntimeGhost, "",
         SyncContract("Bface", true, true, true, HaloLevel::Corner3D, OwnerSyncPolicy::FaceOwner, true),
         FieldValueKind::FaceContravariant2Form},
        {"B_zeta", StaggerLocation::FaceZe, 1, UseRuntimeGhost, "",
         SyncContract("Bface", true, true, true, HaloLevel::Corner3D, OwnerSyncPolicy::FaceOwner, true),
         FieldValueKind::FaceContravariant2Form},
        {"E_xi", StaggerLocation::EdgeXi, 1, UseRuntimeGhost, "",
         SyncContract("Eedge", true, true, true, HaloLevel::Corner3D, OwnerSyncPolicy::EdgeOwner, true),
         FieldValueKind::EdgeCovariant1Form},
        {"E_eta", StaggerLocation::EdgeEt, 1, UseRuntimeGhost, "",
         SyncContract("Eedge", true, true, true, HaloLevel::Corner3D, OwnerSyncPolicy::EdgeOwner, true),
         FieldValueKind::EdgeCovariant1Form},
        {"E_zeta", StaggerLocation::EdgeZe, 1, UseRuntimeGhost, "",
         SyncContract("Eedge", true, true, true, HaloLevel::Corner3D, OwnerSyncPolicy::EdgeOwner, true),
         FieldValueKind::EdgeCovariant1Form},
        {"Ehall_xi", StaggerLocation::EdgeXi, 1, UseRuntimeGhost, "",
         SyncContract("Ehall", true, true, true, HaloLevel::Corner3D, OwnerSyncPolicy::EdgeOwner, true),
         FieldValueKind::EdgeCovariant1Form},
        {"Ehall_eta", StaggerLocation::EdgeEt, 1, UseRuntimeGhost, "",
         SyncContract("Ehall", true, true, true, HaloLevel::Corner3D, OwnerSyncPolicy::EdgeOwner, true),
         FieldValueKind::EdgeCovariant1Form},
        {"Ehall_zeta", StaggerLocation::EdgeZe, 1, UseRuntimeGhost, "",
         SyncContract("Ehall", true, true, true, HaloLevel::Corner3D, OwnerSyncPolicy::EdgeOwner, true),
         FieldValueKind::EdgeCovariant1Form},
        {"Eres_xi", StaggerLocation::EdgeXi, 1, UseRuntimeGhost, "",
         SyncContract("Eres", true, true, true, HaloLevel::Corner3D, OwnerSyncPolicy::EdgeOwner, true),
         FieldValueKind::EdgeCovariant1Form},
        {"Eres_eta", StaggerLocation::EdgeEt, 1, UseRuntimeGhost, "",
         SyncContract("Eres", true, true, true, HaloLevel::Corner3D, OwnerSyncPolicy::EdgeOwner, true),
         FieldValueKind::EdgeCovariant1Form},
        {"Eres_zeta", StaggerLocation::EdgeZe, 1, UseRuntimeGhost, "",
         SyncContract("Eres", true, true, true, HaloLevel::Corner3D, OwnerSyncPolicy::EdgeOwner, true),
         FieldValueKind::EdgeCovariant1Form},
        {"Eface_xi", StaggerLocation::FaceXi, 3, UseRuntimeGhost, "",
         SyncContract("Eface", true, true, true),
         FieldValueKind::CartesianVector},
        {"Eface_eta", StaggerLocation::FaceEt, 3, UseRuntimeGhost, "",
         SyncContract("Eface", true, true, true),
         FieldValueKind::CartesianVector},
        {"Eface_zeta", StaggerLocation::FaceZe, 3, UseRuntimeGhost, "",
         SyncContract("Eface", true, true, true),
         FieldValueKind::CartesianVector},
        {"J_xi", StaggerLocation::EdgeXi, 1, UseRuntimeGhost, "",
         SyncContract("Jedge", true, true, true, HaloLevel::Corner3D, OwnerSyncPolicy::EdgeOwner, true),
         FieldValueKind::EdgeCovariant1Form},
        {"J_eta", StaggerLocation::EdgeEt, 1, UseRuntimeGhost, "",
         SyncContract("Jedge", true, true, true, HaloLevel::Corner3D, OwnerSyncPolicy::EdgeOwner, true),
         FieldValueKind::EdgeCovariant1Form},
        {"J_zeta", StaggerLocation::EdgeZe, 1, UseRuntimeGhost, "",
         SyncContract("Jedge", true, true, true, HaloLevel::Corner3D, OwnerSyncPolicy::EdgeOwner, true),
         FieldValueKind::EdgeCovariant1Form},
        {"J_cell", StaggerLocation::Cell, 3, UseRuntimeGhost, "",
         SyncContract("J_cell", true, true, true),
         FieldValueKind::CartesianVector},
        {"Badd_xi", StaggerLocation::FaceXi, 1, UseRuntimeGhost, "",
         SyncContract("Badd", true, true, true, HaloLevel::Corner3D, OwnerSyncPolicy::FaceOwner, true),
         FieldValueKind::FaceContravariant2Form},
        {"Badd_eta", StaggerLocation::FaceEt, 1, UseRuntimeGhost, "",
         SyncContract("Badd", true, true, true, HaloLevel::Corner3D, OwnerSyncPolicy::FaceOwner, true),
         FieldValueKind::FaceContravariant2Form},
        {"Badd_zeta", StaggerLocation::FaceZe, 1, UseRuntimeGhost, "",
         SyncContract("Badd", true, true, true, HaloLevel::Corner3D, OwnerSyncPolicy::FaceOwner, true),
         FieldValueKind::FaceContravariant2Form},
        {"B_cell", StaggerLocation::Cell, 3, UseRuntimeGhost, "",
         SyncContract("B_cell", true, true, true),
         FieldValueKind::CartesianVector},
        {"Bind_cell", StaggerLocation::Cell, 3, UseRuntimeGhost, "",
         SyncContract("B_cell", true, true, true),
         FieldValueKind::CartesianVector},

        // Fluid-side auxiliary fields.
        {"Na", StaggerLocation::Cell, 1, UseRuntimeGhost, "Fluid"},
        {"Photo_rate", StaggerLocation::Cell, 1, UseRuntimeGhost, "Fluid"},
        {"U_plus", StaggerLocation::Cell, 3, UseRuntimeGhost, "Fluid",
         SyncContract("Uplus", false, false, true), FieldValueKind::CartesianVector},
        {"PV_H", StaggerLocation::Cell, 5, UseRuntimeGhost, "Fluid"},
        {"PV_Na", StaggerLocation::Cell, 5, UseRuntimeGhost, "Fluid"},
        {"Bcell_from_Bface_w", StaggerLocation::Cell, 18, UseRuntimeGhost, ""},
        {"Jcell_from_Jedge_w", StaggerLocation::Cell, 36, UseRuntimeGhost, ""},

        // Work fields with no halo contract.
        {"F_xi", StaggerLocation::FaceXi, 5, 0, "Fluid"},
        {"F_eta", StaggerLocation::FaceEt, 5, 0, "Fluid"},
        {"F_zeta", StaggerLocation::FaceZe, 5, 0, "Fluid"},
        {"RHS_H", StaggerLocation::Cell, 5, 0, "Fluid"},
        {"RHS_Na", StaggerLocation::Cell, 5, 0, "Fluid"},
        {"RHS_B_xi", StaggerLocation::FaceXi, 1, 0, ""},
        {"RHS_B_eta", StaggerLocation::FaceEt, 1, 0, ""},
        {"RHS_B_zeta", StaggerLocation::FaceZe, 1, 0, ""},
        {"RHS_Bres_xi", StaggerLocation::FaceXi, 1, 0, ""},
        {"RHS_Bres_eta", StaggerLocation::FaceEt, 1, 0, ""},
        {"RHS_Bres_zeta", StaggerLocation::FaceZe, 1, 0, ""},
        {"divB", StaggerLocation::Cell, 1, 1, ""},

        // Hall implicit increment and predictor fields.
        {"dE_xi", StaggerLocation::EdgeXi, 1, UseRuntimeGhost, "",
         SyncContract("dE", true, true, true, HaloLevel::Corner3D, OwnerSyncPolicy::EdgeOwner, true),
         FieldValueKind::EdgeCovariant1Form},
        {"dE_eta", StaggerLocation::EdgeEt, 1, UseRuntimeGhost, "",
         SyncContract("dE", true, true, true, HaloLevel::Corner3D, OwnerSyncPolicy::EdgeOwner, true),
         FieldValueKind::EdgeCovariant1Form},
        {"dE_zeta", StaggerLocation::EdgeZe, 1, UseRuntimeGhost, "",
         SyncContract("dE", true, true, true, HaloLevel::Corner3D, OwnerSyncPolicy::EdgeOwner, true),
         FieldValueKind::EdgeCovariant1Form},
        {"dB_xi", StaggerLocation::FaceXi, 1, UseRuntimeGhost, "",
         SyncContract("dB", true, true, true, HaloLevel::Corner3D, OwnerSyncPolicy::FaceOwner, true),
         FieldValueKind::FaceContravariant2Form},
        {"dB_eta", StaggerLocation::FaceEt, 1, UseRuntimeGhost, "",
         SyncContract("dB", true, true, true, HaloLevel::Corner3D, OwnerSyncPolicy::FaceOwner, true),
         FieldValueKind::FaceContravariant2Form},
        {"dB_zeta", StaggerLocation::FaceZe, 1, UseRuntimeGhost, "",
         SyncContract("dB", true, true, true, HaloLevel::Corner3D, OwnerSyncPolicy::FaceOwner, true),
         FieldValueKind::FaceContravariant2Form},
        {"dJ_xi", StaggerLocation::EdgeXi, 1, UseRuntimeGhost, "",
         SyncContract("dJ", true, true, true, HaloLevel::Corner3D, OwnerSyncPolicy::EdgeOwner, true),
         FieldValueKind::EdgeCovariant1Form},
        {"dJ_eta", StaggerLocation::EdgeEt, 1, UseRuntimeGhost, "",
         SyncContract("dJ", true, true, true, HaloLevel::Corner3D, OwnerSyncPolicy::EdgeOwner, true),
         FieldValueKind::EdgeCovariant1Form},
        {"dJ_zeta", StaggerLocation::EdgeZe, 1, UseRuntimeGhost, "",
         SyncContract("dJ", true, true, true, HaloLevel::Corner3D, OwnerSyncPolicy::EdgeOwner, true),
         FieldValueKind::EdgeCovariant1Form},
        {"dJ_cell", StaggerLocation::Cell, 3, UseRuntimeGhost, "",
         SyncContract("dJcell", true, true, true),
         FieldValueKind::CartesianVector},
        {"dEpre_xi", StaggerLocation::EdgeXi, 1, UseRuntimeGhost, "",
         SyncContract("dEpre", true, true, true, HaloLevel::Corner3D, OwnerSyncPolicy::EdgeOwner, true),
         FieldValueKind::EdgeCovariant1Form},
        {"dEpre_eta", StaggerLocation::EdgeEt, 1, UseRuntimeGhost, "",
         SyncContract("dEpre", true, true, true, HaloLevel::Corner3D, OwnerSyncPolicy::EdgeOwner, true),
         FieldValueKind::EdgeCovariant1Form},
        {"dEpre_zeta", StaggerLocation::EdgeZe, 1, UseRuntimeGhost, "",
         SyncContract("dEpre", true, true, true, HaloLevel::Corner3D, OwnerSyncPolicy::EdgeOwner, true),
         FieldValueKind::EdgeCovariant1Form},
    };
}

inline std::vector<SyncGroupSpec> SyncGroups()
{
    std::vector<SyncGroupSpec> groups;
    std::map<std::string, std::size_t> index_by_name;

    for (const auto &spec : FieldSpecs())
    {
        const FieldSyncContract &sync = spec.sync;
        if (sync.group.empty())
            continue;

        auto it = index_by_name.find(sync.group);
        if (it == index_by_name.end())
        {
            SyncGroupSpec group;
            group.name = sync.group;
            group.do_coupling = sync.do_coupling;
            group.do_physical = sync.do_physical;
            group.do_halo = sync.do_halo;
            group.halo_level = sync.halo_level;
            group.fields.push_back(spec.name);

            index_by_name[group.name] = groups.size();
            groups.push_back(group);
            continue;
        }

        SyncGroupSpec &group = groups[it->second];
        group.fields.push_back(spec.name);
        group.do_coupling = group.do_coupling || sync.do_coupling;
        group.do_physical = group.do_physical || sync.do_physical;
        group.do_halo = group.do_halo || sync.do_halo;
        if (static_cast<int>(sync.halo_level) > static_cast<int>(group.halo_level))
            group.halo_level = sync.halo_level;
    }

    return groups;
}

inline void append_unique(std::vector<std::string> &out,
                          std::set<std::string> &seen,
                          const std::string &name)
{
    if (seen.insert(name).second)
        out.push_back(name);
}

inline std::vector<std::string> BoundaryFieldNames()
{
    std::vector<std::string> names;
    std::set<std::string> seen;
    for (const auto &group : SyncGroups())
    {
        if (!group.do_physical)
            continue;
        for (const auto &name : group.fields)
            append_unique(names, seen, name);
    }
    return names;
}

inline std::vector<std::string> CoupledFieldNames()
{
    std::vector<std::string> names;
    std::set<std::string> seen;
    for (const auto &group : SyncGroups())
    {
        if (!group.do_coupling)
            continue;
        for (const auto &name : group.fields)
            append_unique(names, seen, name);
    }
    return names;
}

inline std::vector<FieldHaloRequest> HaloRequests()
{
    // Compatibility helper. Runtime code should prefer Field::halo_requests(),
    // because it has the actual runtime nghost value after registration.
    std::vector<FieldHaloRequest> requests;
    for (const auto &spec : FieldSpecs())
    {
        if (!spec.sync.do_halo &&
            spec.sync.owner_sync == OwnerSyncPolicy::None)
            continue;

        FieldHaloRequest request;
        request.field_name = spec.name;
        request.sync_group = spec.sync.group.empty() ? std::string(spec.name) : spec.sync.group;
        request.location = spec.location;
        request.value_kind = spec.value_kind;
        request.ncomp = spec.ncomp;
        request.nghost = spec.nghost == UseRuntimeGhost ? 0 : spec.nghost;
        request.do_halo = spec.sync.do_halo;
        request.level = spec.sync.halo_level;
        request.owner_sync = spec.sync.owner_sync;
        request.orientation_aware = spec.sync.orientation_aware;
        requests.push_back(request);
    }
    return requests;
}

inline void RegisterFields(Field *fld, int runtime_nghost)
{
    for (const auto &spec : FieldSpecs())
        fld->register_field(spec.descriptor(runtime_nghost));
}

inline void RegisterCouplingChannels(Field *fld,
                                     const TOPO::Topology &topology,
                                     int dimension)
{
    const std::vector<Field::PairKey> directed_pairs = {{"Solid", "Fluid"}, {"Fluid", "Solid"}};
    fld->register_declared_coupling_channels(directed_pairs);

    fld->build_coupling_buffers(topology, dimension);
}

inline void RegisterHaloFields(Field *fld, Halo *halo)
{
    halo->register_halo_fields(fld->halo_requests());
    halo->build_registered_patterns();
}

} // namespace MERCURY_FIELD
