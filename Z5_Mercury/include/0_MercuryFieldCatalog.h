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

    FieldDescriptor descriptor(int runtime_nghost) const
    {
        return FieldDescriptor{
            name,
            location,
            ncomp,
            nghost == UseRuntimeGhost ? runtime_nghost : nghost,
            physics ? physics : "",
            sync};
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
                                      HaloLevel halo_level = HaloLevel::Vertex)
{
    FieldSyncContract c;
    c.group = group ? group : "";
    c.do_coupling = do_coupling;
    c.do_physical = do_physical;
    c.do_halo = do_halo;
    c.halo_level = halo_level;
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
         SyncContract("Ucell", false, true, true)},
        {"U_Na", StaggerLocation::Cell, 5, UseRuntimeGhost, "Fluid",
         SyncContract("Ucell", false, true, true)},

        // Magnetic field, electric field, current, and derived cell fields.
        {"B_xi", StaggerLocation::FaceXi, 1, UseRuntimeGhost, "",
         SyncContract("Bface", true, true, true)},
        {"B_eta", StaggerLocation::FaceEt, 1, UseRuntimeGhost, "",
         SyncContract("Bface", true, true, true)},
        {"B_zeta", StaggerLocation::FaceZe, 1, UseRuntimeGhost, "",
         SyncContract("Bface", true, true, true)},
        {"E_xi", StaggerLocation::EdgeXi, 1, UseRuntimeGhost, "",
         SyncContract("Eedge", true, true, true)},
        {"E_eta", StaggerLocation::EdgeEt, 1, UseRuntimeGhost, "",
         SyncContract("Eedge", true, true, true)},
        {"E_zeta", StaggerLocation::EdgeZe, 1, UseRuntimeGhost, "",
         SyncContract("Eedge", true, true, true)},
        {"Ehall_xi", StaggerLocation::EdgeXi, 1, UseRuntimeGhost, "",
         SyncContract("Ehall", true, true, true)},
        {"Ehall_eta", StaggerLocation::EdgeEt, 1, UseRuntimeGhost, "",
         SyncContract("Ehall", true, true, true)},
        {"Ehall_zeta", StaggerLocation::EdgeZe, 1, UseRuntimeGhost, "",
         SyncContract("Ehall", true, true, true)},
        {"Eface_xi", StaggerLocation::FaceXi, 3, UseRuntimeGhost, "",
         SyncContract("Eface", true, true, true)},
        {"Eface_eta", StaggerLocation::FaceEt, 3, UseRuntimeGhost, "",
         SyncContract("Eface", true, true, true)},
        {"Eface_zeta", StaggerLocation::FaceZe, 3, UseRuntimeGhost, "",
         SyncContract("Eface", true, true, true)},
        {"J_xi", StaggerLocation::EdgeXi, 1, UseRuntimeGhost, "",
         SyncContract("Jedge", true, true, true)},
        {"J_eta", StaggerLocation::EdgeEt, 1, UseRuntimeGhost, "",
         SyncContract("Jedge", true, true, true)},
        {"J_zeta", StaggerLocation::EdgeZe, 1, UseRuntimeGhost, "",
         SyncContract("Jedge", true, true, true)},
        {"J_cell", StaggerLocation::Cell, 3, UseRuntimeGhost, "",
         SyncContract("J_cell", true, true, true)},
        {"Badd_xi", StaggerLocation::FaceXi, 1, UseRuntimeGhost, "",
         SyncContract("Badd", true, true, true)},
        {"Badd_eta", StaggerLocation::FaceEt, 1, UseRuntimeGhost, "",
         SyncContract("Badd", true, true, true)},
        {"Badd_zeta", StaggerLocation::FaceZe, 1, UseRuntimeGhost, "",
         SyncContract("Badd", true, true, true)},
        {"B_cell", StaggerLocation::Cell, 3, UseRuntimeGhost, "",
         SyncContract("B_cell", true, true, true)},
        {"Bind_cell", StaggerLocation::Cell, 3, UseRuntimeGhost, "",
         SyncContract("B_cell", true, true, true)},

        // Fluid-side auxiliary fields.
        {"Na", StaggerLocation::Cell, 1, UseRuntimeGhost, "Fluid"},
        {"Photo_rate", StaggerLocation::Cell, 1, UseRuntimeGhost, "Fluid"},
        {"U_plus", StaggerLocation::Cell, 3, UseRuntimeGhost, "Fluid"},
        {"PV_H", StaggerLocation::Cell, 5, UseRuntimeGhost, "Fluid"},
        {"PV_Na", StaggerLocation::Cell, 5, UseRuntimeGhost, "Fluid"},

        // Work fields with no halo contract.
        {"F_xi", StaggerLocation::FaceXi, 5, 0, "Fluid"},
        {"F_eta", StaggerLocation::FaceEt, 5, 0, "Fluid"},
        {"F_zeta", StaggerLocation::FaceZe, 5, 0, "Fluid"},
        {"RHS_H", StaggerLocation::Cell, 5, 0, "Fluid"},
        {"RHS_Na", StaggerLocation::Cell, 5, 0, "Fluid"},
        {"RHS_B_xi", StaggerLocation::FaceXi, 1, 0, ""},
        {"RHS_B_eta", StaggerLocation::FaceEt, 1, 0, ""},
        {"RHS_B_zeta", StaggerLocation::FaceZe, 1, 0, ""},
        {"divB", StaggerLocation::Cell, 1, 1, ""},

        // Hall implicit increment and predictor fields.
        {"dE_xi", StaggerLocation::EdgeXi, 1, UseRuntimeGhost, "",
         SyncContract("dE", true, true, true)},
        {"dE_eta", StaggerLocation::EdgeEt, 1, UseRuntimeGhost, "",
         SyncContract("dE", true, true, true)},
        {"dE_zeta", StaggerLocation::EdgeZe, 1, UseRuntimeGhost, "",
         SyncContract("dE", true, true, true)},
        {"dB_xi", StaggerLocation::FaceXi, 1, UseRuntimeGhost, "",
         SyncContract("dB", true, true, true)},
        {"dB_eta", StaggerLocation::FaceEt, 1, UseRuntimeGhost, "",
         SyncContract("dB", true, true, true)},
        {"dB_zeta", StaggerLocation::FaceZe, 1, UseRuntimeGhost, "",
         SyncContract("dB", true, true, true)},
        {"dJ_xi", StaggerLocation::EdgeXi, 1, UseRuntimeGhost, "",
         SyncContract("dJ", true, true, true)},
        {"dJ_eta", StaggerLocation::EdgeEt, 1, UseRuntimeGhost, "",
         SyncContract("dJ", true, true, true)},
        {"dJ_zeta", StaggerLocation::EdgeZe, 1, UseRuntimeGhost, "",
         SyncContract("dJ", true, true, true)},
        {"dJ_cell", StaggerLocation::Cell, 3, UseRuntimeGhost, "",
         SyncContract("dJcell", true, true, true)},
        {"dEpre_xi", StaggerLocation::EdgeXi, 1, UseRuntimeGhost, "",
         SyncContract("dEpre", true, true, true)},
        {"dEpre_eta", StaggerLocation::EdgeEt, 1, UseRuntimeGhost, "",
         SyncContract("dEpre", true, true, true)},
        {"dEpre_zeta", StaggerLocation::EdgeZe, 1, UseRuntimeGhost, "",
         SyncContract("dEpre", true, true, true)},
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
    std::vector<FieldHaloRequest> requests;
    std::map<std::string, std::size_t> index_by_name;

    for (const auto &group : SyncGroups())
    {
        if (!group.do_halo)
            continue;

        for (const auto &raw_name : group.fields)
        {
            std::string name(raw_name);
            auto it = index_by_name.find(name);
            if (it == index_by_name.end())
            {
                index_by_name[name] = requests.size();
                FieldHaloRequest request;
                request.field_name = name;
                request.sync_group = name;
                request.level = group.halo_level;
                requests.push_back(request);
                continue;
            }

            HaloLevel &old_level = requests[it->second].level;
            if (static_cast<int>(group.halo_level) > static_cast<int>(old_level))
                old_level = group.halo_level;
        }
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
