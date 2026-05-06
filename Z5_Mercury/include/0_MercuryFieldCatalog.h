#pragma once

#include <map>
#include <set>
#include <string>
#include <vector>

#include "3_field/2_MPCNS_Field.h"
#include "4_halo/1_MPCNS_Halo.h"

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

    FieldDescriptor descriptor(int runtime_nghost) const
    {
        return FieldDescriptor{
            name,
            location,
            ncomp,
            nghost == UseRuntimeGhost ? runtime_nghost : nghost,
            physics ? physics : ""};
    }
};

struct SyncGroupSpec
{
    const char *name;
    std::vector<const char *> fields;

    bool do_coupling;
    bool do_physical;
    bool do_halo;
    HaloLevel halo_level;
};

inline std::vector<FieldSpec> FieldSpecs()
{
    return {
        // Conservative fluid state.
        {"U_H", StaggerLocation::Cell, 5, UseRuntimeGhost, "Fluid"},
        {"U_Na", StaggerLocation::Cell, 5, UseRuntimeGhost, "Fluid"},

        // Magnetic field, electric field, current, and derived cell fields.
        {"B_xi", StaggerLocation::FaceXi, 1, UseRuntimeGhost, ""},
        {"B_eta", StaggerLocation::FaceEt, 1, UseRuntimeGhost, ""},
        {"B_zeta", StaggerLocation::FaceZe, 1, UseRuntimeGhost, ""},
        {"E_xi", StaggerLocation::EdgeXi, 1, UseRuntimeGhost, ""},
        {"E_eta", StaggerLocation::EdgeEt, 1, UseRuntimeGhost, ""},
        {"E_zeta", StaggerLocation::EdgeZe, 1, UseRuntimeGhost, ""},
        {"Ehall_xi", StaggerLocation::EdgeXi, 1, UseRuntimeGhost, ""},
        {"Ehall_eta", StaggerLocation::EdgeEt, 1, UseRuntimeGhost, ""},
        {"Ehall_zeta", StaggerLocation::EdgeZe, 1, UseRuntimeGhost, ""},
        {"Eface_xi", StaggerLocation::FaceXi, 3, UseRuntimeGhost, ""},
        {"Eface_eta", StaggerLocation::FaceEt, 3, UseRuntimeGhost, ""},
        {"Eface_zeta", StaggerLocation::FaceZe, 3, UseRuntimeGhost, ""},
        {"J_xi", StaggerLocation::EdgeXi, 1, UseRuntimeGhost, ""},
        {"J_eta", StaggerLocation::EdgeEt, 1, UseRuntimeGhost, ""},
        {"J_zeta", StaggerLocation::EdgeZe, 1, UseRuntimeGhost, ""},
        {"J_cell", StaggerLocation::Cell, 3, UseRuntimeGhost, ""},
        {"Badd_xi", StaggerLocation::FaceXi, 1, UseRuntimeGhost, ""},
        {"Badd_eta", StaggerLocation::FaceEt, 1, UseRuntimeGhost, ""},
        {"Badd_zeta", StaggerLocation::FaceZe, 1, UseRuntimeGhost, ""},
        {"B_cell", StaggerLocation::Cell, 3, UseRuntimeGhost, ""},
        {"Bind_cell", StaggerLocation::Cell, 3, UseRuntimeGhost, ""},

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
        {"dE_xi", StaggerLocation::EdgeXi, 1, UseRuntimeGhost, ""},
        {"dE_eta", StaggerLocation::EdgeEt, 1, UseRuntimeGhost, ""},
        {"dE_zeta", StaggerLocation::EdgeZe, 1, UseRuntimeGhost, ""},
        {"dB_xi", StaggerLocation::FaceXi, 1, UseRuntimeGhost, ""},
        {"dB_eta", StaggerLocation::FaceEt, 1, UseRuntimeGhost, ""},
        {"dB_zeta", StaggerLocation::FaceZe, 1, UseRuntimeGhost, ""},
        {"dJ_xi", StaggerLocation::EdgeXi, 1, UseRuntimeGhost, ""},
        {"dJ_eta", StaggerLocation::EdgeEt, 1, UseRuntimeGhost, ""},
        {"dJ_zeta", StaggerLocation::EdgeZe, 1, UseRuntimeGhost, ""},
        {"dJ_cell", StaggerLocation::Cell, 3, UseRuntimeGhost, ""},
        {"dEpre_xi", StaggerLocation::EdgeXi, 1, UseRuntimeGhost, ""},
        {"dEpre_eta", StaggerLocation::EdgeEt, 1, UseRuntimeGhost, ""},
        {"dEpre_zeta", StaggerLocation::EdgeZe, 1, UseRuntimeGhost, ""},
    };
}

inline std::vector<SyncGroupSpec> SyncGroups()
{
    return {
        {"Ucell", {"U_H", "U_Na"}, false, true, true, HaloLevel::Vertex},
        {"Jedge", {"J_xi", "J_eta", "J_zeta"}, true, true, true, HaloLevel::Vertex},
        {"Eedge", {"E_xi", "E_eta", "E_zeta"}, true, true, true, HaloLevel::Vertex},
        {"Ehall", {"Ehall_xi", "Ehall_eta", "Ehall_zeta"}, true, true, true, HaloLevel::Vertex},
        {"Bface", {"B_xi", "B_eta", "B_zeta"}, true, true, true, HaloLevel::Vertex},
        {"Eface", {"Eface_xi", "Eface_eta", "Eface_zeta"}, true, true, true, HaloLevel::Vertex},
        {"B_cell", {"B_cell", "Bind_cell"}, true, true, true, HaloLevel::Vertex},
        {"J_cell", {"J_cell"}, true, true, true, HaloLevel::Vertex},
        {"dE", {"dE_xi", "dE_eta", "dE_zeta"}, true, true, true, HaloLevel::Vertex},
        {"dB", {"dB_xi", "dB_eta", "dB_zeta"}, true, true, true, HaloLevel::Vertex},
        {"dJ", {"dJ_xi", "dJ_eta", "dJ_zeta"}, true, true, true, HaloLevel::Vertex},
        {"dJcell", {"dJ_cell"}, true, true, true, HaloLevel::Vertex},
        {"dEpre", {"dEpre_xi", "dEpre_eta", "dEpre_zeta"}, true, true, true, HaloLevel::Vertex},
        {"Badd", {"Badd_xi", "Badd_eta", "Badd_zeta"}, true, true, true, HaloLevel::Vertex},
    };
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
        for (const char *name : group.fields)
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
        for (const char *name : group.fields)
            append_unique(names, seen, name);
    }
    return names;
}

inline std::vector<HaloFieldRequest> HaloRequests()
{
    std::vector<HaloFieldRequest> requests;
    std::map<std::string, std::size_t> index_by_name;

    for (const auto &group : SyncGroups())
    {
        if (!group.do_halo)
            continue;

        for (const char *raw_name : group.fields)
        {
            std::string name(raw_name);
            auto it = index_by_name.find(name);
            if (it == index_by_name.end())
            {
                index_by_name[name] = requests.size();
                requests.push_back(HaloFieldRequest{name, group.halo_level});
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
    const auto names = CoupledFieldNames();
    for (const auto &name : names)
    {
        fld->register_coupling_channel("Solid", "Fluid", name);
        fld->register_coupling_channel("Fluid", "Solid", name);
    }

    fld->build_coupling_buffers(topology, dimension);
}

inline void RegisterHaloFields(Halo *halo)
{
    for (const auto &request : HaloRequests())
        halo->register_halo_field(request.field_name, request.level);

    halo->build_registered_patterns();
}

} // namespace MERCURY_FIELD
