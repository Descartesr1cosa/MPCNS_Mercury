#include "MercurySolver.h"

#include <array>
#include <string>

#include "2_topology/2_MPCNS_Topology.h"
#include "3_field/2_MPCNS_Field.h"
#include "4_halo/1_MPCNS_Halo.h"

namespace
{
struct CouplingField
{
    const char *name;
    StaggerLocation location;
    int ncomp;
};

void RegisterCouplingPair(Field *fld, const char *src, const char *dst,
                          const CouplingField *fields, int nfields, int ngg)
{
    for (int n = 0; n < nfields; ++n)
    {
        const auto &f = fields[n];
        fld->register_coupling_channel(src, dst, f.name, f.location, f.ncomp, ngg);
    }
}
} // namespace

void MercurySolver::RegisterFields(Field *fld, int ngg)
{
    const std::array<FieldDescriptor, 28> fields = {{
        {"U_H", StaggerLocation::Cell, 5, ngg, "Fluid"},
        {"U_Na", StaggerLocation::Cell, 5, ngg, "Fluid"},
        {"B_xi", StaggerLocation::FaceXi, 1, ngg, ""},
        {"B_eta", StaggerLocation::FaceEt, 1, ngg, ""},
        {"B_zeta", StaggerLocation::FaceZe, 1, ngg, ""},
        {"E_xi", StaggerLocation::EdgeXi, 1, ngg, ""},
        {"E_eta", StaggerLocation::EdgeEt, 1, ngg, ""},
        {"E_zeta", StaggerLocation::EdgeZe, 1, ngg, ""},
        {"Ehall_xi", StaggerLocation::EdgeXi, 1, ngg, ""},
        {"Ehall_eta", StaggerLocation::EdgeEt, 1, ngg, ""},
        {"Ehall_zeta", StaggerLocation::EdgeZe, 1, ngg, ""},
        {"Eface_xi", StaggerLocation::FaceXi, 3, ngg, ""},
        {"Eface_eta", StaggerLocation::FaceEt, 3, ngg, ""},
        {"Eface_zeta", StaggerLocation::FaceZe, 3, ngg, ""},
        {"J_xi", StaggerLocation::EdgeXi, 1, ngg, ""},
        {"J_eta", StaggerLocation::EdgeEt, 1, ngg, ""},
        {"J_zeta", StaggerLocation::EdgeZe, 1, ngg, ""},
        {"J_cell", StaggerLocation::Cell, 3, ngg, ""},
        {"Badd_xi", StaggerLocation::FaceXi, 1, ngg, ""},
        {"Badd_eta", StaggerLocation::FaceEt, 1, ngg, ""},
        {"Badd_zeta", StaggerLocation::FaceZe, 1, ngg, ""},
        {"B_cell", StaggerLocation::Cell, 3, ngg, ""},
        {"Bind_cell", StaggerLocation::Cell, 3, ngg, ""},
        {"Na", StaggerLocation::Cell, 1, ngg, "Fluid"},
        {"Photo_rate", StaggerLocation::Cell, 1, ngg, "Fluid"},
        {"U_plus", StaggerLocation::Cell, 3, ngg, "Fluid"},
        {"PV_H", StaggerLocation::Cell, 5, ngg, "Fluid"},
        {"PV_Na", StaggerLocation::Cell, 5, ngg, "Fluid"},
    }};

    for (const auto &field : fields)
        fld->register_field(field);

    const std::array<FieldDescriptor, 6> work_fields = {{
        {"F_xi", StaggerLocation::FaceXi, 5, 0, "Fluid"},
        {"F_eta", StaggerLocation::FaceEt, 5, 0, "Fluid"},
        {"F_zeta", StaggerLocation::FaceZe, 5, 0, "Fluid"},
        {"RHS_H", StaggerLocation::Cell, 5, 0, "Fluid"},
        {"RHS_Na", StaggerLocation::Cell, 5, 0, "Fluid"},
        {"divB", StaggerLocation::Cell, 1, 1, ""},
    }};

    for (const auto &field : work_fields)
        fld->register_field(field);

    const std::array<FieldDescriptor, 3> magnetic_rhs = {{
        {"RHS_B_xi", StaggerLocation::FaceXi, 1, 0, ""},
        {"RHS_B_eta", StaggerLocation::FaceEt, 1, 0, ""},
        {"RHS_B_zeta", StaggerLocation::FaceZe, 1, 0, ""},
    }};

    for (const auto &field : magnetic_rhs)
        fld->register_field(field);

    const std::array<FieldDescriptor, 13> hall_implicit_fields = {{
        {"dE_xi", StaggerLocation::EdgeXi, 1, ngg, ""},
        {"dE_eta", StaggerLocation::EdgeEt, 1, ngg, ""},
        {"dE_zeta", StaggerLocation::EdgeZe, 1, ngg, ""},
        {"dB_xi", StaggerLocation::FaceXi, 1, ngg, ""},
        {"dB_eta", StaggerLocation::FaceEt, 1, ngg, ""},
        {"dB_zeta", StaggerLocation::FaceZe, 1, ngg, ""},
        {"dJ_xi", StaggerLocation::EdgeXi, 1, ngg, ""},
        {"dJ_eta", StaggerLocation::EdgeEt, 1, ngg, ""},
        {"dJ_zeta", StaggerLocation::EdgeZe, 1, ngg, ""},
        {"dJ_cell", StaggerLocation::Cell, 3, ngg, ""},
        {"dEpre_xi", StaggerLocation::EdgeXi, 1, ngg, ""},
        {"dEpre_eta", StaggerLocation::EdgeEt, 1, ngg, ""},
        {"dEpre_zeta", StaggerLocation::EdgeZe, 1, ngg, ""},
    }};

    for (const auto &field : hall_implicit_fields)
        fld->register_field(field);
}

void MercurySolver::RegisterCouplingChannels(Field *fld, const TOPO::Topology &topology, int dimension, int ngg)
{
    const std::array<CouplingField, 21> coupling_fields = {{
        {"B_xi", StaggerLocation::FaceXi, 1},
        {"B_eta", StaggerLocation::FaceEt, 1},
        {"B_zeta", StaggerLocation::FaceZe, 1},
        {"Badd_xi", StaggerLocation::FaceXi, 1},
        {"Badd_eta", StaggerLocation::FaceEt, 1},
        {"Badd_zeta", StaggerLocation::FaceZe, 1},
        {"Eface_xi", StaggerLocation::FaceXi, 3},
        {"Eface_eta", StaggerLocation::FaceEt, 3},
        {"Eface_zeta", StaggerLocation::FaceZe, 3},
        {"E_xi", StaggerLocation::EdgeXi, 1},
        {"E_eta", StaggerLocation::EdgeEt, 1},
        {"E_zeta", StaggerLocation::EdgeZe, 1},
        {"Ehall_xi", StaggerLocation::EdgeXi, 1},
        {"Ehall_eta", StaggerLocation::EdgeEt, 1},
        {"Ehall_zeta", StaggerLocation::EdgeZe, 1},
        {"J_xi", StaggerLocation::EdgeXi, 1},
        {"J_eta", StaggerLocation::EdgeEt, 1},
        {"J_zeta", StaggerLocation::EdgeZe, 1},
        {"J_cell", StaggerLocation::Cell, 3},
        {"B_cell", StaggerLocation::Cell, 3},
        {"Bind_cell", StaggerLocation::Cell, 3},
    }};

    RegisterCouplingPair(fld, "Solid", "Fluid", coupling_fields.data(),
                         static_cast<int>(coupling_fields.size()), ngg);
    RegisterCouplingPair(fld, "Fluid", "Solid", coupling_fields.data(),
                         static_cast<int>(coupling_fields.size()), ngg);

    const std::array<CouplingField, 13> hall_implicit_coupling_fields = {{
        {"dE_xi", StaggerLocation::EdgeXi, 1},
        {"dE_eta", StaggerLocation::EdgeEt, 1},
        {"dE_zeta", StaggerLocation::EdgeZe, 1},
        {"dB_xi", StaggerLocation::FaceXi, 1},
        {"dB_eta", StaggerLocation::FaceEt, 1},
        {"dB_zeta", StaggerLocation::FaceZe, 1},
        {"dJ_xi", StaggerLocation::EdgeXi, 1},
        {"dJ_eta", StaggerLocation::EdgeEt, 1},
        {"dJ_zeta", StaggerLocation::EdgeZe, 1},
        {"dJ_cell", StaggerLocation::Cell, 3},
        {"dEpre_xi", StaggerLocation::EdgeXi, 1},
        {"dEpre_eta", StaggerLocation::EdgeEt, 1},
        {"dEpre_zeta", StaggerLocation::EdgeZe, 1},
    }};

    RegisterCouplingPair(fld, "Solid", "Fluid", hall_implicit_coupling_fields.data(),
                         static_cast<int>(hall_implicit_coupling_fields.size()), ngg);
    RegisterCouplingPair(fld, "Fluid", "Solid", hall_implicit_coupling_fields.data(),
                         static_cast<int>(hall_implicit_coupling_fields.size()), ngg);

    fld->build_coupling_buffers(topology, dimension);
}

void MercurySolver::RegisterHaloFields(Halo *halo)
{
    const std::array<const char *, 20> halo_fields = {{
        "U_H", "U_Na",
        "B_xi", "B_eta", "B_zeta",
        "Badd_xi", "Badd_eta", "Badd_zeta",
        "E_xi", "E_eta", "E_zeta",
        "Ehall_xi", "Ehall_eta", "Ehall_zeta",
        "J_xi", "J_eta", "J_zeta",
        "J_cell", "B_cell", "Bind_cell",
    }};

    for (const auto *name : halo_fields)
        halo->register_halo_field(std::string(name), HaloLevel::Vertex);

    const std::array<const char *, 13> hall_implicit_halo_fields = {{
        "dB_xi", "dB_eta", "dB_zeta",
        "dE_xi", "dE_eta", "dE_zeta",
        "dJ_xi", "dJ_eta", "dJ_zeta",
        "dEpre_xi", "dEpre_eta", "dEpre_zeta",
        "dJ_cell",
    }};

    for (const auto *name : hall_implicit_halo_fields)
        halo->register_halo_field(std::string(name), HaloLevel::Vertex);

    halo->build_registered_patterns();
}
