#include "1_Boundary.h"
#include "0_basic/Error.h"
#include "00_Mercury_Const.h"

void MercuryBoundary::InstallHandlers()
{
    if (!par_)
        ERROR::Abort("InstallHandlers: call Setup first");

    InstallPhysicalHandlers_();
    InstallCouplingHandlers_();
}

void MercuryBoundary::InstallPhysicalHandlers_()
{
    InstallDefaultPhysicalHandlers_();
    InstallFarfieldPhysicalHandlers_();
    InstallCoupledPhysicalHandlers_();
    InstallPolePhysicalHandlers_();
}

void MercuryBoundary::InstallDefaultPhysicalHandlers_()
{
    auto copy = [](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
    {
        BoundaryCore::DefaultPhysicalCopy(U, fld, r, ngh);
    };

    auto nop = [](FieldBlock &, Field *, const BOUND::PhysicalRegion &, int) {};

    for (auto &fn : boundary_fields_)
    {
        RegisterPhysical_(fn, "Outflow", copy);
        RegisterPhysical_(fn, "Pole", copy);
        RegisterPhysical_(fn, "Farfield", copy);
        RegisterPhysical_(fn, "Coupled-Solid", nop);
        RegisterPhysical_(fn, "Coupled-Fluid", nop);
    }
}

void MercuryBoundary::InstallFarfieldPhysicalHandlers_()
{
    RegisterPhysical_("U_H", "Farfield",
                      [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
                      {
                          this->BC_UH_Farfield_H_(U, fld, r, ngh);
                      });

    RegisterPhysical_("U_Na", "Farfield",
                      [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
                      {
                          this->BC_UH_Farfield_Na_(U, fld, r, ngh);
                      });
}

void MercuryBoundary::InstallCoupledPhysicalHandlers_()
{
    RegisterPhysical_("U_H", "Coupled-Solid",
                      [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
                      {
                          this->BC_Solid_Surface_(U, fld, r, ngh);
                      });

    RegisterPhysical_("U_Na", "Coupled-Solid",
                      [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
                      {
                          this->BC_Solid_Surface_(U, fld, r, ngh);
                      });

    RegisterPhysical_("J_cell", "Coupled-Solid",
                      [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
                      {
                          this->BC_Solid_Surface_Jcell(U, fld, r, ngh);
                      });

    auto Eface_zero_xi_ = [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
    {
        if (abs(r.direction) == 1)
            BC_Solid_Surface_Eface_(U, fld, r, ngh);
    };
    auto Eface_zero_eta_ = [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
    {
        if (abs(r.direction) == 2)
            BC_Solid_Surface_Eface_(U, fld, r, ngh);
    };
    auto Eface_zero_zeta_ = [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
    {
        if (abs(r.direction) == 3)
            BC_Solid_Surface_Eface_(U, fld, r, ngh);
    };

    auto Eedge_buffer_xi_ = [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
    {
        // if (abs(r.direction) != 1)
        BC_Solid_Surface_Eedge_(U, fld, r, ngh);
    };
    auto Eedge_buffer_eta_ = [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
    {
        // if (abs(r.direction) != 2)
        BC_Solid_Surface_Eedge_(U, fld, r, ngh);
    };
    auto Eedge_buffer_zeta_ = [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
    {
        // if (abs(r.direction) != 3)
        BC_Solid_Surface_Eedge_(U, fld, r, ngh);
    };

    RegisterPhysical_("Eface_xi", "Coupled-Solid", Eface_zero_xi_);
    RegisterPhysical_("Eface_xi", "Coupled-Fluid", Eface_zero_xi_);
    RegisterPhysical_("Eface_eta", "Coupled-Solid", Eface_zero_eta_);
    RegisterPhysical_("Eface_eta", "Coupled-Fluid", Eface_zero_eta_);
    RegisterPhysical_("Eface_zeta", "Coupled-Solid", Eface_zero_zeta_);
    RegisterPhysical_("Eface_zeta", "Coupled-Fluid", Eface_zero_zeta_);

    RegisterPhysical_("Ehall_xi", "Coupled-Solid", Eedge_buffer_xi_);
    RegisterPhysical_("Ehall_xi", "Coupled-Fluid", Eedge_buffer_xi_);
    RegisterPhysical_("Ehall_eta", "Coupled-Solid", Eedge_buffer_eta_);
    RegisterPhysical_("Ehall_eta", "Coupled-Fluid", Eedge_buffer_eta_);
    RegisterPhysical_("Ehall_zeta", "Coupled-Solid", Eedge_buffer_zeta_);
    RegisterPhysical_("Ehall_zeta", "Coupled-Fluid", Eedge_buffer_zeta_);

    RegisterPhysical_("E_xi", "Coupled-Solid", Eedge_buffer_xi_);
    RegisterPhysical_("E_xi", "Coupled-Fluid", Eedge_buffer_xi_);
    RegisterPhysical_("E_eta", "Coupled-Solid", Eedge_buffer_eta_);
    RegisterPhysical_("E_eta", "Coupled-Fluid", Eedge_buffer_eta_);
    RegisterPhysical_("E_zeta", "Coupled-Solid", Eedge_buffer_zeta_);
    RegisterPhysical_("E_zeta", "Coupled-Fluid", Eedge_buffer_zeta_);

    RegisterPhysical_("J_xi", "Coupled-Solid", Eedge_buffer_xi_);
    RegisterPhysical_("J_xi", "Coupled-Fluid", Eedge_buffer_xi_);
    RegisterPhysical_("J_eta", "Coupled-Solid", Eedge_buffer_eta_);
    RegisterPhysical_("J_eta", "Coupled-Fluid", Eedge_buffer_eta_);
    RegisterPhysical_("J_zeta", "Coupled-Solid", Eedge_buffer_zeta_);
    RegisterPhysical_("J_zeta", "Coupled-Fluid", Eedge_buffer_zeta_);
}

void MercuryBoundary::InstallPolePhysicalHandlers_()
{
    auto Eedge_Pole_xi_ = [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
    {
        const int dir = std::abs(r.direction);

        if (dir == 1)
            this->BC_Pole_Eedge_RegulateKAndCopyGhost_(U, fld, r, ngh);
        else if (dir == 2)
            this->BC_Pole_Eedge_(U, fld, r, ngh);
        else
            BoundaryCore::DefaultPhysicalCopy(U, fld, r, ngh);
    };

    auto Eedge_Pole_eta_ = [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
    {
        const int dir = std::abs(r.direction);

        if (dir == 2)
            this->BC_Pole_Eedge_RegulateKAndCopyGhost_(U, fld, r, ngh);
        else if (dir == 1)
            this->BC_Pole_Eedge_(U, fld, r, ngh);
        else
            BoundaryCore::DefaultPhysicalCopy(U, fld, r, ngh);
    };

    auto Eedge_Pole_xi_zero = [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
    {
        if (abs(r.direction) == 2) // Pole rotational direction is zeta, so norm direction should be ETA
            BC_Pole_Eedge_Zero(U, fld, r, ngh);
        else
            BoundaryCore::DefaultPhysicalCopy(U, fld, r, ngh);
    };

    auto Eedge_Pole_eta_zero = [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
    {
        if (abs(r.direction) == 1) // Pole rotational direction is zeta, so norm direction should be XI
            BC_Pole_Eedge_Zero(U, fld, r, ngh);
        else
            BoundaryCore::DefaultPhysicalCopy(U, fld, r, ngh);
    };

    RegisterPhysical_("E_xi", "Pole", Eedge_Pole_xi_);
    RegisterPhysical_("Ehall_xi", "Pole", Eedge_Pole_xi_zero);
    RegisterPhysical_("J_xi", "Pole", Eedge_Pole_xi_zero);
    RegisterPhysical_("E_eta", "Pole", Eedge_Pole_eta_);
    RegisterPhysical_("Ehall_eta", "Pole", Eedge_Pole_eta_zero);
    RegisterPhysical_("J_eta", "Pole", Eedge_Pole_eta_zero);
    // RegisterPhysical_("E_zeta", "Pole", Eedge_Pole_);

    //-------------------------------------------------------------------------
    // Face
    //-------------------------------------------------------------------------
    RegisterPhysical_("B_xi", "Pole", [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
                      { this->BC_Pole_Bface_Collapse_(U, fld, r, ngh); });
    RegisterPhysical_("B_eta", "Pole", [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
                      { this->BC_Pole_Bface_Collapse_(U, fld, r, ngh); });

    //-------------------------------------------------------------------------
    // Cell
    //-------------------------------------------------------------------------
    RegisterPhysical_("B_cell", "Pole", [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
                      { this->BC_Pole_Bcell_Collapse_(U, fld, r, ngh); });
    RegisterPhysical_("Bind_cell", "Pole", [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
                      { this->BC_Pole_Bcell_Collapse_(U, fld, r, ngh); });

    RegisterPhysical_("J_cell", "Pole", [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
                      { this->BC_Pole_Jcell_Collapse_(U, fld, r, ngh); });

    RegisterPhysical_("U_H", "Pole", [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
                      { this->BC_Pole_Cell_(U, fld, r, ngh); });
    RegisterPhysical_("U_Na", "Pole", [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
                      { this->BC_Pole_Cell_(U, fld, r, ngh); });
}

void MercuryBoundary::InstallCouplingHandlers_()
{
    auto ccopy = [](FieldBlock &Udst, Field *fld, CouplingBufferBlock &buf,
                    const std::string &src, const std::string &dst, const std::string &tag)
    {
        BoundaryCore::DefaultCouplingCopy(Udst, fld, buf, src, dst, tag);
    };

    RegisterCoupling_("Solid", "Fluid", StaggerLocation::FaceXi, "B_xi", "B_xi", ccopy);
    RegisterCoupling_("Solid", "Fluid", StaggerLocation::FaceEt, "B_eta", "B_eta", ccopy);
    RegisterCoupling_("Solid", "Fluid", StaggerLocation::FaceZe, "B_zeta", "B_zeta", ccopy);
    RegisterCoupling_("Fluid", "Solid", StaggerLocation::FaceXi, "B_xi", "B_xi", ccopy);
    RegisterCoupling_("Fluid", "Solid", StaggerLocation::FaceEt, "B_eta", "B_eta", ccopy);
    RegisterCoupling_("Fluid", "Solid", StaggerLocation::FaceZe, "B_zeta", "B_zeta", ccopy);

    RegisterCoupling_("Solid", "Fluid", StaggerLocation::FaceXi, "Badd_xi", "Badd_xi", ccopy);
    RegisterCoupling_("Solid", "Fluid", StaggerLocation::FaceEt, "Badd_eta", "Badd_eta", ccopy);
    RegisterCoupling_("Solid", "Fluid", StaggerLocation::FaceZe, "Badd_zeta", "Badd_zeta", ccopy);
    RegisterCoupling_("Fluid", "Solid", StaggerLocation::FaceXi, "Badd_xi", "Badd_xi", ccopy);
    RegisterCoupling_("Fluid", "Solid", StaggerLocation::FaceEt, "Badd_eta", "Badd_eta", ccopy);
    RegisterCoupling_("Fluid", "Solid", StaggerLocation::FaceZe, "Badd_zeta", "Badd_zeta", ccopy);

    RegisterCoupling_("Solid", "Fluid", StaggerLocation::EdgeXi, "J_xi", "J_xi", ccopy);
    RegisterCoupling_("Solid", "Fluid", StaggerLocation::EdgeEt, "J_eta", "J_eta", ccopy);
    RegisterCoupling_("Solid", "Fluid", StaggerLocation::EdgeZe, "J_zeta", "J_zeta", ccopy);
    RegisterCoupling_("Fluid", "Solid", StaggerLocation::EdgeXi, "J_xi", "J_xi", ccopy);
    RegisterCoupling_("Fluid", "Solid", StaggerLocation::EdgeEt, "J_eta", "J_eta", ccopy);
    RegisterCoupling_("Fluid", "Solid", StaggerLocation::EdgeZe, "J_zeta", "J_zeta", ccopy);

    RegisterCoupling_("Solid", "Fluid", StaggerLocation::FaceXi, "Eface_xi", "Eface_xi", ccopy);
    RegisterCoupling_("Solid", "Fluid", StaggerLocation::FaceEt, "Eface_eta", "Eface_eta", ccopy);
    RegisterCoupling_("Solid", "Fluid", StaggerLocation::FaceZe, "Eface_zeta", "Eface_zeta", ccopy);
    RegisterCoupling_("Fluid", "Solid", StaggerLocation::FaceXi, "Eface_xi", "Eface_xi", ccopy);
    RegisterCoupling_("Fluid", "Solid", StaggerLocation::FaceEt, "Eface_eta", "Eface_eta", ccopy);
    RegisterCoupling_("Fluid", "Solid", StaggerLocation::FaceZe, "Eface_zeta", "Eface_zeta", ccopy);

    RegisterCoupling_("Solid", "Fluid", StaggerLocation::EdgeXi, "E_xi", "E_xi", ccopy);
    RegisterCoupling_("Solid", "Fluid", StaggerLocation::EdgeEt, "E_eta", "E_eta", ccopy);
    RegisterCoupling_("Solid", "Fluid", StaggerLocation::EdgeZe, "E_zeta", "E_zeta", ccopy);
    RegisterCoupling_("Fluid", "Solid", StaggerLocation::EdgeXi, "E_xi", "E_xi", ccopy);
    RegisterCoupling_("Fluid", "Solid", StaggerLocation::EdgeEt, "E_eta", "E_eta", ccopy);
    RegisterCoupling_("Fluid", "Solid", StaggerLocation::EdgeZe, "E_zeta", "E_zeta", ccopy);

    RegisterCoupling_("Solid", "Fluid", StaggerLocation::EdgeXi, "Ehall_xi", "Ehall_xi", ccopy);
    RegisterCoupling_("Solid", "Fluid", StaggerLocation::EdgeEt, "Ehall_eta", "Ehall_eta", ccopy);
    RegisterCoupling_("Solid", "Fluid", StaggerLocation::EdgeZe, "Ehall_zeta", "Ehall_zeta", ccopy);
    RegisterCoupling_("Fluid", "Solid", StaggerLocation::EdgeXi, "Ehall_xi", "Ehall_xi", ccopy);
    RegisterCoupling_("Fluid", "Solid", StaggerLocation::EdgeEt, "Ehall_eta", "Ehall_eta", ccopy);
    RegisterCoupling_("Fluid", "Solid", StaggerLocation::EdgeZe, "Ehall_zeta", "Ehall_zeta", ccopy);

    RegisterCoupling_("Solid", "Fluid", StaggerLocation::Cell, "B_cell", "B_cell", ccopy);
    RegisterCoupling_("Fluid", "Solid", StaggerLocation::Cell, "B_cell", "B_cell", ccopy);

    RegisterCoupling_("Solid", "Fluid", StaggerLocation::Cell, "J_cell", "J_cell", ccopy);
    RegisterCoupling_("Fluid", "Solid", StaggerLocation::Cell, "J_cell", "J_cell", ccopy);

    RegisterCoupling_("Solid", "Fluid", StaggerLocation::Cell, "Bind_cell", "Bind_cell", ccopy);
    RegisterCoupling_("Fluid", "Solid", StaggerLocation::Cell, "Bind_cell", "Bind_cell", ccopy);
}

void MercuryBoundary::InstallDefaultGroups()
{
    AddStandardGroup_("Ucell", {"U_H", "U_Na"}, false);
    AddStandardGroup_("Jedge", {"J_xi", "J_eta", "J_zeta"}, true);
    AddStandardGroup_("Eedge", {"E_xi", "E_eta", "E_zeta"}, true);
    AddStandardGroup_("Ehall", {"Ehall_xi", "Ehall_eta", "Ehall_zeta"}, true);
    AddStandardGroup_("Bface", {"B_xi", "B_eta", "B_zeta"}, true);
    AddStandardGroup_("Eface", {"Eface_xi", "Eface_eta", "Eface_zeta"}, true);
    AddStandardGroup_("B_cell", {"B_cell", "Bind_cell"}, true);
    AddStandardGroup_("J_cell", {"J_cell"}, true);

    AddStandardGroup_("dE", {"dE_xi", "dE_eta", "dE_zeta"}, true);
    AddStandardGroup_("dB", {"dB_xi", "dB_eta", "dB_zeta"}, true);
    AddStandardGroup_("dJ", {"dJ_xi", "dJ_eta", "dJ_zeta"}, true);
    AddStandardGroup_("dJcell", {"dJ_cell"}, true);
    AddStandardGroup_("dEpre", {"dEpre_xi", "dEpre_eta", "dEpre_zeta"}, true);

    AddStandardGroup_("Badd", {"Badd_xi", "Badd_eta", "Badd_zeta"}, true);
}

void MercuryBoundary::Build(bool strict_check)
{
    if (!halo_)
        ERROR::Abort("Build: call Setup first");

    // build halo patterns once
    // halo_->build_registered_patterns();

    if (strict_check)
        bound_.CheckPhysicalHandlers(boundary_fields_);

    built_ = true;
}

void MercuryBoundary::InitBCStateFromParam_()
{
    // ---- Constants ----
    bc_state_.gamma = par_->GetDou_List("constant").data["gamma"];
    double NA = par_->GetDou_List("constant").data["NA"];
    double R_uni = par_->GetDou_List("constant").data["R_uni"];
    double q_e = par_->GetDou_List("constant").data["q_e"];
    double k_Boltz = R_uni / NA;
    double mu_mag = par_->GetDou_List("constant").data["mu_mag"];

    List<double> ini = par_->GetDou_List("INITIAL");
    List<double> ref = par_->GetDou_List("REF");

    // ---- IMF（Tesla）----
    double Bx_phy = ini.data["Bx"];
    double By_phy = ini.data["By"];
    double Bz_phy = ini.data["Bz"];

    double B_ref = ref.data["B_ref"];

    // Nondimensional IMF
    double Bx = Bx_phy / B_ref;
    double By = By_phy / B_ref;
    double Bz = Bz_phy / B_ref;

    // Fluid reference physical quantities
    double L_ref = ref.data["L_ref"];
    double U_ref = ref.data["U"];
    double n_ref = ref.data["n"];
    double T_ref = ref.data["T"];
    double Molecular_mass_ref = ref.data["Molecular_mass"];
    double rho_ref = Molecular_mass_ref / NA * n_ref;

    // ---- nondimensional primitive vars ----
    double rho0 = (ini.data["n"] / n_ref) * (ini.data["Molecular_mass"] / Molecular_mass_ref);

    double c_y = ini.data["c_y"];
    double c_z = ini.data["c_z"];
    double c_x = -std::sqrt(1.0 - c_y * c_y - c_z * c_z);
    double u0 = c_x * ini.data["U"] / U_ref;
    double v0 = c_y * ini.data["U"] / U_ref;
    double w0 = c_z * ini.data["U"] / U_ref;

    double p_ini = ini.data["n"] * k_Boltz * ini.data["T"];
    double p0 = p_ini / (rho_ref * U_ref * U_ref);

    double T0 = ini.data["T"] / T_ref;

    // Solar Wind state
    bc_state_.q_pv_inf[0] = u0;
    bc_state_.q_pv_inf[1] = v0;
    bc_state_.q_pv_inf[2] = w0;
    bc_state_.q_pv_inf[3] = p0;
    bc_state_.q_pv_inf[4] = T0;

    bc_state_.qinf[0] = rho0;
    bc_state_.qinf[1] = rho0 * u0;
    bc_state_.qinf[2] = rho0 * v0;
    bc_state_.qinf[3] = rho0 * w0;
    bc_state_.qinf[4] = 0.5 * rho0 * (u0 * u0 + v0 * v0 + w0 * w0) // Kinetic energy
                        + p0 / (bc_state_.gamma - 1.0);            // Inertial energy

    // IMF
    bc_state_.B_imf[0] = Bx;
    bc_state_.B_imf[1] = By;
    bc_state_.B_imf[2] = Bz;

    // qinfs[5], q_pv_infs[5];// Na+: seed initial state

    bc_state_.q_pv_infs[0] = 0.0;                   // Na seeds are assumed to be static
    bc_state_.q_pv_infs[1] = 0.0;                   // Na seeds are assumed to be static
    bc_state_.q_pv_infs[2] = 0.0;                   // Na seeds are assumed to be static
    bc_state_.q_pv_infs[3] = p0 * rho_small / 23.0; // Very low background pressure
    bc_state_.q_pv_infs[4] = T0;                    // Temperature is the same as inflow

    bc_state_.qinfs[0] = rho_small * rho0;
    bc_state_.qinfs[1] = 0.0;
    bc_state_.qinfs[2] = 0.0;
    bc_state_.qinfs[3] = 0.0;
    bc_state_.qinfs[4] = bc_state_.q_pv_infs[3] / (bc_state_.gamma - 1.0) + 0.5 * bc_state_.qinfs[0] * (bc_state_.q_pv_infs[0] * bc_state_.q_pv_infs[0] + bc_state_.q_pv_infs[1] * bc_state_.q_pv_infs[1] + bc_state_.q_pv_infs[2] * bc_state_.q_pv_infs[2]);
}
