// SolverFields.h
#pragma once

#include <array>

// Core
#include "00_Mercury_Const.h"
#include "3_field/2_MPCNS_Field.h"
#include "4_halo/Halo_EdgeOwner_Type.h"

struct SolverFields
{
    Field *fld = nullptr;

    // ---- geometry ----
    int fid_Jac = -1;
    int fid_pinvGT_Cell = -1; // Cell ncomp=9
    IdTriplet fid_metric;     // (xi,eta,zeta) <- (JDxi,JDet,JDze)
    IdTriplet fid_pinvGT;     // (pinvGT_xi, pinvGT_eta, pinvGT_zeta)  ncomp=9
    IdTriplet fid_pinvAT;     // (pinvAT_xi, pinvAT_eta, pinvAT_zeta)  ncomp=9

    // Face metrics:
    IdTriplet Face_Area;   // Face: |S_xi|  |S_eta| |S_ze| ncomp = 1
    IdTriplet Face_dlstar; // Face: |l*_xi|  |l*_eta| |l*_ze| ncomp = 1
    IdTriplet Face_beta;   // Face: beta=  |l*_|/|S_| ncomp = 1
    // Edge metrics:
    IdTriplet Edge_metric; // Edge: S*_xi  S*_eta S*_ze ncomp = 3
    IdTriplet Edge_Astar;  // Edge: |S*_xi|  |S*_eta| |S*_ze| ncomp = 1
    IdTriplet Edge_dl;     // Edge: |e_xi|  |e_eta| |e_ze| ncomp = 1
    IdTriplet Edge_alpha;  // Edge: beta=  |l*_|/|S_| ncomp = 1
    IdTriplet Edge_dr;     // Edge: dr_xi  dr_eta dr_zeta  ncomp = 3

    // ---- field ids ----
    int fid_U_H = -1;
    int fid_U_Na = -1;
    IdTriplet fid_B;     // (xi,eta,zeta) <- (B_xi,B_eta,B_zeta)
    IdTriplet fid_E;     // (xi,eta,zeta) <- (E_xi,E_eta,E_zeta)
    IdTriplet fid_Eface; // (xi,eta,zeta) <- (E_xi,E_eta,E_zeta)
    IdTriplet fid_Ehall; // (xi,eta,zeta) <- (E_xi,E_eta,E_zeta)
    IdTriplet fid_J;     // (xi,eta,zeta) <- (J_xi,J_eta,J_zeta)

    // ---- auxiliary ----
    int fid_PV_H = -1;
    int fid_PV_Na = -1;
    int fid_Bcell = -1;
    int fid_Bindcell = -1;
    int fid_Jcell = -1;
    IdTriplet fid_Badd;
    int fid_Na = -1;
    int fid_Photo = -1;

    // ---- fluid flux and face E ----
    IdTriplet fid_F; // F_xi/F_eta/F_zeta

    // -------- buffers / time advance ----------
    int fid_RHS_H = -1;
    int fid_RHS_Na = -1;
    IdTriplet fid_RHS_b;
    int fid_U_plus = -1;
    // int fid_old_U = -1;
    int fid_divB = -1;
    // IdTriplet fid_old_Bface; // old_B_xi/eta/zeta

    // int fid_RHS_U = -1;          // RHS (cell,5)
    // IdTriplet fid_RHS_Bface;     // RHS_xi/eta/zeta (face,1)
    // IdTriplet fid_RHShall_Bface; // RHShall_xi/eta/zeta

    IdTriplet fid_dJ;
    IdTriplet fid_dE;
    IdTriplet fid_dB;
    IdTriplet fid_dEpre;
    int fid_dJcell;

    void Init(Field *fld_in)
    {
        fld = fld_in;

        if (!fld)
        {
            std::fprintf(stderr, "[SolverFields] fld is null\n");
            std::abort();
        }

        // ---- geometry ----
        fid_Jac = fld->field_id("Jac");
        fid_metric.xi = fld->field_id("JDxi");
        fid_metric.eta = fld->field_id("JDet");
        fid_metric.zeta = fld->field_id("JDze");
        fid_pinvGT.xi = fld->field_id("pinvGT_xi");
        fid_pinvGT.eta = fld->field_id("pinvGT_eta");
        fid_pinvGT.zeta = fld->field_id("pinvGT_zeta");
        fid_pinvAT.xi = fld->field_id("pinvAT_xi");
        fid_pinvAT.eta = fld->field_id("pinvAT_eta");
        fid_pinvAT.zeta = fld->field_id("pinvAT_zeta");

        fid_pinvGT_Cell = fld->field_id("pinvGT_cell");

        Face_Area.xi = fld->field_id("Area_xi");
        Face_Area.eta = fld->field_id("Area_eta");
        Face_Area.zeta = fld->field_id("Area_zeta");
        Face_dlstar.xi = fld->field_id("dlstar_xi");
        Face_dlstar.eta = fld->field_id("dlstar_eta");
        Face_dlstar.zeta = fld->field_id("dlstar_zeta");
        Face_beta.xi = fld->field_id("beta_xi");
        Face_beta.eta = fld->field_id("beta_eta");
        Face_beta.zeta = fld->field_id("beta_zeta");
        Edge_metric.xi = fld->field_id("Sstar_xi");
        Edge_metric.eta = fld->field_id("Sstar_eta");
        Edge_metric.zeta = fld->field_id("Sstar_zeta");
        Edge_Astar.xi = fld->field_id("Astar_xi");
        Edge_Astar.eta = fld->field_id("Astar_eta");
        Edge_Astar.zeta = fld->field_id("Astar_zeta");
        Edge_dl.xi = fld->field_id("dl_xi");
        Edge_dl.eta = fld->field_id("dl_eta");
        Edge_dl.zeta = fld->field_id("dl_zeta");
        Edge_alpha.xi = fld->field_id("alpha_xi");
        Edge_alpha.eta = fld->field_id("alpha_eta");
        Edge_alpha.zeta = fld->field_id("alpha_zeta");

        Edge_dr.xi = fld->field_id("dr_xi");
        Edge_dr.eta = fld->field_id("dr_eta");
        Edge_dr.zeta = fld->field_id("dr_zeta");

        // ---- field ids ----
        fid_U_H = fld->field_id("U_H");
        fid_U_Na = fld->field_id("U_Na");
        fid_B.xi = fld->field_id("B_xi");
        fid_B.eta = fld->field_id("B_eta");
        fid_B.zeta = fld->field_id("B_zeta");

        fid_E.xi = fld->field_id("E_xi");
        fid_E.eta = fld->field_id("E_eta");
        fid_E.zeta = fld->field_id("E_zeta");

        fid_Ehall.xi = fld->field_id("Ehall_xi");
        fid_Ehall.eta = fld->field_id("Ehall_eta");
        fid_Ehall.zeta = fld->field_id("Ehall_zeta");

        fid_Eface.xi = fld->field_id("Eface_xi");
        fid_Eface.eta = fld->field_id("Eface_eta");
        fid_Eface.zeta = fld->field_id("Eface_zeta");

        fid_J.xi = fld->field_id("J_xi");
        fid_J.eta = fld->field_id("J_eta");
        fid_J.zeta = fld->field_id("J_zeta");

        // ---- auxiliary ----
        fid_PV_H = fld->field_id("PV_H");
        fid_PV_Na = fld->field_id("PV_Na");
        fid_Bcell = fld->field_id("B_cell");
        fid_Jcell = fld->field_id("J_cell");
        fid_Bindcell = fld->field_id("Bind_cell");
        fid_Badd.xi = fld->field_id("Badd_xi");
        fid_Badd.eta = fld->field_id("Badd_eta");
        fid_Badd.zeta = fld->field_id("Badd_zeta");
        fid_Na = fld->field_id("Na");
        fid_Photo = fld->field_id("Photo_rate");

        fid_F.xi = fld->field_id("F_xi");
        fid_F.eta = fld->field_id("F_eta");
        fid_F.zeta = fld->field_id("F_zeta");

        // -------- buffers / time advance ----------
        fid_RHS_H = fld->field_id("RHS_H");
        fid_RHS_Na = fld->field_id("RHS_Na");
        fid_RHS_b.xi = fld->field_id("RHS_B_xi");
        fid_RHS_b.eta = fld->field_id("RHS_B_eta");
        fid_RHS_b.zeta = fld->field_id("RHS_B_zeta");
        fid_U_plus = fld->field_id("U_plus");
        // fid_old_U = fld->field_id("old_U_");
        fid_divB = fld->field_id("divB");
        // fid_old_Bface.xi = fld->field_id("old_B_xi");
        // fid_old_Bface.eta = fld->field_id("old_B_eta");
        // fid_old_Bface.zeta = fld->field_id("old_B_zeta");

        // fid_RHS_U = fld->field_id("RHS");
        // fid_RHS_Bface.xi = fld->field_id("RHS_xi");
        // fid_RHS_Bface.eta = fld->field_id("RHS_eta");
        // fid_RHS_Bface.zeta = fld->field_id("RHS_zeta");

        fid_dJ.xi = fld->field_id("dJ_xi");
        fid_dJ.eta = fld->field_id("dJ_eta");
        fid_dJ.zeta = fld->field_id("dJ_zeta");
        fid_dB.xi = fld->field_id("dB_xi");
        fid_dB.eta = fld->field_id("dB_eta");
        fid_dB.zeta = fld->field_id("dB_zeta");
        fid_dE.xi = fld->field_id("dE_xi");
        fid_dE.eta = fld->field_id("dE_eta");
        fid_dE.zeta = fld->field_id("dE_zeta");
        fid_dEpre.xi = fld->field_id("dEpre_xi");
        fid_dEpre.eta = fld->field_id("dEpre_eta");
        fid_dEpre.zeta = fld->field_id("dEpre_zeta");
        fid_dJcell = fld->field_id("dJ_cell");

        // -------- Check --------
        Validate();
    }

    void Validate() const
    {
        if (!fld)
        {
            std::fprintf(stderr, "[SolverFields] fld is null\n");
            std::abort();
        }

        auto require_id = [](int id, const char *name)
        {
            if (id < 0)
            {
                std::fprintf(stderr, "[SolverFields] missing field id: %s\n", name ? name : "(null)");
                std::abort();
            }
        };

        // ---- geometry ----
        require_id(fid_Jac, "Jac");
        require_id(fid_pinvGT_Cell, "pinvGT_Cell");
        fid_metric.require_all("metric(JDxi/JDet/JDze)");
        fid_pinvGT.require_all("pinvGT(edge)");
        fid_pinvAT.require_all("pinvAT(edge)");

        Face_Area.require_all("Face_Area");
        Face_dlstar.require_all("Face_dlstar");
        Face_beta.require_all("Face_beta");

        Edge_metric.require_all("Edge_metric");
        Edge_Astar.require_all("Edge_Astar");
        Edge_dl.require_all("Edge_dl");
        Edge_alpha.require_all("Edge_alpha");

        Edge_dr.require_all("Edge_dr");

        // ---- field ids ----
        require_id(fid_U_H, "U_H");
        require_id(fid_U_Na, "U_Na");
        fid_B.require_all("B_xi/B_eta/B_zeta");
        fid_E.require_all("E_xi/E_eta/E_zeta");
        fid_Ehall.require_all("Ehall_xi/Ehall_eta/Ehall_zeta");
        fid_Eface.require_all("Eface_xi/Eface_eta/Eface_zeta");
        fid_J.require_all("J_xi/J_eta/J_zeta");

        // ---- primary / auxiliary ----
        require_id(fid_PV_H, "PV_H");
        require_id(fid_PV_Na, "PV_Na");
        require_id(fid_Bcell, "B_cell");
        require_id(fid_Jcell, "J_cell");
        require_id(fid_Bindcell, "Bind_cell");
        fid_Badd.require_all("Badd_xi/Badd_eta/Badd_zeta");
        require_id(fid_Na, "Na");
        require_id(fid_Photo, "Photo_rate");
        fid_F.require_all("Flux(F_xi/F_eta/F_zeta)");

        // ---- buffers / time advance ----
        require_id(fid_RHS_H, "RHS_H");
        require_id(fid_RHS_Na, "RHS_Na");
        fid_RHS_b.require_all("Flux(RHS_B_xi/RHS_B_eta/RHS_B_zeta)");
        require_id(fid_U_plus, "U_plus");
        // require_id(fid_old_U, "old_U_");
        require_id(fid_divB, "divB");
        // fid_old_Bface.require_all("old_B_face(old_B_xi/eta/zeta)");

        // require_id(fid_RHS_U, "RHS");
        // fid_RHS_Bface.require_all("RHS_B_face(RHS_xi/eta/zeta)");
        fid_dJ.require_all("dJ_xi/dJ_eta/dJ_zeta");
        fid_dB.require_all("dB_xi/dB_eta/dB_zeta");
        fid_dE.require_all("dE_xi/dE_eta/dE_zeta");
        fid_dEpre.require_all("dEpre_xi/dEpre_eta/dEpre_zeta");
        require_id(fid_dJcell, "dJ_cell");
    }
};
