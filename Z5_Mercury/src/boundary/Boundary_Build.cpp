#include "1_Boundary.h"
#include "0_basic/Error.h"
#include "00_Mercury_Const.h"

void MercuryBoundary::InstallHandlers()
{
    if (!par_)
        ERROR::Abort("InstallHandlers: call Setup first");

    bool is_Mercury_interior_resis = par_->GetBoo("is_Mercury_resistance");

    auto copy = [](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
    {
        BoundaryCore::DefaultPhysicalCopy(U, fld, r, ngh);
    };

    auto nop = [](FieldBlock &, Field *, const BOUND::PhysicalRegion &, int) {};

    // -------------------------------------------------------------------------
    // 0. Default physical handlers
    //
    // 原则：
    //   Outflow / Farfield / Pole 默认 copy；
    //   Coupled-Solid / Coupled-Fluid 默认 nop。
    //
    // 后面只覆盖真正需要特殊处理的字段。
    // -------------------------------------------------------------------------
    for (const auto &fn : boundary_fields_)
    {
        RegisterPhysical_(fn, "Outflow", copy);
        RegisterPhysical_(fn, "Farfield", copy);
        RegisterPhysical_(fn, "Pole", copy);

        RegisterPhysical_(fn, "Coupled-Solid", nop);
        RegisterPhysical_(fn, "Coupled-Fluid", nop);
    }

    RegisterPhysical_("U_plus", "Outflow", nop);
    RegisterPhysical_("U_plus", "Farfield", nop);
    RegisterPhysical_("U_plus", "Pole", nop);
    RegisterPhysical_("U_plus", "Coupled-Solid", nop);
    RegisterPhysical_("U_plus", "Coupled-Fluid", nop);

    // ============================================================================================
    // Fluid BCs
    // -------------------------------------------------------------------------
    //  Fluid farfield
    // -------------------------------------------------------------------------
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

    // -------------------------------------------------------------------------
    //  Solid wall: fluid variables only
    // -------------------------------------------------------------------------
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
    // -------------------------------------------------------------------------
    //  Pole: fluid cells
    // -------------------------------------------------------------------------
    RegisterPhysical_("U_H", "Pole",
                      [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
                      {
                          this->BC_Pole_Cell_(U, fld, r, ngh);
                      });

    RegisterPhysical_("U_Na", "Pole",
                      [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
                      {
                          this->BC_Pole_Cell_(U, fld, r, ngh);
                      });
    // ============================================================================================

    // ============================================================================================
    // Electromagnetic BCs
    // -------------------------------------------------------------------------
    // Cell: B_cell J_cell
    // -------------------------------------------------------------------------
    RegisterPhysical_("B_cell", "Pole",
                      [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
                      {
                          this->BC_Pole_Bcell_Collapse_(U, fld, r, ngh);
                      });

    RegisterPhysical_("Bind_cell", "Pole",
                      [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
                      {
                          this->BC_Pole_Bcell_Collapse_(U, fld, r, ngh);
                      });
    RegisterPhysical_("B_cell", "Farfield", [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
                      { this->BC_Farfield_Bface(U, fld, r, ngh); });
    RegisterPhysical_("Bind_cell", "Farfield", [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
                      { this->BC_Farfield_Bface(U, fld, r, ngh); });
    // J_cell Pole无需特别处理，在计算的时候就已经处理过了
    // J_cell 壁面层 cell 复制域内相邻一层，虚网格再从该壁面层拷贝
    if (!is_Mercury_interior_resis)
    {
        RegisterPhysical_("J_cell", "Coupled-Solid",
                          [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
                          {
                              this->BC_Solid_Surface_Jcell(U, fld, r, ngh);
                          });
    }
    // -------------------------------------------------------------------------
    // Face: B_xi/eta/zeta Eface_xi/eta/zeta
    // -------------------------------------------------------------------------
    // B_xi/eta/zeta应该严格按照CT求解的，因此绝大部分不需要施加边界条件
    // Pole的正则性要求，需要对Pole上的进行处理
    // ------------------------------------------
    // Pole: face magnetic flux
    // B_xi / B_eta 使用 collapse。
    // B_zeta 目前显式保持 copy。不要在 Bcell handler 里顺手改 B_zeta。
    // ------------------------------------------
    auto Bface_pole = [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
    {
        this->BC_Pole_Bface_Collapse_(U, fld, r, ngh);
    };
    RegisterPhysical_("B_xi", "Pole", Bface_pole);
    RegisterPhysical_("B_eta", "Pole", Bface_pole);
    // 当前策略：B_zeta copy。
    RegisterPhysical_("B_zeta", "Pole", copy);

    RegisterPhysical_("B_xi", "Farfield", nop);
    RegisterPhysical_("B_eta", "Farfield", nop);
    RegisterPhysical_("B_zeta", "Farfield", nop);
    // ------------------------------------------
    // E_face
    // ------------------------------------------
    if (!is_Mercury_interior_resis)
    {
        auto Eface_Wall = [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
        {
            this->BC_Solid_Surface_Eface_(U, fld, r, ngh);
        };
        RegisterPhysical_("Eface_xi", "Coupled-Solid", Eface_Wall);
        RegisterPhysical_("Eface_eta", "Coupled-Solid", Eface_Wall);
        RegisterPhysical_("Eface_zeta", "Coupled-Solid", Eface_Wall);
        RegisterPhysical_("Eface_xi", "Coupled-Fluid", Eface_Wall);
        RegisterPhysical_("Eface_eta", "Coupled-Fluid", Eface_Wall);
        RegisterPhysical_("Eface_zeta", "Coupled-Fluid", Eface_Wall);
    }
    else
    {
        auto Eface_Wall = [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
        {
            this->BC_Solid_Surface_Eface_ghots_zero(U, fld, r, ngh);
        };
        RegisterPhysical_("Eface_xi", "Coupled-Solid", Eface_Wall);
        RegisterPhysical_("Eface_eta", "Coupled-Solid", Eface_Wall);
        RegisterPhysical_("Eface_zeta", "Coupled-Solid", Eface_Wall);
    }

    RegisterPhysical_("Eface_xi", "Farfield", nop);
    RegisterPhysical_("Eface_eta", "Farfield", nop);
    RegisterPhysical_("Eface_zeta", "Farfield", nop);
    // -------------------------------------------------------------------------
    // Edge: E_xi/eta/zeta
    // -------------------------------------------------------------------------
    auto Eedge_1stlayer_to_zero = [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
    {
        // no operators
        // this->BC_Solid_Surface_Eedge_(U, fld, r, ngh);
    };
    if (!is_Mercury_interior_resis)
    {
        RegisterPhysical_("E_xi", "Coupled-Solid", Eedge_1stlayer_to_zero);
        RegisterPhysical_("E_xi", "Coupled-Fluid", Eedge_1stlayer_to_zero);
        RegisterPhysical_("E_eta", "Coupled-Solid", Eedge_1stlayer_to_zero);
        RegisterPhysical_("E_eta", "Coupled-Fluid", Eedge_1stlayer_to_zero);
        RegisterPhysical_("E_zeta", "Coupled-Solid", Eedge_1stlayer_to_zero);
        RegisterPhysical_("E_zeta", "Coupled-Fluid", Eedge_1stlayer_to_zero);
    }

    auto Eedge_boundandghost_to_constant = [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
    {
        this->BC_Farfield_Eedge_set_zerocurl(U, fld, r, ngh);
    };
    RegisterPhysical_("E_xi", "Farfield", Eedge_boundandghost_to_constant);
    RegisterPhysical_("E_eta", "Farfield", Eedge_boundandghost_to_constant);
    RegisterPhysical_("E_zeta", "Farfield", Eedge_boundandghost_to_constant);
    RegisterPhysical_("Eres_xi", "Farfield", Eedge_boundandghost_to_constant);
    RegisterPhysical_("Eres_eta", "Farfield", Eedge_boundandghost_to_constant);
    RegisterPhysical_("Eres_zeta", "Farfield", Eedge_boundandghost_to_constant);
    // ------------------------------------------
    // Pole
    //   E_xi / E_eta 保留你原来的 Pole 处理；
    //   E_zeta 显式保持 copy。
    // ------------------------------------------
    auto E_xi_pole = [this, copy](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
    {
        const int dir = (r.direction < 0) ? -r.direction : r.direction;

        if (dir == 1)
            this->BC_Pole_Eedge_RegulateK_Norm(U, fld, r, ngh); // xi norm
        else if (dir == 2)
            this->BC_Pole_Eedge_Axis(U, fld, r, ngh); // eta norm, xi using Pole -u\timesB\codt dr
        else
            this->BC_Pole_Eedge_Zero_Rotate(U, fld, r, ngh);
    };

    auto E_eta_pole = [this, copy](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
    {
        const int dir = (r.direction < 0) ? -r.direction : r.direction;

        if (dir == 2)
            this->BC_Pole_Eedge_RegulateK_Norm(U, fld, r, ngh); // eta norm
        else if (dir == 1)
            this->BC_Pole_Eedge_Axis(U, fld, r, ngh); // xi norm, eta using Pole -u\timesB\codt dr
        else
            this->BC_Pole_Eedge_Zero_Rotate(U, fld, r, ngh);
    };

    auto E_zeta_pole = [this, copy](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
    {
        this->BC_Pole_Eedge_Zero_Rotate(U, fld, r, ngh);
    };

    auto Eres_xi_pole = [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
    {
        const int dir = std::abs(r.direction);

        if (dir == 1)
            this->BC_Pole_Eedge_RegulateK_Norm(U, fld, r, ngh);
        else
            this->BC_Pole_Eedge_Zero_Rotate(U, fld, r, ngh);
    };

    auto Eres_eta_pole = [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
    {
        const int dir = std::abs(r.direction);

        if (dir == 2)
            this->BC_Pole_Eedge_RegulateK_Norm(U, fld, r, ngh);
        else
            this->BC_Pole_Eedge_Zero_Rotate(U, fld, r, ngh);
    };

    auto Eres_zeta_pole = [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
    {
        this->BC_Pole_Eedge_Zero_Rotate(U, fld, r, ngh);
    };

    RegisterPhysical_("E_xi", "Pole", E_xi_pole);
    RegisterPhysical_("E_eta", "Pole", E_eta_pole);
    RegisterPhysical_("E_zeta", "Pole", E_zeta_pole);
    RegisterPhysical_("Eres_xi", "Pole", Eres_xi_pole);
    RegisterPhysical_("Eres_eta", "Pole", Eres_eta_pole);
    RegisterPhysical_("Eres_zeta", "Pole", Eres_zeta_pole);

    // ------------------------------------------
    //   J_xi / E_eta / J_zeta
    // ------------------------------------------
    // Pole
    auto J_xi_pole = [this, copy](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
    {
        const int dir = (r.direction < 0) ? -r.direction : r.direction;

        if (dir == 1)
            this->BC_Pole_Jedge_RegulateK_Norm(U, fld, r, ngh); // xi norm
        else
            this->BC_Pole_Eedge_Zero_Rotate(U, fld, r, ngh);
    };
    auto J_eta_pole = [this, copy](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
    {
        const int dir = (r.direction < 0) ? -r.direction : r.direction;

        if (dir == 2)
            this->BC_Pole_Jedge_RegulateK_Norm(U, fld, r, ngh); // eta norm
        else
            this->BC_Pole_Eedge_Zero_Rotate(U, fld, r, ngh);
    };
    auto J_zeta_pole = [this, copy](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
    {
        this->BC_Pole_Eedge_Zero_Rotate(U, fld, r, ngh);
    };
    RegisterPhysical_("J_xi", "Pole", J_xi_pole);
    RegisterPhysical_("J_eta", "Pole", J_eta_pole);
    RegisterPhysical_("J_zeta", "Pole", J_zeta_pole);

    // -------------------------------------------------------------------------
    // Coupling handlers
    //
    // 不需要显式注册 DefaultCouplingCopy。
    //
    // BoundaryCore::ApplyCouplingPair_* 中：
    //   if (h) h(...);
    //   else   DefaultCouplingCopy(...);
    //
    // 所以以前那些 RegisterCoupling_(..., ccopy) 全部是冗余的。
    // 只有以后某个 coupling channel 需要特殊处理时，才在这里单独注册。
    // -------------------------------------------------------------------------
}

// void MercuryBoundary::InstallHandlers()
// {
//     if (!par_)
//         ERROR::Abort("InstallHandlers: call Setup first");

//     auto copy = [](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
//     {
//         BoundaryCore::DefaultPhysicalCopy(U, fld, r, ngh);
//     };

//     auto nop = [](FieldBlock &, Field *, const BOUND::PhysicalRegion &, int) {};

//     // 1) 先给 boundary_fields_ 全部注册通用 handler（保证 CheckPhysicalHandlers 能过）
//     for (auto &fn : boundary_fields_)
//     {
//         RegisterPhysical_(fn, "Outflow", copy);
//         RegisterPhysical_(fn, "Pole", copy);
//         RegisterPhysical_(fn, "Farfield", copy);
//         RegisterPhysical_(fn, "Coupled-Solid", nop);
//         RegisterPhysical_(fn, "Coupled-Fluid", nop);
//     }

//     // 2) 覆盖真正需要特殊处理的：
//     RegisterPhysical_("U_H", "Farfield",
//                       [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
//                       {
//                           this->BC_UH_Farfield_H_(U, fld, r, ngh);
//                       });

//     RegisterPhysical_("U_Na", "Farfield",
//                       [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
//                       {
//                           this->BC_UH_Farfield_Na_(U, fld, r, ngh);
//                       });

//     RegisterPhysical_("U_H", "Coupled-Solid",
//                       [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
//                       {
//                           this->BC_Solid_Surface_(U, fld, r, ngh);
//                       });

//     RegisterPhysical_("U_Na", "Coupled-Solid",
//                       [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
//                       {
//                           this->BC_Solid_Surface_(U, fld, r, ngh);
//                       });

//     RegisterPhysical_("J_cell", "Coupled-Solid",
//                       [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
//                       {
//                           this->BC_Solid_Surface_Jcell(U, fld, r, ngh);
//                       });

//     auto Eface_zero_xi_ = [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
//     {
//         if (abs(r.direction) == 1)
//             BC_Solid_Surface_Eface_(U, fld, r, ngh);
//     };
//     auto Eface_zero_eta_ = [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
//     {
//         if (abs(r.direction) == 2)
//             BC_Solid_Surface_Eface_(U, fld, r, ngh);
//     };
//     auto Eface_zero_zeta_ = [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
//     {
//         if (abs(r.direction) == 3)
//             BC_Solid_Surface_Eface_(U, fld, r, ngh);
//     };

//     auto Eedge_buffer_xi_ = [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
//     {
//         // if (abs(r.direction) != 1)
//         BC_Solid_Surface_Eedge_(U, fld, r, ngh);
//     };
//     auto Eedge_buffer_eta_ = [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
//     {
//         // if (abs(r.direction) != 2)
//         BC_Solid_Surface_Eedge_(U, fld, r, ngh);
//     };
//     auto Eedge_buffer_zeta_ = [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
//     {
//         // if (abs(r.direction) != 3)
//         BC_Solid_Surface_Eedge_(U, fld, r, ngh);
//     };

//     RegisterPhysical_("Eface_xi", "Coupled-Solid", Eface_zero_xi_);
//     RegisterPhysical_("Eface_xi", "Coupled-Fluid", Eface_zero_xi_);
//     RegisterPhysical_("Eface_eta", "Coupled-Solid", Eface_zero_eta_);
//     RegisterPhysical_("Eface_eta", "Coupled-Fluid", Eface_zero_eta_);
//     RegisterPhysical_("Eface_zeta", "Coupled-Solid", Eface_zero_zeta_);
//     RegisterPhysical_("Eface_zeta", "Coupled-Fluid", Eface_zero_zeta_);

//     RegisterPhysical_("Ehall_xi", "Coupled-Solid", Eedge_buffer_xi_);
//     RegisterPhysical_("Ehall_xi", "Coupled-Fluid", Eedge_buffer_xi_);
//     RegisterPhysical_("Ehall_eta", "Coupled-Solid", Eedge_buffer_eta_);
//     RegisterPhysical_("Ehall_eta", "Coupled-Fluid", Eedge_buffer_eta_);
//     RegisterPhysical_("Ehall_zeta", "Coupled-Solid", Eedge_buffer_zeta_);
//     RegisterPhysical_("Ehall_zeta", "Coupled-Fluid", Eedge_buffer_zeta_);

//     RegisterPhysical_("E_xi", "Coupled-Solid", Eedge_buffer_xi_);
//     RegisterPhysical_("E_xi", "Coupled-Fluid", Eedge_buffer_xi_);
//     RegisterPhysical_("E_eta", "Coupled-Solid", Eedge_buffer_eta_);
//     RegisterPhysical_("E_eta", "Coupled-Fluid", Eedge_buffer_eta_);
//     RegisterPhysical_("E_zeta", "Coupled-Solid", Eedge_buffer_zeta_);
//     RegisterPhysical_("E_zeta", "Coupled-Fluid", Eedge_buffer_zeta_);

//     RegisterPhysical_("J_xi", "Coupled-Solid", Eedge_buffer_xi_);
//     RegisterPhysical_("J_xi", "Coupled-Fluid", Eedge_buffer_xi_);
//     RegisterPhysical_("J_eta", "Coupled-Solid", Eedge_buffer_eta_);
//     RegisterPhysical_("J_eta", "Coupled-Fluid", Eedge_buffer_eta_);
//     RegisterPhysical_("J_zeta", "Coupled-Solid", Eedge_buffer_zeta_);
//     RegisterPhysical_("J_zeta", "Coupled-Fluid", Eedge_buffer_zeta_);

//     //=============================================================================================
//     // Pole Boundary
//     //-------------------------------------------------------------------------
//     // Edge
//     //-------------------------------------------------------------------------
//     auto Eedge_Pole_xi_ = [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
//     {
//         const int dir = std::abs(r.direction);

//         if (dir == 1)
//             this->BC_Pole_Eedge_RegulateKAndCopyGhost_(U, fld, r, ngh);
//         else if (dir == 2)
//             this->BC_Pole_Eedge_(U, fld, r, ngh);
//         else
//             BoundaryCore::DefaultPhysicalCopy(U, fld, r, ngh);
//     };

//     auto Eedge_Pole_eta_ = [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
//     {
//         const int dir = std::abs(r.direction);

//         if (dir == 2)
//             this->BC_Pole_Eedge_RegulateKAndCopyGhost_(U, fld, r, ngh);
//         else if (dir == 1)
//             this->BC_Pole_Eedge_(U, fld, r, ngh);
//         else
//             BoundaryCore::DefaultPhysicalCopy(U, fld, r, ngh);
//     };

//     auto Eedge_Pole_xi_zero = [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
//     {
//         if (abs(r.direction) == 2) // Pole rotational direction is zeta, so norm direction should be ETA
//             BC_Pole_Eedge_Zero(U, fld, r, ngh);
//         else
//             BoundaryCore::DefaultPhysicalCopy(U, fld, r, ngh);
//     };

//     auto Eedge_Pole_eta_zero = [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
//     {
//         if (abs(r.direction) == 1) // Pole rotational direction is zeta, so norm direction should be XI
//             BC_Pole_Eedge_Zero(U, fld, r, ngh);
//         else
//             BoundaryCore::DefaultPhysicalCopy(U, fld, r, ngh);
//     };

//     RegisterPhysical_("E_xi", "Pole", Eedge_Pole_xi_);
//     RegisterPhysical_("Ehall_xi", "Pole", Eedge_Pole_xi_zero);
//     RegisterPhysical_("J_xi", "Pole", Eedge_Pole_xi_zero);
//     RegisterPhysical_("E_eta", "Pole", Eedge_Pole_eta_);
//     RegisterPhysical_("Ehall_eta", "Pole", Eedge_Pole_eta_zero);
//     RegisterPhysical_("J_eta", "Pole", Eedge_Pole_eta_zero);
//     // RegisterPhysical_("E_zeta", "Pole", Eedge_Pole_);

//     //-------------------------------------------------------------------------
//     // Face
//     //-------------------------------------------------------------------------
//     RegisterPhysical_("B_xi", "Pole", [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
//                       { this->BC_Pole_Bface_Collapse_(U, fld, r, ngh); });
//     RegisterPhysical_("B_eta", "Pole", [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
//                       { this->BC_Pole_Bface_Collapse_(U, fld, r, ngh); });

//     //-------------------------------------------------------------------------
//     // Cell
//     //-------------------------------------------------------------------------
//     RegisterPhysical_("B_cell", "Pole", [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
//                       { this->BC_Pole_Bcell_Collapse_(U, fld, r, ngh); });
//     RegisterPhysical_("Bind_cell", "Pole", [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
//                       { this->BC_Pole_Bcell_Collapse_(U, fld, r, ngh); });

//     RegisterPhysical_("J_cell", "Pole", [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
//                       { this->BC_Pole_Jcell_Collapse_(U, fld, r, ngh); });

//     RegisterPhysical_("U_H", "Pole", [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
//                       { this->BC_Pole_Cell_(U, fld, r, ngh); });
//     RegisterPhysical_("U_Na", "Pole", [this](FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh)
//                       { this->BC_Pole_Cell_(U, fld, r, ngh); });

//     //=============================================================================================
//     // 3) coupling：按你的耦合 channel 注册（先 DefaultCouplingCopy）
//     auto ccopy = [](FieldBlock &Udst, Field *fld, CouplingBufferBlock &buf,
//                     const std::string &src, const std::string &dst, const std::string &tag)
//     {
//         BoundaryCore::DefaultCouplingCopy(Udst, fld, buf, src, dst, tag);
//     };
//     auto cnooper = [](FieldBlock &Udst, Field *fld, CouplingBufferBlock &buf,
//                       const std::string &src, const std::string &dst, const std::string &tag) {};

//     // 例如 B_face 三个方向（你按实际 channel_tag/dst_field 写）
//     RegisterCoupling_("Solid", "Fluid", StaggerLocation::FaceXi, "B_xi", "B_xi", ccopy);
//     RegisterCoupling_("Solid", "Fluid", StaggerLocation::FaceEt, "B_eta", "B_eta", ccopy);
//     RegisterCoupling_("Solid", "Fluid", StaggerLocation::FaceZe, "B_zeta", "B_zeta", ccopy);
//     RegisterCoupling_("Fluid", "Solid", StaggerLocation::FaceXi, "B_xi", "B_xi", ccopy);
//     RegisterCoupling_("Fluid", "Solid", StaggerLocation::FaceEt, "B_eta", "B_eta", ccopy);
//     RegisterCoupling_("Fluid", "Solid", StaggerLocation::FaceZe, "B_zeta", "B_zeta", ccopy);

//     RegisterCoupling_("Solid", "Fluid", StaggerLocation::FaceXi, "Badd_xi", "Badd_xi", ccopy);
//     RegisterCoupling_("Solid", "Fluid", StaggerLocation::FaceEt, "Badd_eta", "Badd_eta", ccopy);
//     RegisterCoupling_("Solid", "Fluid", StaggerLocation::FaceZe, "Badd_zeta", "Badd_zeta", ccopy);
//     RegisterCoupling_("Fluid", "Solid", StaggerLocation::FaceXi, "Badd_xi", "Badd_xi", ccopy);
//     RegisterCoupling_("Fluid", "Solid", StaggerLocation::FaceEt, "Badd_eta", "Badd_eta", ccopy);
//     RegisterCoupling_("Fluid", "Solid", StaggerLocation::FaceZe, "Badd_zeta", "Badd_zeta", ccopy);

//     RegisterCoupling_("Solid", "Fluid", StaggerLocation::EdgeXi, "J_xi", "J_xi", ccopy);
//     RegisterCoupling_("Solid", "Fluid", StaggerLocation::EdgeEt, "J_eta", "J_eta", ccopy);
//     RegisterCoupling_("Solid", "Fluid", StaggerLocation::EdgeZe, "J_zeta", "J_zeta", ccopy);
//     RegisterCoupling_("Fluid", "Solid", StaggerLocation::EdgeXi, "J_xi", "J_xi", ccopy);
//     RegisterCoupling_("Fluid", "Solid", StaggerLocation::EdgeEt, "J_eta", "J_eta", ccopy);
//     RegisterCoupling_("Fluid", "Solid", StaggerLocation::EdgeZe, "J_zeta", "J_zeta", ccopy);

//     RegisterCoupling_("Solid", "Fluid", StaggerLocation::FaceXi, "Eface_xi", "Eface_xi", ccopy);
//     RegisterCoupling_("Solid", "Fluid", StaggerLocation::FaceEt, "Eface_eta", "Eface_eta", ccopy);
//     RegisterCoupling_("Solid", "Fluid", StaggerLocation::FaceZe, "Eface_zeta", "Eface_zeta", ccopy);
//     RegisterCoupling_("Fluid", "Solid", StaggerLocation::FaceXi, "Eface_xi", "Eface_xi", ccopy);
//     RegisterCoupling_("Fluid", "Solid", StaggerLocation::FaceEt, "Eface_eta", "Eface_eta", ccopy);
//     RegisterCoupling_("Fluid", "Solid", StaggerLocation::FaceZe, "Eface_zeta", "Eface_zeta", ccopy);

//     RegisterCoupling_("Solid", "Fluid", StaggerLocation::EdgeXi, "E_xi", "E_xi", ccopy);
//     RegisterCoupling_("Solid", "Fluid", StaggerLocation::EdgeEt, "E_eta", "E_eta", ccopy);
//     RegisterCoupling_("Solid", "Fluid", StaggerLocation::EdgeZe, "E_zeta", "E_zeta", ccopy);
//     RegisterCoupling_("Fluid", "Solid", StaggerLocation::EdgeXi, "E_xi", "E_xi", ccopy);
//     RegisterCoupling_("Fluid", "Solid", StaggerLocation::EdgeEt, "E_eta", "E_eta", ccopy);
//     RegisterCoupling_("Fluid", "Solid", StaggerLocation::EdgeZe, "E_zeta", "E_zeta", ccopy);

//     RegisterCoupling_("Solid", "Fluid", StaggerLocation::EdgeXi, "Ehall_xi", "Ehall_xi", ccopy);
//     RegisterCoupling_("Solid", "Fluid", StaggerLocation::EdgeEt, "Ehall_eta", "Ehall_eta", ccopy);
//     RegisterCoupling_("Solid", "Fluid", StaggerLocation::EdgeZe, "Ehall_zeta", "Ehall_zeta", ccopy);
//     RegisterCoupling_("Fluid", "Solid", StaggerLocation::EdgeXi, "Ehall_xi", "Ehall_xi", ccopy);
//     RegisterCoupling_("Fluid", "Solid", StaggerLocation::EdgeEt, "Ehall_eta", "Ehall_eta", ccopy);
//     RegisterCoupling_("Fluid", "Solid", StaggerLocation::EdgeZe, "Ehall_zeta", "Ehall_zeta", ccopy);

//     RegisterCoupling_("Solid", "Fluid", StaggerLocation::Cell, "B_cell", "B_cell", ccopy);
//     RegisterCoupling_("Fluid", "Solid", StaggerLocation::Cell, "B_cell", "B_cell", ccopy);

//     RegisterCoupling_("Solid", "Fluid", StaggerLocation::Cell, "J_cell", "J_cell", ccopy);
//     RegisterCoupling_("Fluid", "Solid", StaggerLocation::Cell, "J_cell", "J_cell", ccopy);

//     RegisterCoupling_("Solid", "Fluid", StaggerLocation::Cell, "Bind_cell", "Bind_cell", ccopy);
//     RegisterCoupling_("Fluid", "Solid", StaggerLocation::Cell, "Bind_cell", "Bind_cell", ccopy);
// }

void MercuryBoundary::InstallDefaultGroups()
{
    BoundGroup gU;
    gU.name = "Ucell";
    gU.fields = {"U_H", "U_Na"};
    gU.do_coupling = false;
    gU.do_physical = true;
    gU.do_halo = true;
    gU.halo_level = HaloLevel::Vertex;
    AddGroup(gU);

    BoundGroup gJ;
    gJ.name = "Jedge";
    gJ.fields = {"J_xi", "J_eta", "J_zeta"};
    gJ.do_coupling = true;
    gJ.do_physical = true;
    gJ.do_halo = true;
    gJ.halo_level = HaloLevel::Vertex;
    gJ.coupling_pairs = {{"Solid", "Fluid"}, {"Fluid", "Solid"}};
    AddGroup(gJ);

    BoundGroup gE;
    gE.name = "Eedge";
    gE.fields = {"E_xi", "E_eta", "E_zeta"};
    gE.do_coupling = true;
    gE.do_physical = true;
    gE.do_halo = true;
    gE.halo_level = HaloLevel::Vertex;
    gE.coupling_pairs = {{"Solid", "Fluid"}, {"Fluid", "Solid"}};
    AddGroup(gE);

    BoundGroup gEhall;
    gEhall.name = "Ehall";
    gEhall.fields = {"Ehall_xi", "Ehall_eta", "Ehall_zeta"};
    gEhall.do_coupling = true;
    gEhall.do_physical = true;
    gEhall.do_halo = true;
    gEhall.halo_level = HaloLevel::Vertex;
    gEhall.coupling_pairs = {{"Solid", "Fluid"}, {"Fluid", "Solid"}};
    AddGroup(gEhall);

    BoundGroup gEres;
    gEres.name = "Eres";
    gEres.fields = {"Eres_xi", "Eres_eta", "Eres_zeta"};
    gEres.do_coupling = true;
    gEres.do_physical = true;
    gEres.do_halo = true;
    gEres.halo_level = HaloLevel::Vertex;
    gEres.coupling_pairs = {{"Solid", "Fluid"}, {"Fluid", "Solid"}};
    AddGroup(gEres);

    BoundGroup gEres1form;
    gEres1form.name = "Eres1form";
    gEres1form.fields = {"Eres_xi", "Eres_eta", "Eres_zeta"};
    gEres1form.do_coupling = true;
    gEres1form.do_physical = true;
    gEres1form.do_halo = true;
    gEres1form.halo_level = HaloLevel::Vertex;
    gEres1form.coupling_pairs = {{"Solid", "Fluid"}, {"Fluid", "Solid"}};
    AddGroup(gEres1form);

    BoundGroup gB;
    gB.name = "Bface";
    gB.fields = {"B_xi", "B_eta", "B_zeta"};
    gB.do_coupling = true;
    gB.do_physical = true;
    gB.do_halo = true;
    gB.halo_level = HaloLevel::Vertex;
    gB.coupling_pairs = {{"Solid", "Fluid"}, {"Fluid", "Solid"}};
    AddGroup(gB);

    BoundGroup EfaceB;
    EfaceB.name = "Eface";
    EfaceB.fields = {"Eface_xi", "Eface_eta", "Eface_zeta"};
    EfaceB.do_coupling = true;
    EfaceB.do_physical = true;
    EfaceB.do_halo = true;
    EfaceB.halo_level = HaloLevel::Vertex;
    EfaceB.coupling_pairs = {{"Solid", "Fluid"}, {"Fluid", "Solid"}};
    AddGroup(EfaceB);

    BoundGroup gBc;
    gBc.name = "B_cell";
    gBc.fields = {"B_cell", "Bind_cell"};
    gBc.do_coupling = true;
    gBc.do_physical = true;
    gBc.do_halo = true;
    gBc.halo_level = HaloLevel::Vertex;
    gBc.coupling_pairs = {{"Solid", "Fluid"}, {"Fluid", "Solid"}};
    AddGroup(gBc);

    BoundGroup gJc;
    gJc.name = "J_cell";
    gJc.fields = {"J_cell"};
    gJc.do_coupling = true;
    gJc.do_physical = true;
    gJc.do_halo = true;
    gJc.halo_level = HaloLevel::Vertex;
    gJc.coupling_pairs = {{"Solid", "Fluid"}, {"Fluid", "Solid"}};
    AddGroup(gJc);

    BoundGroup gUplus;
    gUplus.name = "Uplus";
    gUplus.fields = {"U_plus"};
    gUplus.do_coupling = true;
    gUplus.do_physical = false;
    gUplus.do_halo = true;
    gUplus.halo_level = HaloLevel::Vertex;
    gUplus.coupling_pairs = {{"Solid", "Fluid"}, {"Fluid", "Solid"}};
    AddGroup(gUplus);

    BoundGroup gdE;
    gdE.name = "dE";
    gdE.fields = {"dE_xi", "dE_eta", "dE_zeta"};
    gdE.do_coupling = true;
    gdE.do_physical = true;
    gdE.do_halo = true;
    gdE.halo_level = HaloLevel::Vertex;
    gdE.coupling_pairs = {{"Solid", "Fluid"}, {"Fluid", "Solid"}};
    AddGroup(gdE);

    BoundGroup gdB;
    gdB.name = "dB";
    gdB.fields = {"dB_xi", "dB_eta", "dB_zeta"};
    gdB.do_coupling = true;
    gdB.do_physical = true;
    gdB.do_halo = true;
    gdB.halo_level = HaloLevel::Vertex;
    gdB.coupling_pairs = {{"Solid", "Fluid"}, {"Fluid", "Solid"}};
    AddGroup(gdB);

    BoundGroup gdJ;
    gdJ.name = "dJ";
    gdJ.fields = {"dJ_xi", "dJ_eta", "dJ_zeta"};
    gdJ.do_coupling = true;
    gdJ.do_physical = true;
    gdJ.do_halo = true;
    gdJ.halo_level = HaloLevel::Vertex;
    gdJ.coupling_pairs = {{"Solid", "Fluid"}, {"Fluid", "Solid"}};
    AddGroup(gdJ);

    BoundGroup gdJcell;
    gdJcell.name = "dJcell";
    gdJcell.fields = {"dJ_cell"};
    gdJcell.do_coupling = true;
    gdJcell.do_physical = true;
    gdJcell.do_halo = true;
    gdJcell.halo_level = HaloLevel::Vertex;
    gdJcell.coupling_pairs = {{"Solid", "Fluid"}, {"Fluid", "Solid"}};
    AddGroup(gdJcell);

    BoundGroup gBadd;
    gBadd.name = "Badd";
    gBadd.fields = {"Badd_xi", "Badd_eta", "Badd_zeta"};
    gBadd.do_coupling = true;
    gBadd.do_physical = true;
    gBadd.do_halo = true;
    gBadd.halo_level = HaloLevel::Vertex;
    gBadd.coupling_pairs = {{"Solid", "Fluid"}, {"Fluid", "Solid"}};
    AddGroup(gBadd);
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
    bc_state_.SetBackground(BuildMercuryBackgroundState(par_));
}
