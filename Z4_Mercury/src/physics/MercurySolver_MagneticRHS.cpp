#include "MercurySolver.h"
#include "0_basic/LayoutTraits.h"
#include <algorithm>
#include <vector>

void MercurySolver::ApplyStationaryWallIdealEMF_()
{
    if (!topo_)
        return;

    const struct EdgeField
    {
        int fid;
        StaggerLocation loc;
        int axis;
    } edge_fields[3] = {
        {fid_.fid_E.xi, StaggerLocation::EdgeXi, 0},
        {fid_.fid_E.eta, StaggerLocation::EdgeEt, 1},
        {fid_.fid_E.zeta, StaggerLocation::EdgeZe, 2}};

    auto zero_wall_face = [&](int block, const Box3 &node_box, int direction)
    {
        if (block < 0 || block >= fld_->num_blocks())
            return;
        const int normal_axis = std::abs(direction) - 1;
        if (normal_axis < 0 || normal_axis > 2)
            return;

        for (const auto &ef : edge_fields)
        {
            // u_wall=0 => ideal E=-u x B=0.  Only edge DOFs tangent
            // to the exact wall face belong to this boundary condition.
            if (ef.axis == normal_axis)
                continue;

            FieldBlock &E = fld_->field(ef.fid, block);
            if (!E.is_allocated())
                continue;

            Box3 wall_edges = LAYOUT::node_box_to_dof_box(ef.loc, node_box);
            const Int3 lo = E.inner_lo();
            const Int3 hi = E.inner_hi();
            wall_edges.lo.i = std::max(wall_edges.lo.i, lo.i);
            wall_edges.lo.j = std::max(wall_edges.lo.j, lo.j);
            wall_edges.lo.k = std::max(wall_edges.lo.k, lo.k);
            wall_edges.hi.i = std::min(wall_edges.hi.i, hi.i);
            wall_edges.hi.j = std::min(wall_edges.hi.j, hi.j);
            wall_edges.hi.k = std::min(wall_edges.hi.k, hi.k);

            for (int i = wall_edges.lo.i; i < wall_edges.hi.i; ++i)
                for (int j = wall_edges.lo.j; j < wall_edges.hi.j; ++j)
                    for (int k = wall_edges.lo.k; k < wall_edges.hi.k; ++k)
                        E(i, j, k, 0) = 0.0;
        }
    };

    for (const auto &patch : topo_->physical_patches)
    {
        if (patch.bc_name == "Solid_Surface" ||
            patch.bc_name == "Coupled-Solid" ||
            patch.bc_name == "Coupled-Fluid")
            zero_wall_face(patch.this_block, patch.this_box_node, patch.direction);
    }

    // Coupling-derived physical patches represent only the `this` side of an
    // InterfacePatch.  For an inner Solid/Fluid interface both block-local
    // aliases live on this rank, so impose the same exact-wall condition on
    // the neighbor side as well.  Parallel interfaces have a reciprocal
    // local patch on the neighbor rank and need only their local `this` side.
    for (const auto &patch : topo_->inner_patches)
    {
        if (!patch.is_coupling ||
            (patch.this_block_name != "Solid" && patch.nb_block_name != "Solid"))
            continue;
        zero_wall_face(patch.this_block, patch.this_box_node, patch.direction);
        zero_wall_face(patch.nb_block, patch.nb_box_node, patch.nb_direction);
    }
    for (const auto &patch : topo_->parallel_patches)
    {
        if (!patch.is_coupling ||
            (patch.this_block_name != "Solid" && patch.nb_block_name != "Solid"))
            continue;
        zero_wall_face(patch.this_block, patch.this_box_node, patch.direction);
    }
}

void MercurySolver::AssembleRHS_Induction_CT_()
{
    const int nb = fld_->num_blocks();

    // 0) E_edge 清零
    for (int ib = 0; ib < nb; ++ib)
    {
        auto &Exi = fld_->field(fid_.fid_E.xi, ib);
        auto &Eeta = fld_->field(fid_.fid_E.eta, ib);
        auto &Eze = fld_->field(fid_.fid_E.zeta, ib);
        if (!Exi.is_allocated())
            continue;

        auto zero_electric = [&](FieldBlock &E)
        {
            Int3 lo = E.get_lo(), hi = E.get_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        E(i, j, k, 0) = 0.0;
        };
        zero_electric(Exi);
        zero_electric(Eeta);
        zero_electric(Eze);
    }

    AddIdealEdgeEMF_();
    ApplyStationaryWallIdealEMF_();
    AddAmbipolarEdgeEMF_();
    AddArtificialResistivityToEdgeEMF_();
    AddLocalArtificialResistivityToEdgeEMF_();
    AssembleSingularEdgeEMF_NonHall_();

#if HALL_EXPLICIT == 1
    AddHallEdgeEMF_();
    mercury_bound_.Sync("Ehall");
    for (int ib=0; ib<nb; ++ib)
    {
        auto add_one=[](FieldBlock &e,FieldBlock &h)
        {
            if (!e.is_allocated() || !h.is_allocated()) return;
            const Int3 lo=e.inner_lo(), hi=e.inner_hi();
            for(int i=lo.i;i<hi.i;++i) for(int j=lo.j;j<hi.j;++j) for(int k=lo.k;k<hi.k;++k)
                e(i,j,k,0)+=h(i,j,k,0);
        };
        add_one(fld_->field(fid_.fid_E.xi,ib),fld_->field(fid_.fid_Ehall.xi,ib));
        add_one(fld_->field(fid_.fid_E.eta,ib),fld_->field(fid_.fid_Ehall.eta,ib));
        add_one(fld_->field(fid_.fid_E.zeta,ib),fld_->field(fid_.fid_Ehall.zeta,ib));
    }
#endif

    mercury_bound_.Sync("Eedge");

    for (int ib = 0; ib < nb; ++ib)
    {
        auto &Exi = fld_->field(fid_.fid_E.xi, ib);
        auto &Eeta = fld_->field(fid_.fid_E.eta, ib);
        auto &Eze = fld_->field(fid_.fid_E.zeta, ib);

        auto &RHSBxi = fld_->field(fid_.fid_RHS_b.xi, ib);
        auto &RHSBeta = fld_->field(fid_.fid_RHS_b.eta, ib);
        auto &RHSBze = fld_->field(fid_.fid_RHS_b.zeta, ib);
        if (!Exi.is_allocated())
            continue;

        CTOperators::CurlEdgeToFace(ib, Exi, Eeta, Eze, RHSBxi, RHSBeta, RHSBze, /*multiper=*/-1.0);
    }
}

void MercurySolver::ResistiveDiffusionSubcycles_()
{
    if (!resist_control.is_Mercury_resistance)
        return;

    const int nsub = std::max(1, resist_control.n_subcycles);
    const double dt_res = dt / static_cast<double>(nsub);
    const int nb = fld_->num_blocks();

    auto zero_one = [](FieldBlock &F)
    {
        if (!F.is_allocated())
            return;

        Int3 lo = F.get_lo();
        Int3 hi = F.get_hi();
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                    for (int m = 0; m < F.descriptor().ncomp; ++m)
                        F(i, j, k, m) = 0.0;
    };

    for (int isub = 0; isub < nsub; ++isub)
    {
        mercury_bound_.Sync("Bface");
        UpdateMagneticDerivedFields_();

        for (int ib = 0; ib < nb; ++ib)
        {
            zero_one(fld_->field(fid_.fid_Eres.xi, ib));
            zero_one(fld_->field(fid_.fid_Eres.eta, ib));
            zero_one(fld_->field(fid_.fid_Eres.zeta, ib));
            zero_one(fld_->field(fid_.fid_RHS_b_res.xi, ib));
            zero_one(fld_->field(fid_.fid_RHS_b_res.eta, ib));
            zero_one(fld_->field(fid_.fid_RHS_b_res.zeta, ib));
        }

        AddResistiveEdgeEMF_To_(fid_.fid_Eres);
        mercury_bound_.Sync("Eres");

        for (int ib = 0; ib < nb; ++ib)
        {
            auto &Exi = fld_->field(fid_.fid_Eres.xi, ib);
            auto &Eeta = fld_->field(fid_.fid_Eres.eta, ib);
            auto &Eze = fld_->field(fid_.fid_Eres.zeta, ib);

            auto &RHSBxi = fld_->field(fid_.fid_RHS_b_res.xi, ib);
            auto &RHSBeta = fld_->field(fid_.fid_RHS_b_res.eta, ib);
            auto &RHSBze = fld_->field(fid_.fid_RHS_b_res.zeta, ib);

            if (!Exi.is_allocated())
                continue;

            CTOperators::CurlEdgeToFace(ib, Exi, Eeta, Eze,
                                        RHSBxi, RHSBeta, RHSBze,
                                        /*multiper=*/-1.0);
        }

        ApplyUpdate_Euler_BfaceOnly_(dt_res, fid_.fid_RHS_b_res);
        mercury_bound_.Sync("Bface");
    }
}
