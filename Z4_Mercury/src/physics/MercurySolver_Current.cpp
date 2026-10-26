#include "MercurySolver.h"

#include "mpcns/discrete/lumped_hodge.hpp"
#include "mpcns/operators/electromagnetic_operators.hpp"
#include "mpcns/operators/pole_current_constraint.hpp"
#include "mpcns/operators/reconstruction_operators.hpp"

void MercurySolver::ReduceEdgeAliasCandidatesToOwners_(const IdTriplet &fid_edge)
{
    if(!topo_) return;
    const int myid=par_->GetInt("myid");
    if(!edge_alias_reduce_cache_ready_)
    {
        edge_alias_reduce_classes_.clear();
        edge_alias_reduce_classes_.reserve(topo_->edges.classes.size());
        for(const auto &cls:topo_->edges.classes)
            if(cls.global_id>=0 && cls.members.size()>1)
                edge_alias_reduce_classes_.push_back(&cls);
        std::sort(edge_alias_reduce_classes_.begin(),edge_alias_reduce_classes_.end(),
                  [](const auto *a,const auto *b)
                  { return a->global_id<b->global_id; });
        edge_alias_reduce_local_sum_.resize(edge_alias_reduce_classes_.size());
        edge_alias_reduce_global_sum_.resize(edge_alias_reduce_classes_.size());

        edge_alias_reduce_local_terms_.clear();
        edge_alias_reduce_owner_writes_.clear();
        for(std::size_t n=0;n<edge_alias_reduce_classes_.size();++n)
        {
            const auto &cls=*edge_alias_reduce_classes_[n];
            for(const auto &m:cls.members)
                if(m.entity.rank==myid)
                    edge_alias_reduce_local_terms_.push_back(
                        {n,m.entity,m.orient_sign});
            if(cls.owner.entity.rank==myid)
                edge_alias_reduce_owner_writes_.push_back(
                    {n,cls.owner.entity,cls.owner.orient_sign,
                     static_cast<double>(cls.members.size())});
        }
        edge_alias_reduce_cache_ready_=true;
    }

    const auto &classes=edge_alias_reduce_classes_;
    const std::size_t nclass=classes.size();
    if(nclass==0) return;
    auto &local_sum=edge_alias_reduce_local_sum_;
    auto &global_sum=edge_alias_reduce_global_sum_;
    std::fill(local_sum.begin(),local_sum.end(),0.0);
    for(const auto &term:edge_alias_reduce_local_terms_)
    {
        const auto &e=term.edge;
        FieldBlock &value=fld_->field(fid_edge.at(static_cast<int>(e.axis)+1),e.block);
        local_sum[term.class_index]+=term.orient_sign*value(e.i,e.j,e.k,0);
    }
    MPI_Allreduce(local_sum.data(),global_sum.data(),static_cast<int>(nclass),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    for(const auto &write:edge_alias_reduce_owner_writes_)
    {
        const auto &e=write.edge;
        FieldBlock &owner=fld_->field(fid_edge.at(static_cast<int>(e.axis)+1),e.block);
        owner(e.i,e.j,e.k,0)=write.orient_sign*
            global_sum[write.class_index]/write.member_count;
    }
}

void MercurySolver::AssembleSingularEdgeCurrent_(const IdTriplet &fid_Bface,
                                                 const IdTriplet &fid_Jedge)
{
    if (!singular_edges_ || singular_edges_->empty()) return;

    auto incidence_sign = [](const TOPO::EntityKey &edge,
                             const TOPO::EntityKey &face) -> int
    {
        using TOPO::EntityAxis;
        if (edge.rank != face.rank || edge.block != face.block)
            return 0;

        if (edge.axis == EntityAxis::Xi)
        {
            if (face.axis == EntityAxis::Eta && face.i == edge.i && face.j == edge.j)
            {
                if (face.k == edge.k - 1) return +1;
                if (face.k == edge.k) return -1;
            }
            if (face.axis == EntityAxis::Zeta && face.i == edge.i && face.k == edge.k)
            {
                if (face.j == edge.j) return +1;
                if (face.j == edge.j - 1) return -1;
            }
        }
        else if (edge.axis == EntityAxis::Eta)
        {
            if (face.axis == EntityAxis::Xi && face.i == edge.i && face.j == edge.j)
            {
                if (face.k == edge.k) return +1;
                if (face.k == edge.k - 1) return -1;
            }
            if (face.axis == EntityAxis::Zeta && face.j == edge.j && face.k == edge.k)
            {
                if (face.i == edge.i - 1) return +1;
                if (face.i == edge.i) return -1;
            }
        }
        else if (edge.axis == EntityAxis::Zeta)
        {
            if (face.axis == EntityAxis::Xi && face.i == edge.i && face.k == edge.k)
            {
                if (face.j == edge.j - 1) return +1;
                if (face.j == edge.j) return -1;
            }
            if (face.axis == EntityAxis::Eta && face.j == edge.j && face.k == edge.k)
            {
                if (face.i == edge.i) return +1;
                if (face.i == edge.i - 1) return -1;
            }
        }
        return 0;
    };

    auto contribution = [&](const METRIC::SingularPhysicalEdge &edge,
                            const METRIC::WeightedIncidentEntity &inc) -> double
    {
        const auto &face = inc.entity;
        const int sign = incidence_sign(inc.source_alias, face);
        if (sign == 0) return 0.0;

        const int dir = static_cast<int>(face.axis) + 1;
        FieldBlock &B = fld_->field(fid_Bface.at(dir), face.block);
        FieldBlock &beta = fld_->field(fid_.Face_beta.at(dir), face.block);
        if (!B.is_allocated() || !beta.is_allocated()) return 0.0;

        // source_orientation maps the local positive edge direction to the
        // canonical physical-edge direction.  local_incident_faces contains
        // each quotient face exactly once, so this is a sum, not an average.
        return edge.inverse_hodge * inc.source_orientation * sign *
               beta(face.i, face.j, face.k, 0) * B(face.i, face.j, face.k, 0);
    };

    singular_edges_->assemble_face_triplet_to_local_owners(
        *fld_,
        {fld_->descriptor(fid_Jedge.xi).name,
         fld_->descriptor(fid_Jedge.eta).name,
         fld_->descriptor(fid_Jedge.zeta).name},
        contribution);
}

void MercurySolver::ComputeRegularEdgeCurrent_(
    const IdTriplet &fid_Bface,
    const IdTriplet &fid_Jedge)
{
    for (int iblk = 0; iblk < fld_->num_blocks(); ++iblk)
    {
        FieldBlock &Bxi = fld_->field(fid_Bface.xi, iblk);
        FieldBlock &Beta = fld_->field(fid_Bface.eta, iblk);
        FieldBlock &Bzeta = fld_->field(fid_Bface.zeta, iblk);

        FieldBlock &Jxi = fld_->field(fid_Jedge.xi, iblk);
        FieldBlock &Jeta = fld_->field(fid_Jedge.eta, iblk);
        FieldBlock &Jzeta = fld_->field(fid_Jedge.zeta, iblk);

        FieldBlock &star_xi = fld_->field(fid_.fid_Jstar.xi, iblk);
        FieldBlock &star_eta = fld_->field(fid_.fid_Jstar.eta, iblk);
        FieldBlock &star_zeta = fld_->field(fid_.fid_Jstar.zeta, iblk);

        if (!Bxi.is_allocated() || !Beta.is_allocated() || !Bzeta.is_allocated() ||
            !Jxi.is_allocated() || !Jeta.is_allocated() || !Jzeta.is_allocated() ||
            !star_xi.is_allocated() || !star_eta.is_allocated() || !star_zeta.is_allocated())
            continue;

        const mpcns::discrete::SnapshotLumpedHodge hodge{
            geometry_snapshot_.block(static_cast<std::size_t>(iblk)).lumped_hodge};
        const mpcns::operators::ElectromagneticOperators electromagnetic{hodge};

        const mpcns::state::PrimalFaceFormView magnetic_flux{
            Bxi, Beta, Bzeta};
        const mpcns::operators::MagneticCurrentWorkspace workspace{
            mpcns::state::DualEdgeFormView{star_xi, star_eta, star_zeta},
            mpcns::state::DualFaceFormView{Jxi, Jeta, Jzeta}};
        const mpcns::state::PrimalEdgeFormView current{
            Jxi, Jeta, Jzeta};

        electromagnetic.compute_current(magnetic_flux, workspace, current);
    }
}

void MercurySolver::Calc_J_Edge()
{
    // B_xi/B_eta/B_zeta are induced face-flux DOFs. Badd is analytical and
    // deliberately stays outside the solved-current operator.
    ComputeRegularEdgeCurrent_(fid_.fid_B, fid_.fid_J);

    // These three operations are intentionally outside the local operator:
    // quotient aliases require a global reduction, extraordinary edges use
    // variable-valence geometry, and the boundary owns the final halo policy.
    ReduceEdgeAliasCandidatesToOwners_(fid_.fid_J);
    AssembleSingularEdgeCurrent_(fid_.fid_B, fid_.fid_J);
    mercury_bound_.Sync("Jedge");
}

void MercurySolver::calc_Jcell()
{
    const int nblock = fld_->num_blocks();

    for (int ib = 0; ib < nblock; ++ib)
    {
        auto &Jcell = fld_->field(fid_.fid_Jcell, ib);
        auto &Jxi = fld_->field(fid_.fid_J.xi, ib);
        auto &Jeta = fld_->field(fid_.fid_J.eta, ib);
        auto &Jzeta = fld_->field(fid_.fid_J.zeta, ib);
        auto &weights = fld_->field(fid_.fid_Jcell_from_Jedge_w, ib);

        const mpcns::state::PrimalEdgeFormView current{Jxi, Jeta, Jzeta};
        mpcns::operators::WeightedCellReconstruction::edge_1form_to_cell_vector(
            current, weights, Jcell);
    }

    if (topo_)
    {
        std::vector<mpcns::operators::PoleCurrentCellPatch> pole_patches;
        pole_patches.reserve(topo_->physical_patches.size());
        for (const auto &p : topo_->physical_patches)
        {
            if (p.bc_name != "Pole" ||
                (std::abs(p.direction) != 1 && std::abs(p.direction) != 2) ||
                p.this_block < 0 || p.this_block >= nblock)
                continue;

            const int ib = p.this_block;
            pole_patches.push_back({
                p.direction,
                p.this_box_node.lo.i, p.this_box_node.lo.j, p.this_box_node.lo.k,
                p.this_box_node.hi.i, p.this_box_node.hi.j, p.this_box_node.hi.k,
                &fld_->field(fid_.fid_Jcell, ib),
                &fld_->field(fid_.fid_J.xi, ib),
                &fld_->field(fid_.fid_J.eta, ib),
                &fld_->field(fid_.fid_J.zeta, ib),
                &fld_->field(fid_.Edge_dr.xi, ib),
                &fld_->field(fid_.Edge_dr.eta, ib),
                &fld_->field(fid_.Edge_dr.zeta, ib)});
        }
        mpcns::operators::PoleCurrentCellConstraint::apply(pole_patches);
    }

    mercury_bound_.Sync("J_cell");

}
