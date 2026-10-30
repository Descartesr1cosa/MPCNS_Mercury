#include "MercurySolver.h"
#include <set>
#include <map>

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

        edge_alias_reduce_local_terms_.clear();
        edge_alias_reduce_owner_writes_.clear();
        edge_alias_reduce_send_peers_.clear();
        edge_alias_reduce_recv_peers_.clear();
        std::map<int,EdgeAliasReducePeer_> send_peers,recv_peers;
        std::size_t nslot=0;
        for(const auto *clsp:edge_alias_reduce_classes_)
        {
            const auto &cls=*clsp;
            const int owner_rank=cls.owner.entity.rank;
            bool mine=false;
            std::set<int> member_ranks;
            for(const auto &m:cls.members)
            {
                member_ranks.insert(m.entity.rank);
                if(m.entity.rank==myid)
                {
                    mine=true;
                    edge_alias_reduce_local_terms_.push_back(
                        {nslot,m.entity,m.orient_sign});
                }
            }
            if(!mine) continue;
            const std::size_t slot=nslot++;
            if(owner_rank!=myid)
            {
                auto &peer=send_peers[owner_rank];
                peer.rank=owner_rank;
                peer.slots.push_back(slot);
                continue;
            }
            const std::size_t w=edge_alias_reduce_owner_writes_.size();
            edge_alias_reduce_owner_writes_.push_back(
                {slot,cls.owner.entity,cls.owner.orient_sign,
                 static_cast<double>(cls.members.size())});
            for(int r:member_ranks)
                if(r!=myid)
                {
                    auto &peer=recv_peers[r];
                    peer.rank=r;
                    peer.writes.push_back(w);
                }
        }
        for(auto &kv:send_peers)
        {
            kv.second.buffer.resize(kv.second.slots.size());
            edge_alias_reduce_send_peers_.push_back(std::move(kv.second));
        }
        for(auto &kv:recv_peers)
        {
            kv.second.buffer.resize(kv.second.writes.size());
            edge_alias_reduce_recv_peers_.push_back(std::move(kv.second));
        }
        edge_alias_reduce_local_sum_.assign(nslot,0.0);
        edge_alias_reduce_owner_sum_.assign(edge_alias_reduce_owner_writes_.size(),0.0);
        edge_alias_reduce_cache_ready_=true;
    }

    auto &local_sum=edge_alias_reduce_local_sum_;
    auto &owner_sum=edge_alias_reduce_owner_sum_;
    std::fill(local_sum.begin(),local_sum.end(),0.0);
    for(const auto &term:edge_alias_reduce_local_terms_)
    {
        const auto &e=term.edge;
        FieldBlock &value=fld_->field(fid_edge.at(static_cast<int>(e.axis)+1),e.block);
        local_sum[term.slot]+=term.orient_sign*value(e.i,e.j,e.k,0);
    }

    constexpr int tag=7351;
    std::vector<MPI_Request> requests;
    requests.reserve(edge_alias_reduce_send_peers_.size()+edge_alias_reduce_recv_peers_.size());
    for(auto &peer:edge_alias_reduce_recv_peers_)
    {
        requests.emplace_back();
        MPI_Irecv(peer.buffer.data(),static_cast<int>(peer.buffer.size()),MPI_DOUBLE,
                  peer.rank,tag,MPI_COMM_WORLD,&requests.back());
    }
    for(auto &peer:edge_alias_reduce_send_peers_)
    {
        for(std::size_t n=0;n<peer.slots.size();++n)
            peer.buffer[n]=local_sum[peer.slots[n]];
        requests.emplace_back();
        MPI_Isend(peer.buffer.data(),static_cast<int>(peer.buffer.size()),MPI_DOUBLE,
                  peer.rank,tag,MPI_COMM_WORLD,&requests.back());
    }
    for(std::size_t w=0;w<edge_alias_reduce_owner_writes_.size();++w)
        owner_sum[w]=local_sum[edge_alias_reduce_owner_writes_[w].slot];
    if(!requests.empty())
        MPI_Waitall(static_cast<int>(requests.size()),requests.data(),MPI_STATUSES_IGNORE);
    for(const auto &peer:edge_alias_reduce_recv_peers_)
        for(std::size_t n=0;n<peer.writes.size();++n)
            owner_sum[peer.writes[n]]+=peer.buffer[n];

    for(std::size_t w=0;w<edge_alias_reduce_owner_writes_.size();++w)
    {
        const auto &write=edge_alias_reduce_owner_writes_[w];
        const auto &e=write.edge;
        FieldBlock &owner=fld_->field(fid_edge.at(static_cast<int>(e.axis)+1),e.block);
        owner(e.i,e.j,e.k,0)=write.orient_sign*owner_sum[w]/write.member_count;
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
    const auto state_b_xi = state_layout_binding_.quantity_id_for_legacy_field_id(fid_Bface.xi);
    const auto state_b_eta = state_layout_binding_.quantity_id_for_legacy_field_id(fid_Bface.eta);
    const auto state_b_zeta = state_layout_binding_.quantity_id_for_legacy_field_id(fid_Bface.zeta);
    const auto state_j_xi = state_layout_binding_.quantity_id_for_legacy_field_id(fid_Jedge.xi);
    const auto state_j_eta = state_layout_binding_.quantity_id_for_legacy_field_id(fid_Jedge.eta);
    const auto state_j_zeta = state_layout_binding_.quantity_id_for_legacy_field_id(fid_Jedge.zeta);

    for (int iblk = 0; iblk < fld_->num_blocks(); ++iblk)
    {
        FieldBlock &Bxi = fld_->field(fid_Bface.xi, iblk);
        FieldBlock &Beta = fld_->field(fid_Bface.eta, iblk);
        FieldBlock &Bzeta = fld_->field(fid_Bface.zeta, iblk);

        FieldBlock &Jxi = fld_->field(fid_Jedge.xi, iblk);
        FieldBlock &Jeta = fld_->field(fid_Jedge.eta, iblk);
        FieldBlock &Jzeta = fld_->field(fid_Jedge.zeta, iblk);

        if (!Bxi.is_allocated() || !Beta.is_allocated() || !Bzeta.is_allocated() ||
            !Jxi.is_allocated() || !Jeta.is_allocated() || !Jzeta.is_allocated())
            continue;
        if (static_cast<std::size_t>(iblk) >= jstar_state_views_.size() ||
            !jstar_state_views_[static_cast<std::size_t>(iblk)].xi.is_bound() ||
            !jstar_state_views_[static_cast<std::size_t>(iblk)].eta.is_bound() ||
            !jstar_state_views_[static_cast<std::size_t>(iblk)].zeta.is_bound())
            throw std::logic_error("MercurySolver Jstar Workspace view is missing for an active block");

        const mpcns::discrete::SnapshotLumpedHodge hodge{
            geometry_snapshot_.block(static_cast<std::size_t>(iblk)).lumped_hodge};
        const mpcns::operators::ElectromagneticOperators electromagnetic{hodge};

        const mpcns::state::PrimalFaceFormStateView magnetic_flux{
            state_view_adapter_->view(state_b_xi, iblk),
            state_view_adapter_->view(state_b_eta, iblk),
            state_view_adapter_->view(state_b_zeta, iblk)};
        const JstarStateBlockViews &jstar = jstar_state_views_[static_cast<std::size_t>(iblk)];
        hodge.apply_star_2(magnetic_flux,
                           mpcns::state::DualEdgeFormStateView{jstar.xi, jstar.eta, jstar.zeta});

        // Jedge intentionally aliases its historical dual-face intermediate
        // and final primal-edge state through these three borrowed views.
        const mpcns::state::StateView j_xi = state_view_adapter_->view(state_j_xi, iblk);
        const mpcns::state::StateView j_eta = state_view_adapter_->view(state_j_eta, iblk);
        const mpcns::state::StateView j_zeta = state_view_adapter_->view(state_j_zeta, iblk);
        electromagnetic.compute_current_from_star_2(
            mpcns::state::DualEdgeFormStateView{jstar.xi, jstar.eta, jstar.zeta},
            mpcns::state::DualFaceFormStateView{j_xi, j_eta, j_zeta},
            mpcns::state::PrimalEdgeFormStateView{j_xi, j_eta, j_zeta});
    }
}

void MercurySolver::SetupJstarWorkspace_()
{
    if (!fld_ || !state_view_adapter_ || !state_workspace_)
        throw std::logic_error("MercurySolver Jstar Workspace requires initialized State services");

    jstar_state_views_.clear();
    jstar_state_views_.resize(static_cast<std::size_t>(fld_->num_blocks()));
    for (int iblk = 0; iblk < fld_->num_blocks(); ++iblk)
    {
        const bool has_xi = fld_->field(fid_.fid_Jstar.xi, iblk).is_allocated();
        const bool has_eta = fld_->field(fid_.fid_Jstar.eta, iblk).is_allocated();
        const bool has_zeta = fld_->field(fid_.fid_Jstar.zeta, iblk).is_allocated();
        if (!has_xi && !has_eta && !has_zeta)
            continue;
        if (!has_xi || !has_eta || !has_zeta)
            throw std::logic_error("MercurySolver Jstar Workspace requires a complete local triplet");

        const auto bind_jstar = [&](mpcns::state::StateHandle handle, int field_id)
        {
            state_workspace_->bind(handle, iblk,
                                   state_view_adapter_->view(handle.quantity(), iblk));
            // This one-time check proves that the production Workspace holds
            // the exact FieldStorage allocation subsequently consumed by the
            // legacy Hodge/DEC operator bridge.
            const mpcns::state::StateView view = state_workspace_->view(handle, iblk);
            FieldBlock &legacy = fld_->field(field_id, iblk);
            const Int3 lo = legacy.get_lo();
            if (&view(lo.i, lo.j, lo.k, 0) != &legacy(lo.i, lo.j, lo.k, 0))
                throw std::logic_error("MercurySolver Jstar Workspace lost zero-copy storage identity");
        };
        bind_jstar(state_jstar_xi_, fid_.fid_Jstar.xi);
        bind_jstar(state_jstar_eta_, fid_.fid_Jstar.eta);
        bind_jstar(state_jstar_zeta_, fid_.fid_Jstar.zeta);
        JstarStateBlockViews &cached = jstar_state_views_[static_cast<std::size_t>(iblk)];
        cached.xi = state_workspace_->view(state_jstar_xi_, iblk);
        cached.eta = state_workspace_->view(state_jstar_eta_, iblk);
        cached.zeta = state_workspace_->view(state_jstar_zeta_, iblk);
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
