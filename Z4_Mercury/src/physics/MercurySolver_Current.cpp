#include "MercurySolver.h"

#include "mpcns/compat/legacy_geometry_field_ids.hpp"
#include "mpcns/constraints/least_squares_3.hpp"
#include "mpcns/discrete/lumped_hodge.hpp"
#include "mpcns/operators/electromagnetic_operators.hpp"
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

        const mpcns::discrete::LumpedHodge hodge{
            mpcns::compat::bind_legacy_lumped_hodge_coefficients(*fld_, iblk)};
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
        constexpr double eps_len = 1.0e-14;

        for (const auto &p : topo_->physical_patches)
        {
            if (p.bc_name != "Pole")
                continue;

            const int dir = std::abs(p.direction);
            if (dir != 1 && dir != 2)
                continue;

            const int ib = p.this_block;
            if (ib < 0 || ib >= nblock)
                continue;

            auto &Jcell = fld_->field(fid_.fid_Jcell, ib);
            auto &Jxi = fld_->field(fid_.fid_J.xi, ib);
            auto &Jeta = fld_->field(fid_.fid_J.eta, ib);
            auto &Jzeta = fld_->field(fid_.fid_J.zeta, ib);
            auto &dr_xi = fld_->field(fid_.Edge_dr.xi, ib);
            auto &dr_eta = fld_->field(fid_.Edge_dr.eta, ib);
            auto &dr_zeta = fld_->field(fid_.Edge_dr.zeta, ib);

            if (!Jcell.is_allocated() || !Jxi.is_allocated() ||
                !Jeta.is_allocated() || !Jzeta.is_allocated() ||
                !dr_xi.is_allocated() || !dr_eta.is_allocated() ||
                !dr_zeta.is_allocated())
                continue;

            const bool high_side = (p.direction > 0);
            const Box3 &node = p.this_box_node;

            auto push_edge = [&](FieldBlock &Je, FieldBlock &dr,
                                 int i, int j, int k,
                                 double A[3][3], double b[3], int &cnt)
            {
                const double dx = dr(i, j, k, 0);
                const double dy = dr(i, j, k, 1);
                const double dz = dr(i, j, k, 2);
                const double L = std::sqrt(dx * dx + dy * dy + dz * dz);

                if (L < eps_len)
                    return;

                const double tau[3] = {dx / L, dy / L, dz / L};
                const double y = Je(i, j, k, 0) / L;

                for (int a = 0; a < 3; ++a)
                {
                    b[a] += y * tau[a];
                    for (int c = 0; c < 3; ++c)
                        A[a][c] += tau[a] * tau[c];
                }

                ++cnt;
            };

            auto write_ring = [&](int ic, int jc,
                                  double A[3][3], double b[3], int cnt,
                                  int klo, int khi)
            {
                if (cnt <= 0)
                    return;

                double Jp[3] = {0.0, 0.0, 0.0};
                if (!mpcns::constraints::solve_regularized_least_squares_3(A, b, Jp))
                    return;

                for (int k = klo; k < khi; ++k)
                {
                    Jcell(ic, jc, k, 0) = Jp[0];
                    Jcell(ic, jc, k, 1) = Jp[1];
                    Jcell(ic, jc, k, 2) = Jp[2];
                }
            };

            if (dir == 1)
            {
                const int ic = high_side ? (node.lo.i - 1) : node.lo.i;
                const int it = high_side ? ic : (ic + 1);

                for (int j = node.lo.j; j < node.hi.j - 1; ++j)
                {
                    double A[3][3] = {
                        {0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0}};
                    double b[3] = {0.0, 0.0, 0.0};
                    int cnt = 0;

                    for (int k = node.lo.k; k < node.hi.k - 1; ++k)
                    {
                        // norm edges on the Pole layer
                        push_edge(Jxi, dr_xi, ic, j, k, A, b, cnt);
                        push_edge(Jxi, dr_xi, ic, j + 1, k, A, b, cnt);
                        push_edge(Jxi, dr_xi, ic, j, k + 1, A, b, cnt);
                        push_edge(Jxi, dr_xi, ic, j + 1, k + 1, A, b, cnt);

                        // axis/rotate edges one layer away from the Pole
                        push_edge(Jeta, dr_eta, it, j, k, A, b, cnt);
                        push_edge(Jeta, dr_eta, it, j, k + 1, A, b, cnt);
                        push_edge(Jzeta, dr_zeta, it, j, k, A, b, cnt);
                        push_edge(Jzeta, dr_zeta, it, j + 1, k, A, b, cnt);
                    }

                    write_ring(ic, j, A, b, cnt, node.lo.k, node.hi.k - 1);
                }
            }
            else
            {
                const int jc = high_side ? (node.lo.j - 1) : node.lo.j;
                const int jt = high_side ? jc : (jc + 1);

                for (int i = node.lo.i; i < node.hi.i - 1; ++i)
                {
                    double A[3][3] = {
                        {0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0}};
                    double b[3] = {0.0, 0.0, 0.0};
                    int cnt = 0;

                    for (int k = node.lo.k; k < node.hi.k - 1; ++k)
                    {
                        // norm edges on the Pole layer
                        push_edge(Jeta, dr_eta, i, jc, k, A, b, cnt);
                        push_edge(Jeta, dr_eta, i + 1, jc, k, A, b, cnt);
                        push_edge(Jeta, dr_eta, i, jc, k + 1, A, b, cnt);
                        push_edge(Jeta, dr_eta, i + 1, jc, k + 1, A, b, cnt);

                        // axis/rotate edges one layer away from the Pole
                        push_edge(Jxi, dr_xi, i, jt, k, A, b, cnt);
                        push_edge(Jxi, dr_xi, i, jt, k + 1, A, b, cnt);
                        push_edge(Jzeta, dr_zeta, i, jt, k, A, b, cnt);
                        push_edge(Jzeta, dr_zeta, i + 1, jt, k, A, b, cnt);
                    }

                    write_ring(i, jc, A, b, cnt, node.lo.k, node.hi.k - 1);
                }
            }
        }
    }

    mercury_bound_.Sync("J_cell");

}
