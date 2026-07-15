#include "LunarSolver.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

void LunarSolver::AddIdealEdgeEMF_()
{
    for (int iblk = 0; iblk < fld_->num_blocks(); iblk++)
    {
        auto &Uplus = fld_->field(fid_.fid_U_plus, iblk);
        auto &UH = fld_->field(fid_.fid_U_H, iblk);
        if (!Uplus.is_allocated() || !UH.is_allocated())
            continue;

        auto &Bxi = fld_->field(fid_.fid_B.xi, iblk);
        auto &Beta = fld_->field(fid_.fid_B.eta, iblk);
        auto &Bzeta = fld_->field(fid_.fid_B.zeta, iblk);
        auto &Badd_xi = fld_->field(fid_.fid_Badd.xi, iblk);
        auto &Badd_eta = fld_->field(fid_.fid_Badd.eta, iblk);
        auto &Badd_zeta = fld_->field(fid_.fid_Badd.zeta, iblk);

        auto &Jac = fld_->field(fid_.fid_Jac, iblk);
        auto &JDxi = fld_->field(fid_.fid_metric.xi, iblk);
        auto &JDet = fld_->field(fid_.fid_metric.eta, iblk);
        auto &JDze = fld_->field(fid_.fid_metric.zeta, iblk);

        AssembleOneDirectionEMF_(1, fld_->field(fid_.fid_Eface.xi, iblk),
                                 Bxi, Beta, Bzeta,
                                 Badd_xi, Badd_eta, Badd_zeta,
                                 Jac, JDxi, JDet, JDze, Uplus);
        AssembleOneDirectionEMF_(2, fld_->field(fid_.fid_Eface.eta, iblk),
                                 Bxi, Beta, Bzeta,
                                 Badd_xi, Badd_eta, Badd_zeta,
                                 Jac, JDxi, JDet, JDze, Uplus);
        AssembleOneDirectionEMF_(3, fld_->field(fid_.fid_Eface.zeta, iblk),
                                 Bxi, Beta, Bzeta,
                                 Badd_xi, Badd_eta, Badd_zeta,
                                 Jac, JDxi, JDet, JDze, Uplus);
    }

    // Eface is a block-local work tensor: the field location is a face 2-form
    // axis, while its three stored components are covariant EMF components in
    // that block's computational basis.  A face-triplet halo exchange can
    // transform the outer face axis, but it cannot also transform those three
    // inner components correctly.  In particular, copying the component index
    // unchanged across a cubic-sphere panel seam corrupts the edge candidate.
    //
    // AssembleOneDirectionEMF_ therefore fills the one ghost layer consumed by
    // AssembleEdgeEMF_FromFaceE_Ideal_ from the already-synchronised U/B fields
    // and the receiver block's own metrics.  Do not exchange Eface itself.

    AssembleEdgeEMF_FromFaceE_Ideal_();
}

void LunarSolver::AssembleSingularEdgeEMF_NonHall_()
{
    if (!singular_edges_ || singular_edges_->empty()) return;

    // One state per unique real sector around each quotient edge:
    // u(3), B(3), rho, p, cell-centred non-Hall EMF 1-form.
    constexpr int nstate=9;
    constexpr double Cuct=0.5;
    const auto &edges=singular_edges_->entries();
    std::vector<std::size_t> offset(edges.size()+1,0);
    for(std::size_t n=0;n<edges.size();++n)
        offset[n+1]=offset[n]+edges[n].ordered_cell_centers.size();
    const std::size_t nsector=offset.back();
    std::vector<double> local(nstate*nsector,0.0),global(nstate*nsector,0.0);
    std::vector<int> local_count(nsector,0),global_count(nsector,0);
    std::vector<double> local_average(edges.size(),0.0),global_average(edges.size(),0.0);

    for(std::size_t n=0;n<edges.size();++n)
        for(const auto &inc:edges[n].local_incident_cells)
        {
            if(inc.sector_index<0) continue;
            const auto &c=inc.entity;
            FieldBlock &u=fld_->field(fid_.fid_U_plus,c.block);
            FieldBlock &b=fld_->field(fid_.fid_Bcell,c.block);
            FieldBlock &uH=fld_->field(fid_.fid_U_H,c.block);
            FieldBlock &pH=fld_->field(fid_.fid_PV_H,c.block);
            FieldBlock &jc=fld_->field(fid_.fid_Jcell,c.block);
            if(!u.is_allocated() || !b.is_allocated() || !uH.is_allocated() ||
               !pH.is_allocated() || !jc.is_allocated()) continue;
            const std::size_t s=offset[n]+static_cast<std::size_t>(inc.sector_index);
            double *q=&local[nstate*s];
            for(int d=0;d<3;++d) { q[d]+=u(c.i,c.j,c.k,d); q[3+d]+=b(c.i,c.j,c.k,d); }
            q[6]+=uH(c.i,c.j,c.k,0);
            q[7]+=pH(c.i,c.j,c.k,3);

            const double ex=-(u(c.i,c.j,c.k,1)*b(c.i,c.j,c.k,2)-
                              u(c.i,c.j,c.k,2)*b(c.i,c.j,c.k,1));
            const double ey=-(u(c.i,c.j,c.k,2)*b(c.i,c.j,c.k,0)-
                              u(c.i,c.j,c.k,0)*b(c.i,c.j,c.k,2));
            const double ez=-(u(c.i,c.j,c.k,0)*b(c.i,c.j,c.k,1)-
                              u(c.i,c.j,c.k,1)*b(c.i,c.j,c.k,0));
            double emf=ex*edges[n].canonical_edge_vector[0]+
                       ey*edges[n].canonical_edge_vector[1]+
                       ez*edges[n].canonical_edge_vector[2];

            // Non-ideal terms retain their physical-cell reconstruction.  The
            // only centered derivative is tangent to the real edge; collapsed
            // transverse ghosts are never used.
            if(ambipolar_control.enabled)
            {
                int im=c.i,jm=c.j,km=c.k,ip=c.i,jp=c.j,kp=c.k;
                const int ax=static_cast<int>(inc.source_alias.axis);
                if(ax==0){--im;++ip;} else if(ax==1){--jm;++jp;} else {--km;++kp;}
                const double dp=0.5*(pH(ip,jp,kp,3)-pH(im,jm,km,3));
                const NumInfo num=Hall_Num_Limiter(uH(c.i,c.j,c.k,0));
                const double pressure_coef=(rho_ref*U_ref)/(q_e*L_ref*B_ref*n_ref);
                emf-=pressure_coef*inc.source_orientation*dp/num.ne_eff;
            }
            const double jedge=jc(c.i,c.j,c.k,0)*edges[n].canonical_edge_vector[0]+
                               jc(c.i,c.j,c.k,1)*edges[n].canonical_edge_vector[1]+
                               jc(c.i,c.j,c.k,2)*edges[n].canonical_edge_vector[2];
            const double jmag=std::sqrt(jc(c.i,c.j,c.k,0)*jc(c.i,c.j,c.k,0)+
                                        jc(c.i,c.j,c.k,1)*jc(c.i,c.j,c.k,1)+
                                        jc(c.i,c.j,c.k,2)*jc(c.i,c.j,c.k,2));
            if(arti_resist_control.eta_max>0.0)
            {
                double a=(jmag-arti_resist_control.J_range_start)/
                    std::max(arti_resist_control.J_range_on-arti_resist_control.J_range_start,1.e-30);
                a=std::max(0.0,std::min(1.0,a));
                emf+=arti_resist_control.eta_max*a*a*(3.0-2.0*a)*jedge;
            }
            if(arti_resist_control.local_enabled && arti_resist_control.local_eta_max>0.0)
            {
                auto &cx=grd_->grids(c.block).dual_x;
                auto &cy=grd_->grids(c.block).dual_y;
                auto &cz=grd_->grids(c.block).dual_z;
                const double dx=cx(c.i+1,c.j+1,c.k+1)-arti_resist_control.local_center[0];
                const double dy=cy(c.i+1,c.j+1,c.k+1)-arti_resist_control.local_center[1];
                const double dz=cz(c.i+1,c.j+1,c.k+1)-arti_resist_control.local_center[2];
                const double rr=std::sqrt(dx*dx+dy*dy+dz*dz);
                if(arti_resist_control.local_r_cutoff<=0.0 || rr<arti_resist_control.local_r_cutoff)
                    emf+=arti_resist_control.local_eta_max*
                         std::exp(-rr/arti_resist_control.local_r_decay)*jedge;
            }
            q[8]+=emf;
            local_average[n]+=inc.weight*emf;
            ++local_count[s];
        }

    MPI_Allreduce(local.data(),global.data(),static_cast<int>(global.size()),
                  MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(local_count.data(),global_count.data(),static_cast<int>(global_count.size()),
                  MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(local_average.data(),global_average.data(),static_cast<int>(global_average.size()),
                  MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    auto dot=[](const double *a,const std::array<double,3> &b)
    { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; };
    auto fast_speed=[&](const double *q,const std::array<double,3> &normal)
    {
        // Use the actual positive density, matching the ideal-MHD CFL policy.
        const double rho=q[6];
        if(!(rho>0.0) || !std::isfinite(rho)) return 0.0;
        const double p=std::max(q[7],1.e-8);
        const double cs2=gamma_*p/rho;
        const double b2=q[3]*q[3]+q[4]*q[4]+q[5]*q[5];
        const double bn=q[3]*normal[0]+q[4]*normal[1]+q[5]*normal[2];
        const double va2=inver_MA2*b2/rho;
        const double van2=inver_MA2*bn*bn/rho;
        const double term=cs2+va2;
        const double disc=std::max(0.0,term*term-4.0*cs2*van2);
        const double cf=std::sqrt(std::max(0.0,0.5*(term+std::sqrt(disc))));
        return std::abs(dot(q,normal))+cf;
    };

    for(std::size_t n=0;n<edges.size();++n)
    {
        const auto &edge=edges[n];
        const std::size_t ns=edge.ordered_cell_centers.size();
        if(ns<3 || edge.owner.rank!=par_->GetInt("myid")) continue;
        int owner_sign=+1;
        for(const auto &alias:edge.aliases) if(alias.owner){owner_sign=alias.orientation;break;}
        const int fid=fid_.fid_E.at(static_cast<int>(edge.owner.axis)+1);
        if(singular_emf_mode_=="cell_average")
        {
            fld_->field(fid,edge.owner.block)(edge.owner.i,edge.owner.j,edge.owner.k,0)=
                owner_sign*global_average[n];
            continue;
        }
        const double L=std::sqrt(edge.canonical_edge_vector[0]*edge.canonical_edge_vector[0]+
                                 edge.canonical_edge_vector[1]*edge.canonical_edge_vector[1]+
                                 edge.canonical_edge_vector[2]*edge.canonical_edge_vector[2]);
        if(!(L>0.0)) continue;
        const std::array<double,3> t{{edge.canonical_edge_vector[0]/L,
                                      edge.canonical_edge_vector[1]/L,
                                      edge.canonical_edge_vector[2]/L}};
        double emf_sum=0.0,weight_sum=0.0;
        for(std::size_t m=0;m<ns;++m)
        {
            const std::size_t ia=offset[n]+m,ib=offset[n]+(m+1)%ns;
            if(global_count[ia]!=1 || global_count[ib]!=1)
                throw std::runtime_error("Singular EMF sector does not map to one unique real cell");
            const double *qa=&global[nstate*ia],*qb=&global[nstate*ib];
            std::array<double,3> normal{{
                edge.ordered_cell_centers[(m+1)%ns][0]-edge.ordered_cell_centers[m][0],
                edge.ordered_cell_centers[(m+1)%ns][1]-edge.ordered_cell_centers[m][1],
                edge.ordered_cell_centers[(m+1)%ns][2]-edge.ordered_cell_centers[m][2]}};
            const double nt=normal[0]*t[0]+normal[1]*t[1]+normal[2]*t[2];
            for(int d=0;d<3;++d) normal[d]-=nt*t[d];
            const double dn=std::sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
            if(!(dn>0.0)) continue;
            for(double &v:normal) v/=dn;
            // Transverse magnetic component carried by this incident face.
            const std::array<double,3> side{{normal[1]*t[2]-normal[2]*t[1],
                                             normal[2]*t[0]-normal[0]*t[2],
                                             normal[0]*t[1]-normal[1]*t[0]}};
            const double dBs=(qb[3]-qa[3])*side[0]+(qb[4]-qa[4])*side[1]+(qb[5]-qa[5])*side[2];
            const double a=std::max(fast_speed(qa,normal),fast_speed(qb,normal));
            const double pair_emf=0.5*(qa[8]+qb[8])+Cuct*0.5*a*dBs*L;
            emf_sum+=dn*pair_emf;
            weight_sum+=dn;
        }
        if(!(weight_sum>0.0)) continue;
        fld_->field(fid,edge.owner.block)(edge.owner.i,edge.owner.j,edge.owner.k,0)=
            owner_sign*emf_sum/weight_sum;
    }
}

//=========================================================================
// Face candidates to the unique shared edge EMF.
// Eface_* stores candidate components (e_xi,e_eta,e_zeta), not Cartesian E.
//
// This version adds a simple UCT transverse upwind correction.
// It keeps the CT structure because the final B update is still curl(E_edge).
void LunarSolver::AssembleEdgeEMF_FromFaceE_Ideal_()
{
    constexpr double jac_floor = 1.0e-30;

    // UCT correction strength.
    // Start with 0.5. If still too weak, try 1.0.
    // If it becomes noisy, reduce to 0.25.
    constexpr double Cuct = 0.5;

    auto shift_cell = [](int axis, int side, int &i, int &j, int &k)
    {
        if (axis == 0)
            i += side;
        else if (axis == 1)
            j += side;
        else
            k += side;
    };

    Int3 sub, sup;

    for (int iblk = 0; iblk < fld_->num_blocks(); iblk++)
    {
        auto &Uplus = fld_->field(fid_.fid_U_plus, iblk);

        auto &Bxi = fld_->field(fid_.fid_B.xi, iblk);
        auto &Beta = fld_->field(fid_.fid_B.eta, iblk);
        auto &Bzeta = fld_->field(fid_.fid_B.zeta, iblk);

        auto &Jac = fld_->field(fid_.fid_Jac, iblk);
        auto &JDxi = fld_->field(fid_.fid_metric.xi, iblk);
        auto &JDet = fld_->field(fid_.fid_metric.eta, iblk);
        auto &JDze = fld_->field(fid_.fid_metric.zeta, iblk);

        if (!Uplus.is_allocated())
            continue;

        auto metric_component = [&](int axis, int i, int j, int k, int comp) -> double
        {
            if (axis == 0)
                return JDxi(i, j, k, comp);
            if (axis == 1)
                return JDet(i, j, k, comp);
            return JDze(i, j, k, comp);
        };

        auto cell_u_contra = [&](int axis, int i, int j, int k) -> double
        {
            int ip = i, jp = j, kp = k;
            shift_cell(axis, 1, ip, jp, kp);

            const double kx = 0.5 * (metric_component(axis, i, j, k, 0) +
                                     metric_component(axis, ip, jp, kp, 0));

            const double ky = 0.5 * (metric_component(axis, i, j, k, 1) +
                                     metric_component(axis, ip, jp, kp, 1));

            const double kz = 0.5 * (metric_component(axis, i, j, k, 2) +
                                     metric_component(axis, ip, jp, kp, 2));

            const double inv_jac = 1.0 / (std::abs(Jac(i, j, k, 0)) + jac_floor);

            return (kx * Uplus(i, j, k, 0) +
                    ky * Uplus(i, j, k, 1) +
                    kz * Uplus(i, j, k, 2)) *
                   inv_jac;
        };

        auto max_abs_u4 = [&](int axis,
                              int i0, int j0, int k0,
                              int i1, int j1, int k1,
                              int i2, int j2, int k2,
                              int i3, int j3, int k3) -> double
        {
            double a = 0.0;
            a = std::max(a, std::abs(cell_u_contra(axis, i0, j0, k0)));
            a = std::max(a, std::abs(cell_u_contra(axis, i1, j1, k1)));
            a = std::max(a, std::abs(cell_u_contra(axis, i2, j2, k2)));
            a = std::max(a, std::abs(cell_u_contra(axis, i3, j3, k3)));
            return a;
        };

        //=============================================================
        // E_xi edge
        //
        // E_xi = central
        //      + 0.5 * a_eta  * d_beta_zeta / d_eta
        //      - 0.5 * a_zeta * d_beta_eta  / d_zeta
        //
        // Here beta means the evolved face 2-form B, not Badd.
        {
            auto &Exi = fld_->field(fid_.fid_E.xi, iblk);
            auto &E_face_eta = fld_->field(fid_.fid_Eface.eta, iblk);
            auto &E_face_zeta = fld_->field(fid_.fid_Eface.zeta, iblk);

            sub = Exi.inner_lo();
            sup = Exi.inner_hi();

            for (int i = sub.i; i < sup.i; i++)
                for (int j = sub.j; j < sup.j; j++)
                    for (int k = sub.k; k < sup.k; k++)
                    {
                        const double Ec = 0.25 * (E_face_eta(i, j, k, 0) +
                                                  E_face_eta(i, j, k - 1, 0) +
                                                  E_face_zeta(i, j, k, 0) +
                                                  E_face_zeta(i, j - 1, k, 0));

                        const double a_eta = max_abs_u4(
                            1,
                            i, j, k,
                            i, j - 1, k,
                            i, j, k - 1,
                            i, j - 1, k - 1);

                        const double a_zeta = max_abs_u4(
                            2,
                            i, j, k,
                            i, j - 1, k,
                            i, j, k - 1,
                            i, j - 1, k - 1);

                        const double dBzeta_eta =
                            Bzeta(i, j, k, 0) -
                            Bzeta(i, j - 1, k, 0);

                        const double dBeta_zeta =
                            Beta(i, j, k, 0) -
                            Beta(i, j, k - 1, 0);

                        Exi(i, j, k, 0) =
                            Ec + Cuct * (0.5 * a_eta * dBzeta_eta - 0.5 * a_zeta * dBeta_zeta);
                    }
        }

        //=============================================================
        // E_eta edge
        //
        // E_eta = central
        //       + 0.5 * a_zeta * d_beta_xi   / d_zeta
        //       - 0.5 * a_xi   * d_beta_zeta / d_xi
        {
            auto &Eeta = fld_->field(fid_.fid_E.eta, iblk);
            auto &E_face_xi = fld_->field(fid_.fid_Eface.xi, iblk);
            auto &E_face_zeta = fld_->field(fid_.fid_Eface.zeta, iblk);

            sub = Eeta.inner_lo();
            sup = Eeta.inner_hi();

            for (int i = sub.i; i < sup.i; i++)
                for (int j = sub.j; j < sup.j; j++)
                    for (int k = sub.k; k < sup.k; k++)
                    {
                        const double Ec = 0.25 * (E_face_xi(i, j, k, 1) +
                                                  E_face_xi(i, j, k - 1, 1) +
                                                  E_face_zeta(i, j, k, 1) +
                                                  E_face_zeta(i - 1, j, k, 1));

                        const double a_zeta = max_abs_u4(
                            2,
                            i, j, k,
                            i - 1, j, k,
                            i, j, k - 1,
                            i - 1, j, k - 1);

                        const double a_xi = max_abs_u4(
                            0,
                            i, j, k,
                            i - 1, j, k,
                            i, j, k - 1,
                            i - 1, j, k - 1);

                        const double dBxi_zeta =
                            Bxi(i, j, k, 0) -
                            Bxi(i, j, k - 1, 0);

                        const double dBzeta_xi =
                            Bzeta(i, j, k, 0) -
                            Bzeta(i - 1, j, k, 0);

                        Eeta(i, j, k, 0) =
                            Ec + Cuct * (0.5 * a_zeta * dBxi_zeta - 0.5 * a_xi * dBzeta_xi);
                    }
        }

        //=============================================================
        // E_zeta edge
        //
        // E_zeta = central
        //        + 0.5 * a_xi  * d_beta_eta / d_xi
        //        - 0.5 * a_eta * d_beta_xi  / d_eta
        {
            auto &Ezeta = fld_->field(fid_.fid_E.zeta, iblk);
            auto &E_face_xi = fld_->field(fid_.fid_Eface.xi, iblk);
            auto &E_face_eta = fld_->field(fid_.fid_Eface.eta, iblk);

            sub = Ezeta.inner_lo();
            sup = Ezeta.inner_hi();

            for (int i = sub.i; i < sup.i; i++)
                for (int j = sub.j; j < sup.j; j++)
                    for (int k = sub.k; k < sup.k; k++)
                    {
                        const double Ec = 0.25 * (E_face_xi(i, j, k, 2) +
                                                  E_face_xi(i, j - 1, k, 2) +
                                                  E_face_eta(i, j, k, 2) +
                                                  E_face_eta(i - 1, j, k, 2));

                        const double a_xi = max_abs_u4(
                            0,
                            i, j, k,
                            i - 1, j, k,
                            i, j - 1, k,
                            i - 1, j - 1, k);

                        const double a_eta = max_abs_u4(
                            1,
                            i, j, k,
                            i - 1, j, k,
                            i, j - 1, k,
                            i - 1, j - 1, k);

                        const double dBeta_xi =
                            Beta(i, j, k, 0) -
                            Beta(i - 1, j, k, 0);

                        const double dBxi_eta =
                            Bxi(i, j, k, 0) -
                            Bxi(i, j - 1, k, 0);

                        Ezeta(i, j, k, 0) =
                            Ec + Cuct * (0.5 * a_xi * dBeta_xi - 0.5 * a_eta * dBxi_eta);
                    }
        }
    }
}

//=========================================================================
// Face candidates to the unique shared edge EMF.
// Eface_* stores candidate components (e_xi,e_eta,e_zeta), not Cartesian E.
// void LunarSolver::AssembleEdgeEMF_FromFaceE_Ideal_()
// {
//     Int3 sub, sup;
//     for (int iblk = 0; iblk < fld_->num_blocks(); iblk++)
//     {
//         {
//             auto &Exi = fld_->field(fid_.fid_E.xi, iblk);
//             auto &E_face_eta = fld_->field(fid_.fid_Eface.eta, iblk);
//             auto &E_face_zeta = fld_->field(fid_.fid_Eface.zeta, iblk);
//             sub = Exi.inner_lo();
//             sup = Exi.inner_hi();
//             for (int i = sub.i; i < sup.i; i++)
//                 for (int j = sub.j; j < sup.j; j++)
//                     for (int k = sub.k; k < sup.k; k++)
//                     {
//                         Exi(i, j, k, 0) = 0.25 * (
//                             E_face_eta(i, j, k, 0) + E_face_eta(i, j, k - 1, 0) +
//                             E_face_zeta(i, j, k, 0) + E_face_zeta(i, j - 1, k, 0));
//                     }
//         }

//         {
//             auto &Eeta = fld_->field(fid_.fid_E.eta, iblk);
//             auto &E_face_xi = fld_->field(fid_.fid_Eface.xi, iblk);
//             auto &E_face_zeta = fld_->field(fid_.fid_Eface.zeta, iblk);
//             sub = Eeta.inner_lo();
//             sup = Eeta.inner_hi();
//             for (int i = sub.i; i < sup.i; i++)
//                 for (int j = sub.j; j < sup.j; j++)
//                     for (int k = sub.k; k < sup.k; k++)
//                     {
//                         Eeta(i, j, k, 0) = 0.25 * (
//                             E_face_xi(i, j, k, 1) + E_face_xi(i, j, k - 1, 1) +
//                             E_face_zeta(i, j, k, 1) + E_face_zeta(i - 1, j, k, 1));
//                     }
//         }

//         {
//             auto &Ezeta = fld_->field(fid_.fid_E.zeta, iblk);
//             auto &E_face_xi = fld_->field(fid_.fid_Eface.xi, iblk);
//             auto &E_face_eta = fld_->field(fid_.fid_Eface.eta, iblk);
//             sub = Ezeta.inner_lo();
//             sup = Ezeta.inner_hi();
//             for (int i = sub.i; i < sup.i; i++)
//                 for (int j = sub.j; j < sup.j; j++)
//                     for (int k = sub.k; k < sup.k; k++)
//                     {
//                         Ezeta(i, j, k, 0) = 0.25 * (
//                             E_face_xi(i, j, k, 2) + E_face_xi(i, j - 1, k, 2) +
//                             E_face_eta(i, j, k, 2) + E_face_eta(i - 1, j, k, 2));
//                     }
//         }
//     }
// }

//=========================================================================
void LunarSolver::AssembleOneDirectionEMF_(
    int dir, FieldBlock &E_face,
    FieldBlock &Bxi, FieldBlock &Beta, FieldBlock &Bzeta,
    FieldBlock &Badd_xi, FieldBlock &Badd_eta, FieldBlock &Badd_zeta,
    FieldBlock &Jac,
    FieldBlock &JDxi, FieldBlock &JDet, FieldBlock &JDze,
    FieldBlock &Uplus)
{
    constexpr double jac_floor = 1.0e-30;

    auto shift_cell = [](int axis, int side, int &i, int &j, int &k)
    {
        if (axis == 0)
            i += side;
        else if (axis == 1)
            j += side;
        else
            k += side;
    };

    auto face_beta_ind = [&](int axis, int i, int j, int k) -> double
    {
        if (axis == 0)
            return Bxi(i, j, k, 0);
        if (axis == 1)
            return Beta(i, j, k, 0);
        return Bzeta(i, j, k, 0);
    };

    auto face_beta_add = [&](int axis, int i, int j, int k) -> double
    {
        if (axis == 0)
            return Badd_xi(i, j, k, 0);
        if (axis == 1)
            return Badd_eta(i, j, k, 0);
        return Badd_zeta(i, j, k, 0);
    };

    auto face_beta_total = [&](int axis, int i, int j, int k) -> double
    {
        return face_beta_ind(axis, i, j, k) + face_beta_add(axis, i, j, k);
    };

    auto cell_beta_ind = [&](int axis, int i, int j, int k) -> double
    {
        int ip = i, jp = j, kp = k;
        shift_cell(axis, 1, ip, jp, kp);
        return 0.5 * (face_beta_ind(axis, i, j, k) + face_beta_ind(axis, ip, jp, kp));
    };

    auto cell_beta_total = [&](int axis, int i, int j, int k) -> double
    {
        int ip = i, jp = j, kp = k;
        shift_cell(axis, 1, ip, jp, kp);
        return 0.5 * (face_beta_total(axis, i, j, k) + face_beta_total(axis, ip, jp, kp));
    };

    auto metric_component = [&](int axis, int i, int j, int k, int comp) -> double
    {
        if (axis == 0)
            return JDxi(i, j, k, comp);
        if (axis == 1)
            return JDet(i, j, k, comp);
        return JDze(i, j, k, comp);
    };

    auto cell_u_contra = [&](int axis, int i, int j, int k) -> double
    {
        int ip = i, jp = j, kp = k;
        shift_cell(axis, 1, ip, jp, kp);

        const double kx = 0.5 * (metric_component(axis, i, j, k, 0) + metric_component(axis, ip, jp, kp, 0));
        const double ky = 0.5 * (metric_component(axis, i, j, k, 1) + metric_component(axis, ip, jp, kp, 1));
        const double kz = 0.5 * (metric_component(axis, i, j, k, 2) + metric_component(axis, ip, jp, kp, 2));
        const double inv_jac = 1.0 / (std::abs(Jac(i, j, k, 0)) + jac_floor);

        return (kx * Uplus(i, j, k, 0) +
                ky * Uplus(i, j, k, 1) +
                kz * Uplus(i, j, k, 2)) *
               inv_jac;
    };

    auto flux = [&](int beta_axis, int flux_axis,
                    int iL, int jL, int kL,
                    int iR, int jR, int kR,
                    double beta_flux_face) -> double
    {
        const double betaTotalL = cell_beta_total(beta_axis, iL, jL, kL);
        const double betaTotalR = cell_beta_total(beta_axis, iR, jR, kR);
        const double betaIndL = cell_beta_ind(beta_axis, iL, jL, kL);
        const double betaIndR = cell_beta_ind(beta_axis, iR, jR, kR);

        const double uFluxL = cell_u_contra(flux_axis, iL, jL, kL);
        const double uFluxR = cell_u_contra(flux_axis, iR, jR, kR);
        const double uBetaL = cell_u_contra(beta_axis, iL, jL, kL);
        const double uBetaR = cell_u_contra(beta_axis, iR, jR, kR);

        const double GL = uFluxL * betaTotalL - uBetaL * beta_flux_face;
        const double GR = uFluxR * betaTotalR - uBetaR * beta_flux_face;
        const double radius = std::max(std::abs(uFluxL), std::abs(uFluxR));

        return 0.5 * (GL + GR) - 0.5 * radius * (betaIndR - betaIndL);
    };

    const int flux_axis = dir - 1;
    if (flux_axis < 0 || flux_axis > 2)
        throw std::runtime_error("AssembleOneDirectionEMF_: invalid direction");

    Int3 sub = E_face.inner_lo();
    Int3 sup = E_face.inner_hi();

    // Edge assembly reads the two face candidates adjacent to an inner edge,
    // hence it needs one locally evaluated Eface ghost layer.  The primary
    // state and geometry fields carry a wider runtime halo (normally four), so
    // this expansion leaves the flux stencil inside their allocated ranges.
    const Int3 alloc_lo = E_face.get_lo();
    const Int3 alloc_hi = E_face.get_hi();
    sub.i = std::max(sub.i - 1, alloc_lo.i + 1);
    sub.j = std::max(sub.j - 1, alloc_lo.j + 1);
    sub.k = std::max(sub.k - 1, alloc_lo.k + 1);
    sup.i = std::min(sup.i + 1, alloc_hi.i - 1);
    sup.j = std::min(sup.j + 1, alloc_hi.j - 1);
    sup.k = std::min(sup.k + 1, alloc_hi.k - 1);

    for (int i = sub.i; i < sup.i; ++i)
        for (int j = sub.j; j < sup.j; ++j)
            for (int k = sub.k; k < sup.k; ++k)
            {
                int iL = i, jL = j, kL = k;
                int iR = i, jR = j, kR = k;
                shift_cell(flux_axis, -1, iL, jL, kL);

                const double beta_flux_face = face_beta_total(flux_axis, i, j, k);
                double e[3] = {0.0, 0.0, 0.0};

                if (flux_axis == 0)
                {
                    const double G_eta_xi = flux(1, 0, iL, jL, kL, iR, jR, kR, beta_flux_face);
                    const double G_zeta_xi = flux(2, 0, iL, jL, kL, iR, jR, kR, beta_flux_face);
                    e[1] = G_zeta_xi;
                    e[2] = -G_eta_xi;
                }
                else if (flux_axis == 1)
                {
                    const double G_xi_eta = flux(0, 1, iL, jL, kL, iR, jR, kR, beta_flux_face);
                    const double G_zeta_eta = flux(2, 1, iL, jL, kL, iR, jR, kR, beta_flux_face);
                    e[0] = -G_zeta_eta;
                    e[2] = G_xi_eta;
                }
                else
                {
                    const double G_xi_zeta = flux(0, 2, iL, jL, kL, iR, jR, kR, beta_flux_face);
                    const double G_eta_zeta = flux(1, 2, iL, jL, kL, iR, jR, kR, beta_flux_face);
                    e[0] = G_eta_zeta;
                    e[1] = -G_xi_zeta;
                }

                E_face(i, j, k, 0) = e[0];
                E_face(i, j, k, 1) = e[1];
                E_face(i, j, k, 2) = e[2];
            }
}
