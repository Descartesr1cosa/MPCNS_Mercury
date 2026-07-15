#include "LunarSolver.h"

#include "1_grid/1_MPCNS_Grid.h"
#include "0_basic/MPI_WRAPPER.h"

void LunarSolver::CanonicalizeSharedFaceGeometry_()
{
    if (!topo_ || topo_->faces.classes.empty())
        return;

    const int myid = par_->GetInt("myid");
    std::vector<const TOPO::EquivClass *> shared;
    for (const auto &cls : topo_->faces.classes)
        if (cls.global_id >= 0 && cls.members.size() > 1)
            shared.push_back(&cls);
    std::sort(shared.begin(), shared.end(), [](const auto *a, const auto *b)
              { return a->global_id < b->global_id; });
    if (shared.empty())
        return;

    std::vector<double> local_owner(3 * shared.size(), 0.0);
    std::vector<double> global_owner(3 * shared.size(), 0.0);
    for (std::size_t n = 0; n < shared.size(); ++n)
    {
        const auto &owner = shared[n]->owner;
        const auto &f = owner.entity;
        if (f.rank != myid)
            continue;
        FieldBlock &S = fld_->field(fid_.fid_metric.at(static_cast<int>(f.axis) + 1), f.block);
        for (int q = 0; q < 3; ++q)
            local_owner[3 * n + q] = owner.orient_sign * S(f.i, f.j, f.k, q);
    }
    MPI_Allreduce(local_owner.data(), global_owner.data(),
                  static_cast<int>(global_owner.size()), MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);

    double mismatch_local = 0.0;
    for (std::size_t n = 0; n < shared.size(); ++n)
    {
        const auto &canonical = global_owner;
        const double cmag = std::sqrt(canonical[3*n]*canonical[3*n] +
                                      canonical[3*n+1]*canonical[3*n+1] +
                                      canonical[3*n+2]*canonical[3*n+2]);
        for (const auto &member : shared[n]->members)
        {
            const auto &f = member.entity;
            if (f.rank != myid)
                continue;
            const int ax = static_cast<int>(f.axis);
            FieldBlock &S = fld_->field(fid_.fid_metric.at(ax + 1), f.block);
            double d2 = 0.0;
            for (int q = 0; q < 3; ++q)
            {
                const double old_canonical = member.orient_sign * S(f.i, f.j, f.k, q);
                const double d = old_canonical - canonical[3*n+q];
                d2 += d*d;
                S(f.i, f.j, f.k, q) = member.orient_sign * canonical[3*n+q];
            }
            mismatch_local = std::max(mismatch_local,
                                      std::sqrt(d2) / std::max(cmag, 1.0e-300));

            FieldBlock &area = fld_->field(fid_.Face_Area.at(ax + 1), f.block);
            FieldBlock &dlstar = fld_->field(fid_.Face_dlstar.at(ax + 1), f.block);
            FieldBlock &beta = fld_->field(fid_.Face_beta.at(ax + 1), f.block);
            const double amag = cmag;
            area(f.i, f.j, f.k, 0) = amag;
            beta(f.i, f.j, f.k, 0) = amag > 1.0e-300
                                         ? dlstar(f.i, f.j, f.k, 0) / amag
                                         : 0.0;
        }
    }

    // The closure residual is the exact constant-vector free-stream residual
    // of the finite-volume metric, independent of the physical flux function.
    double closure_local = 0.0;
    for (int ib = 0; ib < fld_->num_blocks(); ++ib)
    {
        FieldBlock &jac = fld_->field(fid_.fid_Jac, ib);
        FieldBlock &Sx = fld_->field(fid_.fid_metric.xi, ib);
        FieldBlock &Se = fld_->field(fid_.fid_metric.eta, ib);
        FieldBlock &Sz = fld_->field(fid_.fid_metric.zeta, ib);
        const Int3 lo = jac.inner_lo(), hi = jac.inner_hi();
        for (int i=lo.i;i<hi.i;++i) for (int j=lo.j;j<hi.j;++j) for (int k=lo.k;k<hi.k;++k)
        {
            double r[3] = {0.0,0.0,0.0};
            double denom = 0.0;
            for (int q=0;q<3;++q)
            {
                r[q] = Sx(i+1,j,k,q)-Sx(i,j,k,q) +
                       Se(i,j+1,k,q)-Se(i,j,k,q) +
                       Sz(i,j,k+1,q)-Sz(i,j,k,q);
            }
            auto add_norm = [&](FieldBlock &S,int ii,int jj,int kk)
            {
                double a=0.0; for(int q=0;q<3;++q) a+=S(ii,jj,kk,q)*S(ii,jj,kk,q);
                denom+=std::sqrt(a);
            };
            add_norm(Sx,i,j,k); add_norm(Sx,i+1,j,k);
            add_norm(Se,i,j,k); add_norm(Se,i,j+1,k);
            add_norm(Sz,i,j,k); add_norm(Sz,i,j,k+1);
            const double rn=std::sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
            closure_local=std::max(closure_local,rn/std::max(denom,1.0e-300));
        }
    }
    double diag_local[2] = {mismatch_local, closure_local};
    double diag_global[2] = {0.0,0.0};
    PARALLEL::mpi_max(diag_local, diag_global, 2);
    if (myid == 0)
        std::cout << "[Lunar][Geometry] shared_faces=" << shared.size()
                  << " alias_rel_before=" << diag_global[0]
                  << " cell_closure_rel_after=" << diag_global[1] << "\n";
}

void LunarSolver::SetupHallFaceScratch_()
{
    const int nb = fld_->num_blocks();
    hall_face_scratch_.clear();
    hall_face_scratch_.resize(nb);

    for (int ib = 0; ib < nb; ++ib)
    {
        auto &Bc = fld_->field(fid_.fid_Bcell, ib);
        if (!Bc.is_allocated())
            continue;

        const Int3 lo = Bc.get_lo();
        const Int3 hi = Bc.get_hi();
        const int ghost = -lo.i;
        if (lo.j != -ghost || lo.k != -ghost)
            throw std::runtime_error("SetupHallFaceScratch_: inconsistent ghost indexing.");

        auto &buf = hall_face_scratch_[ib];
        buf.Ehc.SetSize(hi.i - lo.i, hi.j - lo.j, hi.k - lo.k, ghost, 3);
        buf.beta.SetSize(hi.i - lo.i, hi.j - lo.j, hi.k - lo.k, ghost);
    }
}

void LunarSolver::SetupCellReconstructionWeights_()
{
    constexpr double eps = 1.0e-30;
    constexpr double eig_tol = 1.0e-10;

    auto eigen_sym3 = [](double A[3][3], double V[3][3], double d[3])
    {
        for (int a=0;a<3;++a) for (int b=0;b<3;++b) V[a][b]=(a==b)?1.0:0.0;
        for (int it=0;it<18;++it)
        {
            for (int p=0;p<2;++p) for (int q=p+1;q<3;++q)
            {
                const double apq=A[p][q];
                if (std::abs(apq) <= 1.0e-16*(std::abs(A[p][p])+std::abs(A[q][q])+1.0)) continue;
                const double phi=0.5*std::atan2(2.0*apq,A[q][q]-A[p][p]);
                const double c=std::cos(phi), s=std::sin(phi);
                const double app=A[p][p], aqq=A[q][q];
                for(int k=0;k<3;++k) if(k!=p&&k!=q)
                {
                    const double akp=A[k][p], akq=A[k][q];
                    A[k][p]=A[p][k]=c*akp-s*akq;
                    A[k][q]=A[q][k]=s*akp+c*akq;
                }
                A[p][p]=c*c*app-2.0*s*c*apq+s*s*aqq;
                A[q][q]=s*s*app+2.0*s*c*apq+c*c*aqq;
                A[p][q]=A[q][p]=0.0;
                for(int k=0;k<3;++k)
                {
                    const double vkp=V[k][p], vkq=V[k][q];
                    V[k][p]=c*vkp-s*vkq;
                    V[k][q]=s*vkp+c*vkq;
                }
            }
        }
        d[0]=A[0][0]; d[1]=A[1][1]; d[2]=A[2][2];
    };

    // Face 2-form fluxes -> Cartesian Bcell (six oriented cell faces).
    for (int ib=0;ib<fld_->num_blocks();++ib)
    {
        auto &Bc=fld_->field(fid_.fid_Bcell,ib);
        auto &W=fld_->field(fid_.fid_Bcell_from_Bface_w,ib);
        auto &Sxi=fld_->field(fid_.fid_metric.xi,ib);
        auto &Set=fld_->field(fid_.fid_metric.eta,ib);
        auto &Sze=fld_->field(fid_.fid_metric.zeta,ib);
        if(!Bc.is_allocated()||!W.is_allocated()||!Sxi.is_allocated()||!Set.is_allocated()||!Sze.is_allocated()) continue;
        auto &b=grd_->grids(ib);
        const Int3 lo=Bc.inner_lo(), hi=Bc.inner_hi();
        for(int i=lo.i;i<hi.i;++i) for(int j=lo.j;j<hi.j;++j) for(int k=lo.k;k<hi.k;++k)
        {
            double s[6][3];
            for(int a=0;a<3;++a)
            {
                s[0][a]=-Sxi(i,j,k,a); s[1][a]=Sxi(i+1,j,k,a);
                s[2][a]=-Set(i,j,k,a);  s[3][a]=Set(i,j+1,k,a);
                s[4][a]=-Sze(i,j,k,a);  s[5][a]=Sze(i,j,k+1,a);
            }
            auto avg4=[](double a,double b0,double c,double d){return 0.25*(a+b0+c+d);};
            auto avg8=[](double a0,double a1,double a2,double a3,double a4,double a5,double a6,double a7)
            { return 0.125*(a0+a1+a2+a3+a4+a5+a6+a7); };
            const double xc=avg8(b.x(i,j,k),b.x(i+1,j,k),b.x(i,j+1,k),b.x(i+1,j+1,k),b.x(i,j,k+1),b.x(i+1,j,k+1),b.x(i,j+1,k+1),b.x(i+1,j+1,k+1));
            const double yc=avg8(b.y(i,j,k),b.y(i+1,j,k),b.y(i,j+1,k),b.y(i+1,j+1,k),b.y(i,j,k+1),b.y(i+1,j,k+1),b.y(i,j+1,k+1),b.y(i+1,j+1,k+1));
            const double zc=avg8(b.z(i,j,k),b.z(i+1,j,k),b.z(i,j+1,k),b.z(i+1,j+1,k),b.z(i,j,k+1),b.z(i+1,j,k+1),b.z(i,j+1,k+1),b.z(i+1,j+1,k+1));
            double fc[6][3];
            auto setfc=[&](int n,int i0,int j0,int k0,int i1,int j1,int k1,int i2,int j2,int k2,int i3,int j3,int k3)
            {
                fc[n][0]=avg4(b.x(i0,j0,k0),b.x(i1,j1,k1),b.x(i2,j2,k2),b.x(i3,j3,k3));
                fc[n][1]=avg4(b.y(i0,j0,k0),b.y(i1,j1,k1),b.y(i2,j2,k2),b.y(i3,j3,k3));
                fc[n][2]=avg4(b.z(i0,j0,k0),b.z(i1,j1,k1),b.z(i2,j2,k2),b.z(i3,j3,k3));
            };
            setfc(0,i,j,k,i,j+1,k,i,j,k+1,i,j+1,k+1);
            setfc(1,i+1,j,k,i+1,j+1,k,i+1,j,k+1,i+1,j+1,k+1);
            setfc(2,i,j,k,i+1,j,k,i,j,k+1,i+1,j,k+1);
            setfc(3,i,j+1,k,i+1,j+1,k,i,j+1,k+1,i+1,j+1,k+1);
            setfc(4,i,j,k,i+1,j,k,i,j+1,k,i+1,j+1,k);
            setfc(5,i,j,k+1,i+1,j,k+1,i,j+1,k+1,i+1,j+1,k+1);
            double smag[6],mu[6],sref=0.0;
            for(int n=0;n<6;++n)
            {
                smag[n]=std::sqrt(s[n][0]*s[n][0]+s[n][1]*s[n][1]+s[n][2]*s[n][2]);
                sref=std::max(sref,smag[n]);
                const double dx=fc[n][0]-xc,dy=fc[n][1]-yc,dz=fc[n][2]-zc;
                mu[n]=1.0/std::max(std::sqrt(dx*dx+dy*dy+dz*dz),1.0e-14);
            }
            double M[3][3]={{0,0,0},{0,0,0},{0,0,0}};
            bool active[6];
            for(int n=0;n<6;++n)
            {
                active[n]=smag[n]>1.0e-10*sref;
                if(!active[n]) continue;
                for(int a=0;a<3;++a) for(int c=0;c<3;++c) M[a][c]+=mu[n]*s[n][a]*s[n][c];
            }
            double A[3][3]; for(int a=0;a<3;++a) for(int c=0;c<3;++c) A[a][c]=M[a][c];
            double V[3][3],lam[3]; eigen_sym3(A,V,lam);
            const double lmaxeig=std::max(lam[0],std::max(lam[1],lam[2]));
            for(int m=0;m<18;++m) W(i,j,k,m)=0.0;
            if(lmaxeig<=eps) continue;
            double Minv[3][3]={{0,0,0},{0,0,0},{0,0,0}};
            for(int q=0;q<3;++q)
            {
                const double inv=(lam[q]>eig_tol*lmaxeig)?1.0/lam[q]:0.0;
                for(int a=0;a<3;++a) for(int c=0;c<3;++c) Minv[a][c]+=V[a][q]*inv*V[c][q];
            }
            for(int n=0;n<6;++n) if(active[n]) for(int a=0;a<3;++a)
                W(i,j,k,6*a+n)=mu[n]*(Minv[a][0]*s[n][0]+Minv[a][1]*s[n][1]+Minv[a][2]*s[n][2]);
        }
    }

    // Edge 1-form line integrals -> Cartesian Jcell (twelve cell edges).
    for (int ib=0;ib<fld_->num_blocks();++ib)
    {
        auto &Jc=fld_->field(fid_.fid_Jcell,ib);
        auto &W=fld_->field(fid_.fid_Jcell_from_Jedge_w,ib);
        if(!Jc.is_allocated()||!W.is_allocated()) continue;
        auto &b=grd_->grids(ib);
        const Int3 lo=Jc.inner_lo(), hi=Jc.inner_hi();
        for(int i=lo.i;i<hi.i;++i) for(int j=lo.j;j<hi.j;++j) for(int k=lo.k;k<hi.k;++k)
        {
            double l[12][3];
            auto edge=[&](int n,int i0,int j0,int k0,int i1,int j1,int k1)
            {
                l[n][0]=b.x(i1,j1,k1)-b.x(i0,j0,k0);
                l[n][1]=b.y(i1,j1,k1)-b.y(i0,j0,k0);
                l[n][2]=b.z(i1,j1,k1)-b.z(i0,j0,k0);
            };
            edge(0,i,j,k,i+1,j,k); edge(1,i,j+1,k,i+1,j+1,k);
            edge(2,i,j,k+1,i+1,j,k+1); edge(3,i,j+1,k+1,i+1,j+1,k+1);
            edge(4,i,j,k,i,j+1,k); edge(5,i+1,j,k,i+1,j+1,k);
            edge(6,i,j,k+1,i,j+1,k+1); edge(7,i+1,j,k+1,i+1,j+1,k+1);
            edge(8,i,j,k,i,j,k+1); edge(9,i+1,j,k,i+1,j,k+1);
            edge(10,i,j+1,k,i,j+1,k+1); edge(11,i+1,j+1,k,i+1,j+1,k+1);

            double M[3][3]={{0,0,0},{0,0,0},{0,0,0}};
            double lmax=0.0;
            for(int n=0;n<12;++n)
            {
                const double ln=std::sqrt(l[n][0]*l[n][0]+l[n][1]*l[n][1]+l[n][2]*l[n][2]);
                lmax=std::max(lmax,ln);
            }
            for(int m=0;m<36;++m) W(i,j,k,m)=0.0;
            if(lmax<=eps) continue;
            bool active[12];
            for(int n=0;n<12;++n)
            {
                const double ln=std::sqrt(l[n][0]*l[n][0]+l[n][1]*l[n][1]+l[n][2]*l[n][2]);
                active[n]=ln>1.0e-10*lmax;
                if(!active[n]) continue;
                for(int a=0;a<3;++a) for(int c=0;c<3;++c) M[a][c]+=l[n][a]*l[n][c];
            }
            double A[3][3]; for(int a=0;a<3;++a) for(int c=0;c<3;++c) A[a][c]=M[a][c];
            double V[3][3],lam[3]; eigen_sym3(A,V,lam);
            const double lmaxeig=std::max(lam[0],std::max(lam[1],lam[2]));
            if(lmaxeig<=eps) continue;
            double Minv[3][3]={{0,0,0},{0,0,0},{0,0,0}};
            for(int q=0;q<3;++q)
            {
                const double inv=(lam[q]>eig_tol*lmaxeig)?1.0/lam[q]:0.0;
                for(int a=0;a<3;++a) for(int c=0;c<3;++c) Minv[a][c]+=V[a][q]*inv*V[c][q];
            }
            for(int n=0;n<12;++n) if(active[n])
                for(int a=0;a<3;++a)
                    W(i,j,k,12*a+n)=Minv[a][0]*l[n][0]+Minv[a][1]*l[n][1]+Minv[a][2]*l[n][2];
        }
    }

    // ------------------------------------------------------------------
    // Consistent local face Hodge M2 (MFD construction for flux DOFs).
    //
    // S contains the six outward face-area vectors and X the vectors from
    // cell center to face center.  For an exact polyhedron S^T X = V I and
    // the consistency condition for face-flux cochains is
    //
    //                       M2 S = X.
    //
    // X (X^T S)^-1 X^T supplies the consistent rank-three part.  A positive
    // projector on the complement of range(S) supplies stabilization, making
    // the complete 6x6 local matrix symmetric positive definite without a
    // global solve.  The global application is assembled face-wise later.
    // ------------------------------------------------------------------
    double m2_consistency_local=0.0;
    double m2_symmetry_local=0.0;
    for(int ib=0;ib<fld_->num_blocks();++ib)
    {
        auto &M2=fld_->field(fid_.fid_Hodge_M2_cell_consistent,ib);
        auto &W=fld_->field(fid_.fid_Bcell_from_Bface_w,ib);
        auto &Jac=fld_->field(fid_.fid_Jac,ib);
        auto &Sxi=fld_->field(fid_.fid_metric.xi,ib);
        auto &Set=fld_->field(fid_.fid_metric.eta,ib);
        auto &Sze=fld_->field(fid_.fid_metric.zeta,ib);
        auto &b=grd_->grids(ib);
        // Build the same local geometry operator in ghosts.  Only the first
        // physical ghost layer is used later; seam ghosts are never assembled
        // as cells and therefore cannot double-count a neighbouring block.
        const Int3 lo=M2.get_lo(),hi=M2.get_hi();
        for(int i=lo.i;i<hi.i;++i) for(int j=lo.j;j<hi.j;++j) for(int k=lo.k;k<hi.k;++k)
        {
            double S[6][3];
            for(int a=0;a<3;++a)
            {
                S[0][a]=-Sxi(i,j,k,a); S[1][a]=Sxi(i+1,j,k,a);
                S[2][a]=-Set(i,j,k,a); S[3][a]=Set(i,j+1,k,a);
                S[4][a]=-Sze(i,j,k,a); S[5][a]=Sze(i,j,k+1,a);
            }
            auto avg4=[](double a,double b0,double c,double d){return 0.25*(a+b0+c+d);};
            auto avg8=[](double a0,double a1,double a2,double a3,double a4,double a5,double a6,double a7)
            { return 0.125*(a0+a1+a2+a3+a4+a5+a6+a7); };
            const double xc[3]={
                avg8(b.x(i,j,k),b.x(i+1,j,k),b.x(i,j+1,k),b.x(i+1,j+1,k),b.x(i,j,k+1),b.x(i+1,j,k+1),b.x(i,j+1,k+1),b.x(i+1,j+1,k+1)),
                avg8(b.y(i,j,k),b.y(i+1,j,k),b.y(i,j+1,k),b.y(i+1,j+1,k),b.y(i,j,k+1),b.y(i+1,j,k+1),b.y(i,j+1,k+1),b.y(i+1,j+1,k+1)),
                avg8(b.z(i,j,k),b.z(i+1,j,k),b.z(i,j+1,k),b.z(i+1,j+1,k),b.z(i,j,k+1),b.z(i+1,j,k+1),b.z(i,j+1,k+1),b.z(i+1,j+1,k+1))};
            double fc[6][3];
            auto setfc=[&](int n,int i0,int j0,int k0,int i1,int j1,int k1,int i2,int j2,int k2,int i3,int j3,int k3)
            {
                fc[n][0]=avg4(b.x(i0,j0,k0),b.x(i1,j1,k1),b.x(i2,j2,k2),b.x(i3,j3,k3));
                fc[n][1]=avg4(b.y(i0,j0,k0),b.y(i1,j1,k1),b.y(i2,j2,k2),b.y(i3,j3,k3));
                fc[n][2]=avg4(b.z(i0,j0,k0),b.z(i1,j1,k1),b.z(i2,j2,k2),b.z(i3,j3,k3));
            };
            setfc(0,i,j,k,i,j+1,k,i,j,k+1,i,j+1,k+1);
            setfc(1,i+1,j,k,i+1,j+1,k,i+1,j,k+1,i+1,j+1,k+1);
            setfc(2,i,j,k,i+1,j,k,i,j,k+1,i+1,j,k+1);
            setfc(3,i,j+1,k,i+1,j+1,k,i,j+1,k+1,i+1,j+1,k+1);
            setfc(4,i,j,k,i+1,j,k,i,j+1,k,i+1,j+1,k);
            setfc(5,i,j,k+1,i+1,j,k+1,i,j+1,k+1,i+1,j+1,k+1);
            double X[6][3];
            for(int n=0;n<6;++n) for(int a=0;a<3;++a) X[n][a]=fc[n][a]-xc[a];

            auto sym_pinv3=[&](double G[3][3],double Gi[3][3])->bool
            {
                double A[3][3]; for(int a=0;a<3;++a) for(int c=0;c<3;++c) A[a][c]=G[a][c];
                double V[3][3],lam[3]; eigen_sym3(A,V,lam);
                const double lmax=std::max(lam[0],std::max(lam[1],lam[2]));
                for(int a=0;a<3;++a) for(int c=0;c<3;++c) Gi[a][c]=0.0;
                if(!(lmax>eps)) return false;
                bool full=true;
                for(int q=0;q<3;++q)
                {
                    if(!(lam[q]>eig_tol*lmax)){full=false;continue;}
                    for(int a=0;a<3;++a) for(int c=0;c<3;++c)
                        Gi[a][c]+=V[a][q]*V[c][q]/lam[q];
                }
                return full;
            };

            double STS[3][3]={{0,0,0},{0,0,0},{0,0,0}};
            double STX[3][3]={{0,0,0},{0,0,0},{0,0,0}};
            for(int a=0;a<3;++a) for(int c=0;c<3;++c) for(int n=0;n<6;++n)
            {
                STS[a][c]+=S[n][a]*S[n][c];
                STX[a][c]+=S[n][a]*X[n][c];
            }
            double STSi[3][3];
            const bool goodS=sym_pinv3(STS,STSi);

            // The arithmetic center of a warped quadrilateral is not its
            // exact area centroid, so raw geometry generally misses the
            // polyhedral identity S^T X = V I by O(h^2).  Apply the minimum
            // Frobenius-norm correction in range(S); this is local geometric
            // closure and changes neither mesh nodes nor topology.
            if(goodS)
            {
                const double vol=std::abs(Jac(i,j,k,0));
                double R[3][3];
                for(int a=0;a<3;++a) for(int c=0;c<3;++c)
                    R[a][c]=((a==c)?vol:0.0)-STX[a][c];
                for(int n=0;n<6;++n) for(int c=0;c<3;++c)
                {
                    double dx=0.0;
                    for(int a=0;a<3;++a) for(int d=0;d<3;++d)
                        dx+=S[n][a]*STSi[a][d]*R[d][c];
                    X[n][c]+=dx;
                }
            }

            double G[3][3]={{0,0,0},{0,0,0},{0,0,0}};
            for(int a=0;a<3;++a) for(int c=0;c<3;++c) for(int n=0;n<6;++n)
                G[a][c]+=0.5*(X[n][a]*S[n][c]+S[n][a]*X[n][c]);
            double Gi[3][3];
            const bool goodG=sym_pinv3(G,Gi);
            double M[6][6]={{0}};
            if(goodG)
                for(int n=0;n<6;++n) for(int m=0;m<6;++m)
                    for(int a=0;a<3;++a) for(int c=0;c<3;++c)
                        M[n][m]+=X[n][a]*Gi[a][c]*X[m][c];
            else
            {
                const double vol=std::abs(Jac(i,j,k,0));
                for(int n=0;n<6;++n) for(int m=0;m<6;++m)
                    for(int a=0;a<3;++a)
                        M[n][m]+=vol*W(i,j,k,6*a+n)*W(i,j,k,6*a+m);
            }

            double P[6][6]={{0}};
            for(int n=0;n<6;++n) for(int m=0;m<6;++m)
            {
                P[n][m]=(n==m)?1.0:0.0;
                if(goodS)
                    for(int a=0;a<3;++a) for(int c=0;c<3;++c)
                        P[n][m]-=S[n][a]*STSi[a][c]*S[m][c];
            }
            double gamma=0.0;
            for(int n=0;n<6;++n) gamma+=M[n][n];
            gamma=std::max(gamma/3.0,1.0e-14/std::max(std::cbrt(std::abs(Jac(i,j,k,0))),1.0e-14));
            for(int n=0;n<6;++n) for(int m=0;m<6;++m)
            {
                M[n][m]+=gamma*P[n][m];
                M2(i,j,k,6*n+m)=M[n][m];
            }
            const Int3 rlo=Jac.inner_lo(),rhi=Jac.inner_hi();
            if(i>=rlo.i&&i<rhi.i&&j>=rlo.j&&j<rhi.j&&k>=rlo.k&&k<rhi.k)
                for(int n=0;n<6;++n) for(int m=0;m<6;++m)
                    m2_symmetry_local=std::max(m2_symmetry_local,
                                               std::abs(M[n][m]-M[m][n]));
            double err2=0.0,ref2=0.0;
            for(int n=0;n<6;++n) for(int a=0;a<3;++a)
            {
                double ms=0.0; for(int m=0;m<6;++m) ms+=M[n][m]*S[m][a];
                const double d=ms-X[n][a]; err2+=d*d; ref2+=X[n][a]*X[n][a];
            }
            if(i>=rlo.i&&i<rhi.i&&j>=rlo.j&&j<rhi.j&&k>=rlo.k&&k<rhi.k)
                m2_consistency_local=std::max(m2_consistency_local,
                                              std::sqrt(err2/std::max(ref2,1.0e-300)));
        }
    }
    double dlocal[2]={m2_consistency_local,m2_symmetry_local},dglobal[2]={0,0};
    PARALLEL::mpi_max(dlocal,dglobal,2);
    int myid=0; PARALLEL::mpi_rank(&myid);
    if(myid==0)
        std::cout << "[Lunar][Hodge] local_M2 consistency_rel=" << dglobal[0]
                  << " symmetry_abs=" << dglobal[1] << "\n";
}
