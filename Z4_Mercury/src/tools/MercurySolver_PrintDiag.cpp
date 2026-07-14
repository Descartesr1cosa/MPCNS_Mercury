#include <limits>
#include <algorithm>
#include <cstdio>
#include <array>
#include <cmath>
#include "MercurySolver.h"
#include "1_grid/1_MPCNS_Grid.h"

// 统计 rho/p (H, Na) 与 Bx/By/Bz 的范围，并打印（rank0）
void MercurySolver::PrintMinMaxDiagnostics_()
{
    constexpr int NR = 14;
    constexpr double R0 = 0.70;
    constexpr double DR = 0.05;
    std::array<double, NR> jsum_l{}, jmax_l{}, count_l{}, zero_l{};
    std::array<double, NR> jesum_l{}, jemax_l{}, jecount_l{}, jezero_l{};
    double bind_face_max_l = 0.0;
    // --- local init ---
    double rhoH_min_l = std::numeric_limits<double>::infinity();
    double rhoH_max_l = -std::numeric_limits<double>::infinity();
    double pH_min_l = std::numeric_limits<double>::infinity();
    double pH_max_l = -std::numeric_limits<double>::infinity();

    double rhoNa_min_l = std::numeric_limits<double>::infinity();
    double rhoNa_max_l = -std::numeric_limits<double>::infinity();
    double pNa_min_l = std::numeric_limits<double>::infinity();
    double pNa_max_l = -std::numeric_limits<double>::infinity();

    double bx_min_l = std::numeric_limits<double>::infinity();
    double bx_max_l = -std::numeric_limits<double>::infinity();
    double by_min_l = std::numeric_limits<double>::infinity();
    double by_max_l = -std::numeric_limits<double>::infinity();
    double bz_min_l = std::numeric_limits<double>::infinity();
    double bz_max_l = -std::numeric_limits<double>::infinity();
    double b2_min_l = std::numeric_limits<double>::infinity();
    double b2_max_l = -std::numeric_limits<double>::infinity();
    double divb_max_l = -std::numeric_limits<double>::infinity();

    const int nblock = fld_->num_blocks();
    // Per-block shell statistics expose panel/coupling failures that disappear
    // in a global radial average.  The shell is the conducting Mercury layer.
    std::vector<double> shell_jsum_l(nblock, 0.0), shell_jmax_l(nblock, 0.0);
    std::vector<double> shell_count_l(nblock, 0.0), shell_jzero_l(nblock, 0.0);
    std::vector<double> shell_wzero_l(nblock, 0.0);
    constexpr int NSHELL = 16;
    std::array<double, NSHELL> shell_r_sum_l{}, shell_r_max_l{}, shell_r_n_l{};
    std::array<double, NSHELL> shell_r_jzero_l{}, shell_r_edgezero_l{};
    for (int ib = 0; ib < nblock; ++ib)
    {
        auto &UH = fld_->field(fid_.fid_U_H, ib);
        auto &UNa = fld_->field(fid_.fid_U_Na, ib);
        auto &PVH = fld_->field(fid_.fid_PV_H, ib);
        auto &PVNa = fld_->field(fid_.fid_PV_Na, ib);
        auto &Bcel = fld_->field(fid_.fid_Bcell, ib);
        auto &Jcel = fld_->field(fid_.fid_Jcell, ib);
        auto &divB = fld_->field(fid_.fid_divB, ib);

        if (!UH.is_allocated() || !UNa.is_allocated() || !PVH.is_allocated() ||
            !PVNa.is_allocated() || !Bcel.is_allocated() || !Jcel.is_allocated())
            continue;

        const Int3 lo = UH.inner_lo();
        const Int3 hi = UH.inner_hi();

        for (int k = lo.k; k < hi.k; ++k)
            for (int j = lo.j; j < hi.j; ++j)
                for (int i = lo.i; i < hi.i; ++i)
                {
                    const double rhoH = UH(i, j, k, 0);
                    const double pH = PVH(i, j, k, 3);

                    const double rhoNa = UNa(i, j, k, 0);
                    const double pNa = PVNa(i, j, k, 3);

                    const double bx = Bcel(i, j, k, 0);
                    const double by = Bcel(i, j, k, 1);
                    const double bz = Bcel(i, j, k, 2);
                    const double divb = fabs(divB(i, j, k, 0));

                    auto &blk = grd_->grids(ib);
                    const double x = blk.dual_x(i + 1, j + 1, k + 1);
                    const double y = blk.dual_y(i + 1, j + 1, k + 1);
                    const double z = blk.dual_z(i + 1, j + 1, k + 1);
                    const double radius = std::sqrt(x*x + y*y + z*z);
                    const double jx = Jcel(i, j, k, 0);
                    const double jy = Jcel(i, j, k, 1);
                    const double jz = Jcel(i, j, k, 2);
                    const double jmag = std::sqrt(jx*jx + jy*jy + jz*jz);
                    const int ir = static_cast<int>(std::floor((radius - R0) / DR));
                    if (ir >= 0 && ir < NR)
                    {
                        jsum_l[ir] += jmag;
                        jmax_l[ir] = std::max(jmax_l[ir], jmag);
                        count_l[ir] += 1.0;
                        if (jmag < 1.0e-12)
                            zero_l[ir] += 1.0;
                    }

                    rhoH_min_l = std::min(rhoH_min_l, rhoH);
                    rhoH_max_l = std::max(rhoH_max_l, rhoH);
                    pH_min_l = std::min(pH_min_l, pH);
                    pH_max_l = std::max(pH_max_l, pH);

                    rhoNa_min_l = std::min(rhoNa_min_l, rhoNa);
                    rhoNa_max_l = std::max(rhoNa_max_l, rhoNa);
                    pNa_min_l = std::min(pNa_min_l, pNa);
                    pNa_max_l = std::max(pNa_max_l, pNa);

                    bx_min_l = std::min(bx_min_l, bx);
                    bx_max_l = std::max(bx_max_l, bx);
                    by_min_l = std::min(by_min_l, by);
                    by_max_l = std::max(by_max_l, by);
                    bz_min_l = std::min(bz_min_l, bz);
                    bz_max_l = std::max(bz_max_l, bz);

                    b2_min_l = std::min(b2_min_l, 0.5 * (bx * bx + by * by + bz * bz));
                    b2_max_l = std::max(b2_max_l, 0.5 * (bx * bx + by * by + bz * bz));

                    divb_max_l = std::max(divb_max_l, divb);
                }
    }

    // The solid Mercury blocks do not allocate fluid fields, so they must not
    // be gated by UH/PVH availability as the general fluid diagnostics are.
    for (int ib = 0; ib < nblock; ++ib)
    {
        auto &Jcel = fld_->field(fid_.fid_Jcell, ib);
        auto &Wj = fld_->field(fid_.fid_Jcell_from_Jedge_w, ib);
        if (!Jcel.is_allocated()) continue;
        auto &blk = grd_->grids(ib);
        const Int3 lo=Jcel.inner_lo(), hi=Jcel.inner_hi();
        for (int i=lo.i;i<hi.i;++i) for (int j=lo.j;j<hi.j;++j) for (int k=lo.k;k<hi.k;++k)
        {
            const double x=blk.dual_x(i+1,j+1,k+1);
            const double y=blk.dual_y(i+1,j+1,k+1);
            const double z=blk.dual_z(i+1,j+1,k+1);
            const double radius=std::sqrt(x*x+y*y+z*z);
            if (radius < 0.84 || radius >= 1.0) continue;
            const double jx=Jcel(i,j,k,0), jy=Jcel(i,j,k,1), jz=Jcel(i,j,k,2);
            const double jmag=std::sqrt(jx*jx+jy*jy+jz*jz);
            const double sedge[12] = {
                fld_->field(fid_.fid_J.xi,ib)(i,j,k,0), fld_->field(fid_.fid_J.xi,ib)(i,j+1,k,0),
                fld_->field(fid_.fid_J.xi,ib)(i,j,k+1,0), fld_->field(fid_.fid_J.xi,ib)(i,j+1,k+1,0),
                fld_->field(fid_.fid_J.eta,ib)(i,j,k,0), fld_->field(fid_.fid_J.eta,ib)(i+1,j,k,0),
                fld_->field(fid_.fid_J.eta,ib)(i,j,k+1,0), fld_->field(fid_.fid_J.eta,ib)(i+1,j,k+1,0),
                fld_->field(fid_.fid_J.zeta,ib)(i,j,k,0), fld_->field(fid_.fid_J.zeta,ib)(i+1,j,k,0),
                fld_->field(fid_.fid_J.zeta,ib)(i,j+1,k,0), fld_->field(fid_.fid_J.zeta,ib)(i+1,j+1,k,0)};
            double sedge_max=0.0;
            for(double q:sedge) sedge_max=std::max(sedge_max,std::abs(q));
            const int is=std::min(NSHELL-1,std::max(0,static_cast<int>((radius-0.84)/0.01)));
            shell_r_sum_l[is]+=jmag;
            shell_r_max_l[is]=std::max(shell_r_max_l[is],jmag);
            shell_r_n_l[is]+=1.0;
            if(jmag<1.0e-12) shell_r_jzero_l[is]+=1.0;
            if(sedge_max<1.0e-12) shell_r_edgezero_l[is]+=1.0;
            shell_jsum_l[ib]+=jmag;
            shell_jmax_l[ib]=std::max(shell_jmax_l[ib],jmag);
            shell_count_l[ib]+=1.0;
            if(jmag<1.0e-12) shell_jzero_l[ib]+=1.0;
            double wnorm=0.0;
            if(Wj.is_allocated()) for(int m=0;m<36;++m)
                wnorm=std::max(wnorm,std::abs(Wj(i,j,k,m)));
            if(wnorm<1.0e-30) shell_wzero_l[ib]+=1.0;
        }
    }

    for (int ib = 0; ib < nblock; ++ib)
    {
        auto accumulate_edge = [&](FieldBlock &J, FieldBlock &dl, int axis)
        {
            if (!J.is_allocated() || !dl.is_allocated())
                return;
            auto &blk = grd_->grids(ib);
            const Int3 lo = J.inner_lo();
            const Int3 hi = J.inner_hi();
            const int di = axis == 0 ? 1 : 0;
            const int dj = axis == 1 ? 1 : 0;
            const int dk = axis == 2 ? 1 : 0;
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const double x = 0.5*(blk.x(i,j,k)+blk.x(i+di,j+dj,k+dk));
                        const double y = 0.5*(blk.y(i,j,k)+blk.y(i+di,j+dj,k+dk));
                        const double z = 0.5*(blk.z(i,j,k)+blk.z(i+di,j+dj,k+dk));
                        const double radius = std::sqrt(x*x+y*y+z*z);
                        const int ir = static_cast<int>(std::floor((radius-R0)/DR));
                        if (ir < 0 || ir >= NR)
                            continue;
                        const double jmag = std::abs(J(i,j,k,0)) /
                                            std::max(std::abs(dl(i,j,k,0)), 1.0e-30);
                        jesum_l[ir] += jmag;
                        jemax_l[ir] = std::max(jemax_l[ir], jmag);
                        jecount_l[ir] += 1.0;
                        if (jmag < 1.0e-12)
                            jezero_l[ir] += 1.0;
                    }
        };
        accumulate_edge(fld_->field(fid_.fid_J.xi, ib), fld_->field(fid_.Edge_dl.xi, ib), 0);
        accumulate_edge(fld_->field(fid_.fid_J.eta, ib), fld_->field(fid_.Edge_dl.eta, ib), 1);
        accumulate_edge(fld_->field(fid_.fid_J.zeta, ib), fld_->field(fid_.Edge_dl.zeta, ib), 2);

        auto face_max = [&](FieldBlock &Bind)
        {
            if (!Bind.is_allocated()) return;
            const Int3 lo=Bind.inner_lo(), hi=Bind.inner_hi();
            for(int i=lo.i;i<hi.i;++i) for(int j=lo.j;j<hi.j;++j) for(int k=lo.k;k<hi.k;++k)
                bind_face_max_l=std::max(bind_face_max_l,std::abs(Bind(i,j,k,0)));
        };
        face_max(fld_->field(fid_.fid_B.xi,ib));
        face_max(fld_->field(fid_.fid_B.eta,ib));
        face_max(fld_->field(fid_.fid_B.zeta,ib));
    }

    // --- MPI reduction ---
    double mins_l[8] = {rhoH_min_l, pH_min_l, rhoNa_min_l, pNa_min_l, bx_min_l, by_min_l, bz_min_l, b2_min_l};
    double mins_g[8];

    double maxs_l[9] = {rhoH_max_l, pH_max_l, rhoNa_max_l, pNa_max_l, bx_max_l, by_max_l, bz_max_l, b2_max_l, divb_max_l};
    double maxs_g[9];

    PARALLEL::mpi_min(mins_l, mins_g, 8);
    PARALLEL::mpi_max(maxs_l, maxs_g, 9);
    std::array<double, NR> jsum_g{}, jmax_g{}, count_g{}, zero_g{};
    PARALLEL::mpi_sum(jsum_l.data(), jsum_g.data(), NR);
    PARALLEL::mpi_max(jmax_l.data(), jmax_g.data(), NR);
    PARALLEL::mpi_sum(count_l.data(), count_g.data(), NR);
    PARALLEL::mpi_sum(zero_l.data(), zero_g.data(), NR);
    std::array<double, NR> jesum_g{}, jemax_g{}, jecount_g{}, jezero_g{};
    PARALLEL::mpi_sum(jesum_l.data(), jesum_g.data(), NR);
    PARALLEL::mpi_max(jemax_l.data(), jemax_g.data(), NR);
    PARALLEL::mpi_sum(jecount_l.data(), jecount_g.data(), NR);
    PARALLEL::mpi_sum(jezero_l.data(), jezero_g.data(), NR);
    double bind_face_max_g=0.0;
    PARALLEL::mpi_max(&bind_face_max_l,&bind_face_max_g,1);
    std::array<double,NSHELL> shell_r_sum_g{}, shell_r_max_g{}, shell_r_n_g{};
    std::array<double,NSHELL> shell_r_jzero_g{}, shell_r_edgezero_g{};
    PARALLEL::mpi_sum(shell_r_sum_l.data(),shell_r_sum_g.data(),NSHELL);
    PARALLEL::mpi_max(shell_r_max_l.data(),shell_r_max_g.data(),NSHELL);
    PARALLEL::mpi_sum(shell_r_n_l.data(),shell_r_n_g.data(),NSHELL);
    PARALLEL::mpi_sum(shell_r_jzero_l.data(),shell_r_jzero_g.data(),NSHELL);
    PARALLEL::mpi_sum(shell_r_edgezero_l.data(),shell_r_edgezero_g.data(),NSHELL);

    const int myid = par_->GetInt("myid");
    if (myid == 0)
    {
        // -----------------------------
        // NEW: compute |B|_max from Bmag_max (=0.5|B|^2)
        // -----------------------------
        const double Bmag_max = maxs_g[7]; // 0.5|B|^2
        const double Babs_max = std::sqrt(std::max(0.0, 2.0 * Bmag_max));

        // -----------------------------
        // NEW: Omega0_s and dt*Omega_max
        // Omega0_s = (q_e * B_ref * L_ref) / (m_s * U_ref)
        // (dimensionless w.r.t. t' = t * U_ref/L_ref)
        // -----------------------------
        const auto &C = par_->GetDou_List("constant").data;
        const auto &R = par_->GetDou_List("REF").data;

        const double NA = C.at("NA");
        const double q_e = C.at("q_e");
        const double Lref = R.at("L_ref");
        const double Uref = R.at("U");
        const double Bref = R.at("B_ref");

        const double mH = par_->GetDou("mole_mass1") / NA;  // kg/particle
        const double mNa = par_->GetDou("mole_mass2") / NA; // kg/particle

        const double Omega0_H = (q_e * Bref * Lref) / (mH * Uref);
        const double Omega0_Na = (q_e * Bref * Lref) / (mNa * Uref);

        const double OmegaH_max = Omega0_H * Babs_max;
        const double OmegaNa_max = Omega0_Na * Babs_max;
        const double Omega_max = std::max(OmegaH_max, OmegaNa_max);

        const double dtOmega = dt * Omega_max;

        // 与输出对齐："[  Diag  ] " 长度是 11，所以这里用 11 个空格缩进
        // std::printf("           rhoH=[%.3e, %.3e]  pH=[%.3e, %.3e]  rhoNa=[%.3e, %.3e]  pNa=[%.3e, %.3e]\n",
        //             mins_g[0], maxs_g[0], mins_g[1], maxs_g[1], mins_g[2], maxs_g[2], mins_g[3], maxs_g[3]);
        std::printf("\n           rhoH=[%.3e, %.3e]  pH=[%.3e, %.3e]\n",
                    mins_g[0], maxs_g[0], mins_g[1], maxs_g[1]);
        std::printf("           rhoN=[%.3e, %.3e]  pN=[%.3e, %.3e]\n",
                    mins_g[2], maxs_g[2], mins_g[3], maxs_g[3]);

        std::printf("           Bx  =[%.3e, %.3e]  By=[%.3e, %.3e]  Bz   =[%.3e, %.3e]  \n",
                    mins_g[4], maxs_g[4], mins_g[5], maxs_g[5], mins_g[6], maxs_g[6]);
        std::printf("           Pmag=[%.3e, %.3e]  |B|max=%.3e         divB_max =%.3e\n\n",
                    mins_g[7], maxs_g[7], Babs_max, maxs_g[8]);

        std::printf("           radial Jmag:       r-range       mean          max       zero-frac\n");
        for (int ir = 0; ir < NR; ++ir)
        {
            if (count_g[ir] <= 0.0)
                continue;
            const double rlo = R0 + ir * DR;
            std::printf("                            [%.2f,%.2f)  %.3e   %.3e   %.3e\n",
                        rlo, rlo + DR, jsum_g[ir] / count_g[ir], jmax_g[ir],
                        zero_g[ir] / count_g[ir]);
        }
        std::printf("\n");
        std::printf("           radial |Jedge|/|dl|: r-range       mean          max       zero-frac\n");
        for (int ir=0;ir<NR;++ir)
        {
            if (jecount_g[ir] <= 0.0) continue;
            const double rlo=R0+ir*DR;
            std::printf("                            [%.2f,%.2f)  %.3e   %.3e   %.3e\n",
                        rlo,rlo+DR,jesum_g[ir]/jecount_g[ir],jemax_g[ir],
                        jezero_g[ir]/jecount_g[ir]);
        }
        std::printf("           max |Bind_face flux| = %.6e\n\n",bind_face_max_g);
        std::printf("           solid Jcell detail, 0.84<=r<1:\n");
        for(int is=0;is<NSHELL;++is)
        {
            if(shell_r_n_g[is]<=0.0) continue;
            const double rlo=0.84+0.01*is;
            std::printf("             [%.2f,%.2f): mean=%.3e max=%.3e Jzero=%.3e all12edgezero=%.3e\n",
                        rlo,rlo+0.01,shell_r_sum_g[is]/shell_r_n_g[is],shell_r_max_g[is],
                        shell_r_jzero_g[is]/shell_r_n_g[is],shell_r_edgezero_g[is]/shell_r_n_g[is]);
        }
        std::printf("\n");

        // // -----------------------------
        // // NEW: stiffness diagnostics
        // // -----------------------------
        // std::printf("           |B|max=%.3e  Omega0_H=%.3e  Omega0_Na=%.3e\n",
        //             Babs_max, Omega0_H, Omega0_Na);
        // std::printf("           OmegaH_max=%.3e  OmegaNa_max=%.3e  dt*Omega_max=%.3e\n",
        //             OmegaH_max, OmegaNa_max, dtOmega);

        std::fflush(stdout);
    }

    // Local block counts differ between MPI ranks, so print them in rank order
    // instead of reducing variable-length arrays.
    int nrank = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
    for (int rank = 0; rank < nrank; ++rank)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        if (myid != rank) continue;
        for (int ib=0; ib<nblock; ++ib)
        {
            if (shell_count_l[ib] <= 0.0) continue;
            std::printf("[ShellBlock] rank=%d local_block=%d meanJ=%.6e maxJ=%.6e zero=%.6e Wzero=%.6e n=%.0f\n",
                        myid, ib, shell_jsum_l[ib]/shell_count_l[ib], shell_jmax_l[ib],
                        shell_jzero_l[ib]/shell_count_l[ib],
                        shell_wzero_l[ib]/shell_count_l[ib], shell_count_l[ib]);
        }
        std::fflush(stdout);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}
