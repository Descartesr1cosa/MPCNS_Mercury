#include "LunarSolver.h"

#include <mpi.h>
#include <set>
#include <tuple>
#include <vector>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>
void LunarSolver::DebugFindExtremaInner(const std::vector<int> &fids,
                                          const std::vector<std::string> &names,
                                          bool print_min,
                                          bool print_max)
{
    static_assert(std::is_trivially_copyable<DebugItem>::value,
                  "DebugItem must be trivially copyable for MPI byte transfer.");

    const int myid = par_->GetInt("myid");
    int nrank = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);

    auto resolve_fields = [&]() -> std::vector<std::pair<int, std::string>>
    {
        std::vector<std::pair<int, std::string>> out;
        std::set<int> seen;

        for (const auto &nm : names)
        {
            if (!fld_->has_field(nm))
            {
                if (myid == 0)
                    std::cout << "[DebugFindExtremaInner] skip unknown field name: " << nm << "\n";
                continue;
            }
            const int fid = fld_->field_id(nm);
            if (seen.insert(fid).second)
                out.push_back({fid, nm});
        }

        for (int fid : fids)
        {
            if (fid < 0)
                continue;
            if (seen.insert(fid).second)
            {
                std::ostringstream oss;
                oss << "fid=" << fid;
                out.push_back({fid, oss.str()});
            }
        }
        return out;
    };

    auto get_xyz = [&](int fid, int ib, int i, int j, int k,
                       double &x, double &y, double &z)
    {
        auto avg4 = [](double a, double b, double c, double d) -> double
        {
            return 0.25 * (a + b + c + d);
        };
        auto avg8 = [](double a0, double a1, double a2, double a3,
                       double a4, double a5, double a6, double a7) -> double
        {
            return 0.125 * (a0 + a1 + a2 + a3 + a4 + a5 + a6 + a7);
        };

        auto &gx = grd_->grids(ib).x;
        auto &gy = grd_->grids(ib).y;
        auto &gz = grd_->grids(ib).z;

        const auto &desc = fld_->descriptor(fid);

        switch (desc.location)
        {
        case StaggerLocation::Cell:
            x = avg8(gx(i, j, k), gx(i + 1, j, k), gx(i, j + 1, k), gx(i + 1, j + 1, k),
                     gx(i, j, k + 1), gx(i + 1, j, k + 1), gx(i, j + 1, k + 1), gx(i + 1, j + 1, k + 1));
            y = avg8(gy(i, j, k), gy(i + 1, j, k), gy(i, j + 1, k), gy(i + 1, j + 1, k),
                     gy(i, j, k + 1), gy(i + 1, j, k + 1), gy(i, j + 1, k + 1), gy(i + 1, j + 1, k + 1));
            z = avg8(gz(i, j, k), gz(i + 1, j, k), gz(i, j + 1, k), gz(i + 1, j + 1, k),
                     gz(i, j, k + 1), gz(i + 1, j, k + 1), gz(i, j + 1, k + 1), gz(i + 1, j + 1, k + 1));
            break;

        case StaggerLocation::Node:
            x = gx(i, j, k);
            y = gy(i, j, k);
            z = gz(i, j, k);
            break;

        case StaggerLocation::FaceXi:
            x = avg4(gx(i, j, k), gx(i, j + 1, k), gx(i, j, k + 1), gx(i, j + 1, k + 1));
            y = avg4(gy(i, j, k), gy(i, j + 1, k), gy(i, j, k + 1), gy(i, j + 1, k + 1));
            z = avg4(gz(i, j, k), gz(i, j + 1, k), gz(i, j, k + 1), gz(i, j + 1, k + 1));
            break;

        case StaggerLocation::FaceEt:
            x = avg4(gx(i, j, k), gx(i + 1, j, k), gx(i, j, k + 1), gx(i + 1, j, k + 1));
            y = avg4(gy(i, j, k), gy(i + 1, j, k), gy(i, j, k + 1), gy(i + 1, j, k + 1));
            z = avg4(gz(i, j, k), gz(i + 1, j, k), gz(i, j, k + 1), gz(i + 1, j, k + 1));
            break;

        case StaggerLocation::FaceZe:
            x = avg4(gx(i, j, k), gx(i + 1, j, k), gx(i, j + 1, k), gx(i + 1, j + 1, k));
            y = avg4(gy(i, j, k), gy(i + 1, j, k), gy(i, j + 1, k), gy(i + 1, j + 1, k));
            z = avg4(gz(i, j, k), gz(i + 1, j, k), gz(i, j + 1, k), gz(i + 1, j + 1, k));
            break;

        case StaggerLocation::EdgeXi:
            x = 0.5 * (gx(i, j, k) + gx(i + 1, j, k));
            y = 0.5 * (gy(i, j, k) + gy(i + 1, j, k));
            z = 0.5 * (gz(i, j, k) + gz(i + 1, j, k));
            break;

        case StaggerLocation::EdgeEt:
            x = 0.5 * (gx(i, j, k) + gx(i, j + 1, k));
            y = 0.5 * (gy(i, j, k) + gy(i, j + 1, k));
            z = 0.5 * (gz(i, j, k) + gz(i, j + 1, k));
            break;

        case StaggerLocation::EdgeZe:
            x = 0.5 * (gx(i, j, k) + gx(i, j, k + 1));
            y = 0.5 * (gy(i, j, k) + gy(i, j, k + 1));
            z = 0.5 * (gz(i, j, k) + gz(i, j, k + 1));
            break;

        default:
            x = y = z = 0.0;
            break;
        }
    };

    auto fields = resolve_fields();
    if (fields.empty())
    {
        if (myid == 0)
            std::cout << "[DebugFindExtremaInner] no valid fields selected.\n";
        MPI_Barrier(MPI_COMM_WORLD);
        return;
    }

    const int nblock = fld_->num_blocks();

    for (const auto &it_field : fields)
    {
        const int fid = it_field.first;
        const std::string label = it_field.second;

        const auto &desc = fld_->descriptor(fid);
        const int ncomp = desc.ncomp;

        for (int m = 0; m < ncomp; ++m)
        {
            DebugItem loc_min, loc_max;
            loc_min.val = std::numeric_limits<double>::max();
            loc_max.val = -std::numeric_limits<double>::max();
            loc_min.fid = fid;
            loc_max.fid = fid;
            loc_min.comp = m;
            loc_max.comp = m;
            loc_min.rank = myid;
            loc_max.rank = myid;
            loc_min.aux[0] = 0; // valid flag
            loc_max.aux[0] = 0;

            for (int ib = 0; ib < nblock; ++ib)
            {
                auto &F = fld_->field(fid, ib);
                if (!F.is_allocated())
                    continue;

                Int3 lo = F.inner_lo();
                Int3 hi = F.inner_hi();

                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                        {
                            const double v = F(i, j, k, m);

                            if (!loc_min.aux[0] || v < loc_min.val)
                            {
                                loc_min.aux[0] = 1;
                                loc_min.val = v;
                                loc_min.blk = ib;
                                loc_min.i = i;
                                loc_min.j = j;
                                loc_min.k = k;
                                get_xyz(fid, ib, i, j, k,
                                        loc_min.xyz[0], loc_min.xyz[1], loc_min.xyz[2]);
                            }

                            if (!loc_max.aux[0] || v > loc_max.val)
                            {
                                loc_max.aux[0] = 1;
                                loc_max.val = v;
                                loc_max.blk = ib;
                                loc_max.i = i;
                                loc_max.j = j;
                                loc_max.k = k;
                                get_xyz(fid, ib, i, j, k,
                                        loc_max.xyz[0], loc_max.xyz[1], loc_max.xyz[2]);
                            }
                        }
            }

            std::vector<DebugItem> all_min(nrank), all_max(nrank);
            MPI_Allgather(&loc_min, sizeof(DebugItem), MPI_BYTE,
                          all_min.data(), sizeof(DebugItem), MPI_BYTE, MPI_COMM_WORLD);
            MPI_Allgather(&loc_max, sizeof(DebugItem), MPI_BYTE,
                          all_max.data(), sizeof(DebugItem), MPI_BYTE, MPI_COMM_WORLD);

            DebugItem gmin, gmax;
            gmin.aux[0] = 0;
            gmax.aux[0] = 0;

            for (const auto &q : all_min)
            {
                if (!q.aux[0])
                    continue;
                if (!gmin.aux[0] || q.val < gmin.val)
                {
                    gmin = q;
                }
            }

            for (const auto &q : all_max)
            {
                if (!q.aux[0])
                    continue;
                if (!gmax.aux[0] || q.val > gmax.val)
                {
                    gmax = q;
                }
            }

            if (myid == 0)
            {
                std::ostringstream oss;
                oss << std::scientific << std::setprecision(16);

                if (print_min && gmin.aux[0])
                {
                    oss << "\n[DebugExtrema][MIN] field=" << label
                        << " comp=" << m
                        << " val=" << gmin.val
                        << " owner_rank=" << gmin.rank
                        << " block=" << gmin.blk
                        << " idx=(" << gmin.i << "," << gmin.j << "," << gmin.k << ")"
                        << " xyz=(" << gmin.xyz[0] << "," << gmin.xyz[1] << "," << gmin.xyz[2] << ")";
                }

                if (print_max && gmax.aux[0])
                {
                    oss << "\n[DebugExtrema][MAX] field=" << label
                        << " comp=" << m
                        << " val=" << gmax.val
                        << " owner_rank=" << gmax.rank
                        << " block=" << gmax.blk
                        << " idx=(" << gmax.i << "," << gmax.j << "," << gmax.k << ")"
                        << " xyz=(" << gmax.xyz[0] << "," << gmax.xyz[1] << "," << gmax.xyz[2] << ")";
                }

                oss << "\n";
                std::cout << oss.str() << std::flush;
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

void LunarSolver::DebugDumpPointFields(int query_rank,
                                         int blk, int i, int j, int k,
                                         const std::vector<int> &fids,
                                         const std::vector<std::string> &names)
{
    const int myid = par_->GetInt("myid");

    auto resolve_fields = [&]() -> std::vector<std::pair<int, std::string>>
    {
        std::vector<std::pair<int, std::string>> out;
        std::set<int> seen;

        for (const auto &nm : names)
        {
            if (!fld_->has_field(nm))
            {
                if (myid == query_rank)
                    std::cout << "[DebugDumpPointFields] skip unknown field name: " << nm << "\n";
                continue;
            }
            const int fid = fld_->field_id(nm);
            if (seen.insert(fid).second)
                out.push_back({fid, nm});
        }

        for (int fid : fids)
        {
            if (fid < 0)
                continue;
            if (seen.insert(fid).second)
            {
                std::ostringstream oss;
                oss << "fid=" << fid;
                out.push_back({fid, oss.str()});
            }
        }

        return out;
    };

    auto loc_name = [](StaggerLocation loc) -> const char *
    {
        switch (loc)
        {
        case StaggerLocation::Cell:
            return "Cell";
        case StaggerLocation::Node:
            return "Node";
        case StaggerLocation::FaceXi:
            return "FaceXi";
        case StaggerLocation::FaceEt:
            return "FaceEt";
        case StaggerLocation::FaceZe:
            return "FaceZe";
        case StaggerLocation::EdgeXi:
            return "EdgeXi";
        case StaggerLocation::EdgeEt:
            return "EdgeEt";
        case StaggerLocation::EdgeZe:
            return "EdgeZe";
        default:
            return "Unknown";
        }
    };

    auto get_xyz = [&](int fid, int ib, int ii, int jj, int kk,
                       double &x, double &y, double &z)
    {
        auto avg4 = [](double a, double b, double c, double d) -> double
        {
            return 0.25 * (a + b + c + d);
        };
        auto avg8 = [](double a0, double a1, double a2, double a3,
                       double a4, double a5, double a6, double a7) -> double
        {
            return 0.125 * (a0 + a1 + a2 + a3 + a4 + a5 + a6 + a7);
        };

        auto &gx = grd_->grids(ib).x;
        auto &gy = grd_->grids(ib).y;
        auto &gz = grd_->grids(ib).z;

        const auto &desc = fld_->descriptor(fid);

        switch (desc.location)
        {
        case StaggerLocation::Cell:
            x = avg8(gx(ii, jj, kk), gx(ii + 1, jj, kk), gx(ii, jj + 1, kk), gx(ii + 1, jj + 1, kk),
                     gx(ii, jj, kk + 1), gx(ii + 1, jj, kk + 1), gx(ii, jj + 1, kk + 1), gx(ii + 1, jj + 1, kk + 1));
            y = avg8(gy(ii, jj, kk), gy(ii + 1, jj, kk), gy(ii, jj + 1, kk), gy(ii + 1, jj + 1, kk),
                     gy(ii, jj, kk + 1), gy(ii + 1, jj, kk + 1), gy(ii, jj + 1, kk + 1), gy(ii + 1, jj + 1, kk + 1));
            z = avg8(gz(ii, jj, kk), gz(ii + 1, jj, kk), gz(ii, jj + 1, kk), gz(ii + 1, jj + 1, kk),
                     gz(ii, jj, kk + 1), gz(ii + 1, jj, kk + 1), gz(ii, jj + 1, kk + 1), gz(ii + 1, jj + 1, kk + 1));
            break;

        case StaggerLocation::Node:
            x = gx(ii, jj, kk);
            y = gy(ii, jj, kk);
            z = gz(ii, jj, kk);
            break;

        case StaggerLocation::FaceXi:
            x = avg4(gx(ii, jj, kk), gx(ii, jj + 1, kk), gx(ii, jj, kk + 1), gx(ii, jj + 1, kk + 1));
            y = avg4(gy(ii, jj, kk), gy(ii, jj + 1, kk), gy(ii, jj, kk + 1), gy(ii, jj + 1, kk + 1));
            z = avg4(gz(ii, jj, kk), gz(ii, jj + 1, kk), gz(ii, jj, kk + 1), gz(ii, jj + 1, kk + 1));
            break;

        case StaggerLocation::FaceEt:
            x = avg4(gx(ii, jj, kk), gx(ii + 1, jj, kk), gx(ii, jj, kk + 1), gx(ii + 1, jj, kk + 1));
            y = avg4(gy(ii, jj, kk), gy(ii + 1, jj, kk), gy(ii, jj, kk + 1), gy(ii + 1, jj, kk + 1));
            z = avg4(gz(ii, jj, kk), gz(ii + 1, jj, kk), gz(ii, jj, kk + 1), gz(ii + 1, jj, kk + 1));
            break;

        case StaggerLocation::FaceZe:
            x = avg4(gx(ii, jj, kk), gx(ii + 1, jj, kk), gx(ii, jj + 1, kk), gx(ii + 1, jj + 1, kk));
            y = avg4(gy(ii, jj, kk), gy(ii + 1, jj, kk), gy(ii, jj + 1, kk), gy(ii + 1, jj + 1, kk));
            z = avg4(gz(ii, jj, kk), gz(ii + 1, jj, kk), gz(ii, jj + 1, kk), gz(ii + 1, jj + 1, kk));
            break;

        case StaggerLocation::EdgeXi:
            x = 0.5 * (gx(ii, jj, kk) + gx(ii + 1, jj, kk));
            y = 0.5 * (gy(ii, jj, kk) + gy(ii + 1, jj, kk));
            z = 0.5 * (gz(ii, jj, kk) + gz(ii + 1, jj, kk));
            break;

        case StaggerLocation::EdgeEt:
            x = 0.5 * (gx(ii, jj, kk) + gx(ii, jj + 1, kk));
            y = 0.5 * (gy(ii, jj, kk) + gy(ii, jj + 1, kk));
            z = 0.5 * (gz(ii, jj, kk) + gz(ii, jj + 1, kk));
            break;

        case StaggerLocation::EdgeZe:
            x = 0.5 * (gx(ii, jj, kk) + gx(ii, jj, kk + 1));
            y = 0.5 * (gy(ii, jj, kk) + gy(ii, jj, kk + 1));
            z = 0.5 * (gz(ii, jj, kk) + gz(ii, jj, kk + 1));
            break;

        default:
            x = y = z = 0.0;
            break;
        }
    };

    MPI_Barrier(MPI_COMM_WORLD);

    if (myid == query_rank)
    {
        auto fields = resolve_fields();

        std::cout << "\n[DebugDumpPointFields] rank=" << query_rank
                  << " block=" << blk
                  << " idx=(" << i << "," << j << "," << k << ")\n";

        for (const auto &it_field : fields)
        {
            const int fid = it_field.first;
            const std::string label = it_field.second;

            if (blk < 0 || blk >= fld_->num_blocks())
            {
                std::cout << "  field=" << label << " : invalid block\n";
                continue;
            }

            auto &F = fld_->field(fid, blk);
            const auto &desc = fld_->descriptor(fid);

            if (!F.is_allocated())
            {
                std::cout << "  field=" << label << " loc=" << loc_name(desc.location)
                          << " : not allocated\n";
                continue;
            }

            Int3 lo = F.get_lo();
            Int3 hi = F.get_hi();

            if (i < lo.i || i >= hi.i ||
                j < lo.j || j >= hi.j ||
                k < lo.k || k >= hi.k)
            {
                std::cout << "  field=" << label << " loc=" << loc_name(desc.location)
                          << " : out_of_range, valid=[("
                          << lo.i << "," << lo.j << "," << lo.k << "),("
                          << hi.i << "," << hi.j << "," << hi.k << "))\n";
                continue;
            }

            double x, y, z;
            get_xyz(fid, blk, i, j, k, x, y, z);

            std::ostringstream oss;
            oss << std::scientific << std::setprecision(16);
            oss << "  field=" << label
                << " loc=" << loc_name(desc.location)
                << " xyz=(" << x << "," << y << "," << z << ")"
                << " : ";

            for (int m = 0; m < desc.ncomp; ++m)
            {
                oss << F(i, j, k, m);
                if (m + 1 < desc.ncomp)
                    oss << ", ";
            }
            oss << "\n";

            std::cout << oss.str();
        }
        std::cout.flush();
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

void LunarSolver::DebugDumpPointPartners(int query_rank,
                                           int blk, int i, int j, int k,
                                           const std::vector<int> &fids,
                                           const std::vector<std::string> &names,
                                           int ngh,
                                           bool include_topo,
                                           bool include_halo,
                                           bool include_physical)
{
    static_assert(std::is_trivially_copyable<DebugItem>::value,
                  "DebugItem must be trivially copyable for MPI byte transfer.");

    const int myid = par_->GetInt("myid");
    int nrank = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);

    auto resolve_fields = [&]() -> std::vector<std::pair<int, std::string>>
    {
        std::vector<std::pair<int, std::string>> out;
        std::set<int> seen;

        for (const auto &nm : names)
        {
            if (!fld_->has_field(nm))
            {
                if (myid == query_rank)
                    std::cout << "[DebugDumpPointPartners] skip unknown field name: " << nm << "\n";
                continue;
            }
            const int fid = fld_->field_id(nm);
            if (seen.insert(fid).second)
                out.push_back({fid, nm});
        }

        for (int fid : fids)
        {
            if (fid < 0)
                continue;
            if (seen.insert(fid).second)
            {
                std::ostringstream oss;
                oss << "fid=" << fid;
                out.push_back({fid, oss.str()});
            }
        }

        return out;
    };

    auto get_comp = [](const Int3 &x, int d) -> int
    {
        if (d == 0)
            return x.i;
        if (d == 1)
            return x.j;
        return x.k;
    };

    auto set_comp = [](Int3 &x, int d, int v)
    {
        if (d == 0)
            x.i = v;
        else if (d == 1)
            x.j = v;
        else
            x.k = v;
    };

    auto in_box = [](const Box3 &box, const Int3 &p) -> bool
    {
        return (box.lo.i <= p.i && p.i < box.hi.i) &&
               (box.lo.j <= p.j && p.j < box.hi.j) &&
               (box.lo.k <= p.k && p.k < box.hi.k);
    };

    auto loc_name = [](StaggerLocation loc) -> const char *
    {
        switch (loc)
        {
        case StaggerLocation::Cell:
            return "Cell";
        case StaggerLocation::Node:
            return "Node";
        case StaggerLocation::FaceXi:
            return "FaceXi";
        case StaggerLocation::FaceEt:
            return "FaceEt";
        case StaggerLocation::FaceZe:
            return "FaceZe";
        case StaggerLocation::EdgeXi:
            return "EdgeXi";
        case StaggerLocation::EdgeEt:
            return "EdgeEt";
        case StaggerLocation::EdgeZe:
            return "EdgeZe";
        default:
            return "Unknown";
        }
    };

    auto rel_name = [](int rel) -> const char *
    {
        switch (rel)
        {
        case 0:
            return "SELF";
        case 1:
            return "TOPO";
        case 2:
            return "HALO";
        case 3:
            return "PHYSICAL";
        default:
            return "UNKNOWN";
        }
    };

    auto patch_name = [](int tag) -> const char *
    {
        switch (tag)
        {
        case 0:
            return "SELF";
        case 1:
            return "FACE-INNER";
        case 2:
            return "FACE-PAR";
        case 3:
            return "EDGE-INNER";
        case 4:
            return "EDGE-PAR";
        case 5:
            return "VERT-INNER";
        case 6:
            return "VERT-PAR";
        case 7:
            return "PHYSICAL";
        default:
            return "UNKNOWN";
        }
    };

    auto delta_on_axis = [](StaggerLocation loc, int ax) -> int
    {
        switch (loc)
        {
        case StaggerLocation::Cell:
            return 1;
        case StaggerLocation::Node:
            return 0;
        case StaggerLocation::FaceXi:
            return (ax == 0 ? 0 : 1);
        case StaggerLocation::FaceEt:
            return (ax == 1 ? 0 : 1);
        case StaggerLocation::FaceZe:
            return (ax == 2 ? 0 : 1);
        case StaggerLocation::EdgeXi:
            return (ax == 0 ? 1 : 0);
        case StaggerLocation::EdgeEt:
            return (ax == 1 ? 1 : 0);
        case StaggerLocation::EdgeZe:
            return (ax == 2 ? 1 : 0);
        default:
            return 1;
        }
    };

    auto make_dof_box_from_node_box = [&](StaggerLocation loc, const Box3 &node_box) -> Box3
    {
        Box3 box;
        box.lo = node_box.lo;

        switch (loc)
        {
        case StaggerLocation::Cell:
            box.hi = {node_box.hi.i - 1, node_box.hi.j - 1, node_box.hi.k - 1};
            break;
        case StaggerLocation::Node:
            box.hi = {node_box.hi.i, node_box.hi.j, node_box.hi.k};
            break;
        case StaggerLocation::FaceXi:
            box.hi = {node_box.hi.i, node_box.hi.j - 1, node_box.hi.k - 1};
            break;
        case StaggerLocation::FaceEt:
            box.hi = {node_box.hi.i - 1, node_box.hi.j, node_box.hi.k - 1};
            break;
        case StaggerLocation::FaceZe:
            box.hi = {node_box.hi.i - 1, node_box.hi.j - 1, node_box.hi.k};
            break;
        case StaggerLocation::EdgeXi:
            box.hi = {node_box.hi.i - 1, node_box.hi.j, node_box.hi.k};
            break;
        case StaggerLocation::EdgeEt:
            box.hi = {node_box.hi.i, node_box.hi.j - 1, node_box.hi.k};
            break;
        case StaggerLocation::EdgeZe:
            box.hi = {node_box.hi.i, node_box.hi.j, node_box.hi.k - 1};
            break;
        default:
            box.hi = node_box.hi;
            break;
        }

        return box;
    };

    auto map_dof_index = [&](const Int3 &p_local,
                             StaggerLocation loc,
                             const TOPO::IndexTransform &tr) -> Int3
    {
        int lv[3] = {p_local.i, p_local.j, p_local.k};
        Int3 p_nb{0, 0, 0};

        for (int a = 0; a < 3; ++a)
        {
            const int b = tr.perm[a];
            const int s = tr.sign[a];
            const int off = get_comp(tr.offset, a);
            const int delta = delta_on_axis(loc, a);

            int v = s * lv[a] + off;
            if (s < 0)
                v -= delta;

            set_comp(p_nb, b, v);
        }
        return p_nb;
    };

    auto get_xyz = [&](int fid, int ib, int ii, int jj, int kk,
                       double &x, double &y, double &z)
    {
        auto avg4 = [](double a, double b, double c, double d) -> double
        {
            return 0.25 * (a + b + c + d);
        };
        auto avg8 = [](double a0, double a1, double a2, double a3,
                       double a4, double a5, double a6, double a7) -> double
        {
            return 0.125 * (a0 + a1 + a2 + a3 + a4 + a5 + a6 + a7);
        };

        auto &gx = grd_->grids(ib).x;
        auto &gy = grd_->grids(ib).y;
        auto &gz = grd_->grids(ib).z;

        const auto &desc = fld_->descriptor(fid);

        switch (desc.location)
        {
        case StaggerLocation::Cell:
            x = avg8(gx(ii, jj, kk), gx(ii + 1, jj, kk), gx(ii, jj + 1, kk), gx(ii + 1, jj + 1, kk),
                     gx(ii, jj, kk + 1), gx(ii + 1, jj, kk + 1), gx(ii, jj + 1, kk + 1), gx(ii + 1, jj + 1, kk + 1));
            y = avg8(gy(ii, jj, kk), gy(ii + 1, jj, kk), gy(ii, jj + 1, kk), gy(ii + 1, jj + 1, kk),
                     gy(ii, jj, kk + 1), gy(ii + 1, jj, kk + 1), gy(ii, jj + 1, kk + 1), gy(ii + 1, jj + 1, kk + 1));
            z = avg8(gz(ii, jj, kk), gz(ii + 1, jj, kk), gz(ii, jj + 1, kk), gz(ii + 1, jj + 1, kk),
                     gz(ii, jj, kk + 1), gz(ii + 1, jj, kk + 1), gz(ii, jj + 1, kk + 1), gz(ii + 1, jj + 1, kk + 1));
            break;

        case StaggerLocation::Node:
            x = gx(ii, jj, kk);
            y = gy(ii, jj, kk);
            z = gz(ii, jj, kk);
            break;

        case StaggerLocation::FaceXi:
            x = avg4(gx(ii, jj, kk), gx(ii, jj + 1, kk), gx(ii, jj, kk + 1), gx(ii, jj + 1, kk + 1));
            y = avg4(gy(ii, jj, kk), gy(ii, jj + 1, kk), gy(ii, jj, kk + 1), gy(ii, jj + 1, kk + 1));
            z = avg4(gz(ii, jj, kk), gz(ii, jj + 1, kk), gz(ii, jj, kk + 1), gz(ii, jj + 1, kk + 1));
            break;

        case StaggerLocation::FaceEt:
            x = avg4(gx(ii, jj, kk), gx(ii + 1, jj, kk), gx(ii, jj, kk + 1), gx(ii + 1, jj, kk + 1));
            y = avg4(gy(ii, jj, kk), gy(ii + 1, jj, kk), gy(ii, jj, kk + 1), gy(ii + 1, jj, kk + 1));
            z = avg4(gz(ii, jj, kk), gz(ii + 1, jj, kk), gz(ii, jj, kk + 1), gz(ii + 1, jj, kk + 1));
            break;

        case StaggerLocation::FaceZe:
            x = avg4(gx(ii, jj, kk), gx(ii + 1, jj, kk), gx(ii, jj + 1, kk), gx(ii + 1, jj + 1, kk));
            y = avg4(gy(ii, jj, kk), gy(ii + 1, jj, kk), gy(ii, jj + 1, kk), gy(ii + 1, jj + 1, kk));
            z = avg4(gz(ii, jj, kk), gz(ii + 1, jj, kk), gz(ii, jj + 1, kk), gz(ii + 1, jj + 1, kk));
            break;

        case StaggerLocation::EdgeXi:
            x = 0.5 * (gx(ii, jj, kk) + gx(ii + 1, jj, kk));
            y = 0.5 * (gy(ii, jj, kk) + gy(ii + 1, jj, kk));
            z = 0.5 * (gz(ii, jj, kk) + gz(ii + 1, jj, kk));
            break;

        case StaggerLocation::EdgeEt:
            x = 0.5 * (gx(ii, jj, kk) + gx(ii, jj + 1, kk));
            y = 0.5 * (gy(ii, jj, kk) + gy(ii, jj + 1, kk));
            z = 0.5 * (gz(ii, jj, kk) + gz(ii, jj + 1, kk));
            break;

        case StaggerLocation::EdgeZe:
            x = 0.5 * (gx(ii, jj, kk) + gx(ii, jj, kk + 1));
            y = 0.5 * (gy(ii, jj, kk) + gy(ii, jj, kk + 1));
            z = 0.5 * (gz(ii, jj, kk) + gz(ii, jj, kk + 1));
            break;

        default:
            x = y = z = 0.0;
            break;
        }
    };

    auto fields = resolve_fields();
    if (fields.empty())
    {
        if (myid == query_rank)
            std::cout << "[DebugDumpPointPartners] no valid fields selected.\n";
        MPI_Barrier(MPI_COMM_WORLD);
        return;
    }

    for (const auto &it_field : fields)
    {
        const int fid = it_field.first;
        const std::string label = it_field.second;
        const auto &desc = fld_->descriptor(fid);

        std::vector<DebugItem> queries;
        std::set<std::tuple<int, int, int, int, int, int, int>> uniq;

        if (myid == query_rank)
        {
            const Int3 p0{i, j, k};

            auto push_query = [&](int rank_t, int blk_t, int i_t, int j_t, int k_t,
                                  int relation, int patch_tag, int aux0 = 0)
            {
                auto key = std::make_tuple(rank_t, blk_t, i_t, j_t, k_t, relation, patch_tag);
                if (!uniq.insert(key).second)
                    return;

                DebugItem q;
                q.rank = rank_t;
                q.blk = blk_t;
                q.i = i_t;
                q.j = j_t;
                q.k = k_t;
                q.src_rank = query_rank;
                q.src_blk = blk;
                q.src_i = i;
                q.src_j = j;
                q.src_k = k;
                q.fid = fid;
                q.relation = relation;
                q.patch_tag = patch_tag;
                q.aux[0] = aux0;
                queries.push_back(q);
            };

            push_query(query_rank, blk, i, j, k, 0, 0);

            auto scan_patch_list = [&](const auto &plist, int patch_tag)
            {
                for (const auto &p : plist)
                {
                    if (p.this_rank != query_rank)
                        continue;
                    if (p.this_block != blk)
                        continue;

                    const Box3 dof_box = make_dof_box_from_node_box(desc.location, p.this_box_node);
                    const bool in_exact = in_box(dof_box, p0);

                    bool in_ext = false;
                    if (ngh > 0)
                    {
                        Box3 ext_box = dof_box;
                        Int3 nspan{p.this_box_node.hi.i - p.this_box_node.lo.i,
                                   p.this_box_node.hi.j - p.this_box_node.lo.j,
                                   p.this_box_node.hi.k - p.this_box_node.lo.k};

                        if (nspan.i == 1)
                        {
                            ext_box.lo.i -= ngh;
                            ext_box.hi.i += ngh;
                        }
                        if (nspan.j == 1)
                        {
                            ext_box.lo.j -= ngh;
                            ext_box.hi.j += ngh;
                        }
                        if (nspan.k == 1)
                        {
                            ext_box.lo.k -= ngh;
                            ext_box.hi.k += ngh;
                        }

                        in_ext = in_box(ext_box, p0);
                    }

                    if (include_topo && in_exact)
                    {
                        const Int3 p_nb = map_dof_index(p0, desc.location, p.trans);
                        push_query(p.nb_rank, p.nb_block, p_nb.i, p_nb.j, p_nb.k, 1, patch_tag);
                    }

                    if (include_halo && ngh > 0 && !in_exact && in_ext)
                    {
                        const Int3 p_nb = map_dof_index(p0, desc.location, p.trans);
                        push_query(p.nb_rank, p.nb_block, p_nb.i, p_nb.j, p_nb.k, 2, patch_tag);
                    }
                }
            };

            scan_patch_list(topo_->inner_patches, 1);
            scan_patch_list(topo_->parallel_patches, 2);
            scan_patch_list(topo_->inner_edge_patches, 3);
            scan_patch_list(topo_->parallel_edge_patches, 4);
            scan_patch_list(topo_->inner_vertex_patches, 5);
            scan_patch_list(topo_->parallel_vertex_patches, 6);

            if (include_physical)
            {
                for (const auto &p : topo_->physical_patches)
                {
                    if (p.this_rank != query_rank)
                        continue;
                    if (p.this_block != blk)
                        continue;

                    const Box3 dof_box = make_dof_box_from_node_box(desc.location, p.this_box_node);
                    const bool in_exact = in_box(dof_box, p0);

                    bool in_ext = false;
                    if (ngh > 0)
                    {
                        Box3 ext_box = dof_box;
                        Int3 nspan{p.this_box_node.hi.i - p.this_box_node.lo.i,
                                   p.this_box_node.hi.j - p.this_box_node.lo.j,
                                   p.this_box_node.hi.k - p.this_box_node.lo.k};

                        if (nspan.i == 1)
                        {
                            ext_box.lo.i -= ngh;
                            ext_box.hi.i += ngh;
                        }
                        if (nspan.j == 1)
                        {
                            ext_box.lo.j -= ngh;
                            ext_box.hi.j += ngh;
                        }
                        if (nspan.k == 1)
                        {
                            ext_box.lo.k -= ngh;
                            ext_box.hi.k += ngh;
                        }

                        in_ext = in_box(ext_box, p0);
                    }

                    if (in_exact || (ngh > 0 && in_ext))
                    {
                        std::cout << "[DebugDumpPointPartners][PHYSICAL] field=" << label
                                  << " loc=" << loc_name(desc.location)
                                  << " src(rank=" << query_rank
                                  << ", blk=" << blk
                                  << ", idx=(" << i << "," << j << "," << k << "))"
                                  << " touches bc=\"" << p.bc_name
                                  << "\" dir=" << p.direction
                                  << " ngh=" << ngh << "\n";
                    }
                }
            }
        }

        int nquery = (myid == query_rank) ? static_cast<int>(queries.size()) : 0;
        MPI_Bcast(&nquery, 1, MPI_INT, query_rank, MPI_COMM_WORLD);

        if (myid != query_rank)
            queries.resize(nquery);

        if (nquery > 0)
        {
            MPI_Bcast(reinterpret_cast<void *>(queries.data()),
                      nquery * static_cast<int>(sizeof(DebugItem)),
                      MPI_BYTE,
                      query_rank,
                      MPI_COMM_WORLD);
        }

        if (myid == query_rank)
        {
            std::cout << "\n[DebugDumpPointPartners] field=" << label
                      << " loc=" << loc_name(desc.location)
                      << " src(rank=" << query_rank
                      << ", blk=" << blk
                      << ", idx=(" << i << "," << j << "," << k << "))"
                      << " ngh=" << ngh
                      << " include_topo=" << include_topo
                      << " include_halo=" << include_halo
                      << "\n";
            std::cout.flush();
        }

        for (int r = 0; r < nrank; ++r)
        {
            MPI_Barrier(MPI_COMM_WORLD);

            if (myid != r)
                continue;

            for (const auto &q : queries)
            {
                if (q.rank != myid)
                    continue;

                std::ostringstream oss;
                oss << std::scientific << std::setprecision(16);

                if (q.blk < 0 || q.blk >= fld_->num_blocks())
                {
                    oss << "[rank " << myid << "] "
                        << "[" << rel_name(q.relation) << "/" << patch_name(q.patch_tag) << "] "
                        << "field=" << label
                        << " target blk invalid: " << q.blk
                        << " <- src(rank=" << q.src_rank
                        << ", blk=" << q.src_blk
                        << ", idx=(" << q.src_i << "," << q.src_j << "," << q.src_k << "))\n";
                    std::cout << oss.str();
                    continue;
                }

                auto &F = fld_->field(fid, q.blk);
                if (!F.is_allocated())
                {
                    oss << "[rank " << myid << "] "
                        << "[" << rel_name(q.relation) << "/" << patch_name(q.patch_tag) << "] "
                        << "field=" << label
                        << " blk=" << q.blk
                        << " not allocated"
                        << " <- src(rank=" << q.src_rank
                        << ", blk=" << q.src_blk
                        << ", idx=(" << q.src_i << "," << q.src_j << "," << q.src_k << "))\n";
                    std::cout << oss.str();
                    continue;
                }

                Int3 lo = F.get_lo();
                Int3 hi = F.get_hi();

                if (q.i < lo.i || q.i >= hi.i ||
                    q.j < lo.j || q.j >= hi.j ||
                    q.k < lo.k || q.k >= hi.k)
                {
                    oss << "[rank " << myid << "] "
                        << "[" << rel_name(q.relation) << "/" << patch_name(q.patch_tag) << "] "
                        << "field=" << label
                        << " blk=" << q.blk
                        << " idx=(" << q.i << "," << q.j << "," << q.k << ") "
                        << "out_of_range, valid=[("
                        << lo.i << "," << lo.j << "," << lo.k << "),("
                        << hi.i << "," << hi.j << "," << hi.k << "))"
                        << " <- src(rank=" << q.src_rank
                        << ", blk=" << q.src_blk
                        << ", idx=(" << q.src_i << "," << q.src_j << "," << q.src_k << "))\n";
                    std::cout << oss.str();
                    continue;
                }

                double x, y, z;
                get_xyz(fid, q.blk, q.i, q.j, q.k, x, y, z);

                oss << "[rank " << myid << "] "
                    << "[" << rel_name(q.relation) << "/" << patch_name(q.patch_tag) << "] "
                    << "field=" << label
                    << " loc=" << loc_name(desc.location)
                    << " blk=" << q.blk
                    << " idx=(" << q.i << "," << q.j << "," << q.k << ")"
                    << " xyz=(" << x << "," << y << "," << z << ")"
                    << " <- src(rank=" << q.src_rank
                    << ", blk=" << q.src_blk
                    << ", idx=(" << q.src_i << "," << q.src_j << "," << q.src_k << "))"
                    << " : ";

                for (int m = 0; m < desc.ncomp; ++m)
                {
                    oss << F(q.i, q.j, q.k, m);
                    if (m + 1 < desc.ncomp)
                        oss << ", ";
                }
                oss << "\n";

                std::cout << oss.str();
            }

            std::cout.flush();
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }
}

void LunarSolver::Debug_TestJOperator_Manufactured(int test_id)
{
    constexpr int SAFE_PAD = 4;

    int mpi_rank = 0;
    int mpi_size = 1;
#ifdef MPI_VERSION
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#endif

    std::ostringstream fname;
    fname << "jtest_rank" << std::setw(4) << std::setfill('0') << mpi_rank << ".log";
    std::ofstream fout(fname.str(), std::ios::out);
    std::ostream &out = fout.is_open() ? fout : std::cout;

    out << "\n====================================================================================================\n";
    out << " Debug_TestJOperator_Manufactured START\n";
    out << " rank = " << mpi_rank << " / " << mpi_size << "\n";
    out << " SAFE_PAD = " << SAFE_PAD << "\n";
    out << " This test compares RAW_MODE(beta=1,no-alpha) and PRODUCTION_MODE(Calc_J_Edge).\n";
    out << "====================================================================================================\n";

    auto dot3 = [](const std::array<double, 3> &a,
                   const std::array<double, 3> &b) -> double
    {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    };

    auto add3 = [](const std::array<double, 3> &a,
                   const std::array<double, 3> &b) -> std::array<double, 3>
    {
        return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
    };

    auto sub3 = [](const std::array<double, 3> &a,
                   const std::array<double, 3> &b) -> std::array<double, 3>
    {
        return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
    };

    auto scale3 = [](double s,
                     const std::array<double, 3> &a) -> std::array<double, 3>
    {
        return {s * a[0], s * a[1], s * a[2]};
    };

    auto norm3 = [&](const std::array<double, 3> &a) -> double
    {
        return std::sqrt(dot3(a, a));
    };

    auto avg4 = [&](const std::array<double, 3> &a,
                    const std::array<double, 3> &b,
                    const std::array<double, 3> &c,
                    const std::array<double, 3> &d) -> std::array<double, 3>
    {
        return scale3(0.25, add3(add3(a, b), add3(c, d)));
    };

    auto avg8 = [&](const std::array<double, 3> &a,
                    const std::array<double, 3> &b,
                    const std::array<double, 3> &c,
                    const std::array<double, 3> &d,
                    const std::array<double, 3> &e,
                    const std::array<double, 3> &f,
                    const std::array<double, 3> &g,
                    const std::array<double, 3> &h) -> std::array<double, 3>
    {
        return scale3(0.125, add3(add3(add3(a, b), add3(c, d)),
                                  add3(add3(e, f), add3(g, h))));
    };

    auto analytic_b = [](int tid, double x, double y, double z) -> std::array<double, 3>
    {
        switch (tid)
        {
        case 0:
            return {1.0, 0.0, 0.0};
        case 1:
            return {0.0, 1.0, 0.0};
        case 2:
            return {0.0, 0.0, 1.0};

        case 3:
            return {0.0, -0.5 * z, 0.5 * y}; // curl = (1,0,0)
        case 4:
            return {0.5 * z, 0.0, -0.5 * x}; // curl = (0,1,0)
        case 5:
            return {-0.5 * y, 0.5 * x, 0.0}; // curl = (0,0,1)

        default:
            return {0.0, 0.0, 0.0};
        }
    };

    auto analytic_J = [](int tid, double, double, double) -> std::array<double, 3>
    {
        switch (tid)
        {
        case 0:
        case 1:
        case 2:
            return {0.0, 0.0, 0.0};

        case 3:
            return {1.0, 0.0, 0.0};
        case 4:
            return {0.0, 1.0, 0.0};
        case 5:
            return {0.0, 0.0, 1.0};

        default:
            return {0.0, 0.0, 0.0};
        }
    };

    auto test_name = [](int tid) -> std::string
    {
        switch (tid)
        {
        case 0:
            return "constant_bx";
        case 1:
            return "constant_by";
        case 2:
            return "constant_bz";
        case 3:
            return "linear_curl_Jx";
        case 4:
            return "linear_curl_Jy";
        case 5:
            return "linear_curl_Jz";
        default:
            return "unknown";
        }
    };

    // =============================================================================================
    // Geometry helpers
    //
    // 这里假设 grd_->grids(ib).x/y/z 是 node 坐标。
    // 如果你的节点坐标接口不是这个，只改 node_xyz 这个 lambda。
    // =============================================================================================
    auto node_xyz = [&](int ib, int i, int j, int k) -> std::array<double, 3>
    {
        auto &g = grd_->grids(ib);
        return {g.x(i, j, k), g.y(i, j, k), g.z(i, j, k)};
    };

    auto face_xi_center = [&](int ib, int i, int j, int k) -> std::array<double, 3>
    {
        return avg4(node_xyz(ib, i, j, k),
                    node_xyz(ib, i, j + 1, k),
                    node_xyz(ib, i, j, k + 1),
                    node_xyz(ib, i, j + 1, k + 1));
    };

    auto face_eta_center = [&](int ib, int i, int j, int k) -> std::array<double, 3>
    {
        return avg4(node_xyz(ib, i, j, k),
                    node_xyz(ib, i + 1, j, k),
                    node_xyz(ib, i, j, k + 1),
                    node_xyz(ib, i + 1, j, k + 1));
    };

    auto face_zeta_center = [&](int ib, int i, int j, int k) -> std::array<double, 3>
    {
        return avg4(node_xyz(ib, i, j, k),
                    node_xyz(ib, i + 1, j, k),
                    node_xyz(ib, i, j + 1, k),
                    node_xyz(ib, i + 1, j + 1, k));
    };

    auto cell_center = [&](int ib, int i, int j, int k) -> std::array<double, 3>
    {
        return avg8(node_xyz(ib, i, j, k),
                    node_xyz(ib, i + 1, j, k),
                    node_xyz(ib, i, j + 1, k),
                    node_xyz(ib, i + 1, j + 1, k),
                    node_xyz(ib, i, j, k + 1),
                    node_xyz(ib, i + 1, j, k + 1),
                    node_xyz(ib, i, j + 1, k + 1),
                    node_xyz(ib, i + 1, j + 1, k + 1));
    };

    auto edge_xi_mid_dr = [&](int ib, int i, int j, int k,
                              std::array<double, 3> &xm,
                              std::array<double, 3> &dr)
    {
        auto x0 = node_xyz(ib, i, j, k);
        auto x1 = node_xyz(ib, i + 1, j, k);
        dr = sub3(x1, x0);
        xm = add3(x0, scale3(0.5, dr));
    };

    auto edge_eta_mid_dr = [&](int ib, int i, int j, int k,
                               std::array<double, 3> &xm,
                               std::array<double, 3> &dr)
    {
        auto x0 = node_xyz(ib, i, j, k);
        auto x1 = node_xyz(ib, i, j + 1, k);
        dr = sub3(x1, x0);
        xm = add3(x0, scale3(0.5, dr));
    };

    auto edge_zeta_mid_dr = [&](int ib, int i, int j, int k,
                                std::array<double, 3> &xm,
                                std::array<double, 3> &dr)
    {
        auto x0 = node_xyz(ib, i, j, k);
        auto x1 = node_xyz(ib, i, j, k + 1);
        dr = sub3(x1, x0);
        xm = add3(x0, scale3(0.5, dr));
    };

    auto is_safe_interior = [&](const Int3 &lo, const Int3 &hi,
                                int i, int j, int k) -> bool
    {
        return (i >= lo.i + SAFE_PAD && i < hi.i - SAFE_PAD &&
                j >= lo.j + SAFE_PAD && j < hi.j - SAFE_PAD &&
                k >= lo.k + SAFE_PAD && k < hi.k - SAFE_PAD);
    };

    struct MaxRec
    {
        double val = -1.0;
        int ib = -1;
        int i = 0, j = 0, k = 0;
        std::array<double, 3> xyz{0.0, 0.0, 0.0};
        std::array<double, 3> num{0.0, 0.0, 0.0};
        std::array<double, 3> exact{0.0, 0.0, 0.0};

        void update(double e, int iblk, int ii, int jj, int kk,
                    const std::array<double, 3> &x,
                    const std::array<double, 3> &n,
                    const std::array<double, 3> &q)
        {
            if (e > val)
            {
                val = e;
                ib = iblk;
                i = ii;
                j = jj;
                k = kk;
                xyz = x;
                num = n;
                exact = q;
            }
        }
    };

    auto print_maxrec = [&](const std::string &label, const MaxRec &r)
    {
        out << std::scientific << std::setprecision(8);
        out << "  " << label << "\n";
        out << "    max_err = " << r.val
            << "  block=" << r.ib
            << "  idx=(" << r.i << "," << r.j << "," << r.k << ")"
            << "  xyz=(" << r.xyz[0] << "," << r.xyz[1] << "," << r.xyz[2] << ")\n";
        out << "    num   = (" << r.num[0] << "," << r.num[1] << "," << r.num[2] << ")\n";
        out << "    exact = (" << r.exact[0] << "," << r.exact[1] << "," << r.exact[2] << ")\n";
    };

    // =============================================================================================
    // Fill analytic Bface: B_face = b · A_face
    // =============================================================================================
    auto fill_analytic_Bface = [&](int tid)
    {
        for (int ib = 0; ib < fld_->num_blocks(); ++ib)
        {
            auto &Bxi = fld_->field(fid_.fid_B.xi, ib);
            auto &Bet = fld_->field(fid_.fid_B.eta, ib);
            auto &Bze = fld_->field(fid_.fid_B.zeta, ib);

            auto &Axi = fld_->field(fid_.fid_metric.xi, ib);
            auto &Aet = fld_->field(fid_.fid_metric.eta, ib);
            auto &Aze = fld_->field(fid_.fid_metric.zeta, ib);

            if (!Bxi.is_allocated() || !Bet.is_allocated() || !Bze.is_allocated())
                continue;
            if (!Axi.is_allocated() || !Aet.is_allocated() || !Aze.is_allocated())
                continue;

            {
                Int3 lo = Bxi.inner_lo();
                Int3 hi = Bxi.inner_hi();
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                        {
                            auto xc = face_xi_center(ib, i, j, k);
                            auto b = analytic_b(tid, xc[0], xc[1], xc[2]);
                            std::array<double, 3> A{Axi(i, j, k, 0), Axi(i, j, k, 1), Axi(i, j, k, 2)};
                            Bxi(i, j, k, 0) = dot3(b, A);
                        }
            }

            {
                Int3 lo = Bet.inner_lo();
                Int3 hi = Bet.inner_hi();
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                        {
                            auto xc = face_eta_center(ib, i, j, k);
                            auto b = analytic_b(tid, xc[0], xc[1], xc[2]);
                            std::array<double, 3> A{Aet(i, j, k, 0), Aet(i, j, k, 1), Aet(i, j, k, 2)};
                            Bet(i, j, k, 0) = dot3(b, A);
                        }
            }

            {
                Int3 lo = Bze.inner_lo();
                Int3 hi = Bze.inner_hi();
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                        {
                            auto xc = face_zeta_center(ib, i, j, k);
                            auto b = analytic_b(tid, xc[0], xc[1], xc[2]);
                            std::array<double, 3> A{Aze(i, j, k, 0), Aze(i, j, k, 1), Aze(i, j, k, 2)};
                            Bze(i, j, k, 0) = dot3(b, A);
                        }
            }
        }
    };

    auto clear_J_fields = [&]()
    {
        for (int ib = 0; ib < fld_->num_blocks(); ++ib)
        {
            auto &Jc = fld_->field(fid_.fid_Jcell, ib);
            if (Jc.is_allocated())
            {
                Int3 lo = Jc.inner_lo();
                Int3 hi = Jc.inner_hi();
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                            for (int m = 0; m < 3; ++m)
                                Jc(i, j, k, m) = 0.0;
            }

            auto &Jxi = fld_->field(fid_.fid_J.xi, ib);
            auto &Jet = fld_->field(fid_.fid_J.eta, ib);
            auto &Jze = fld_->field(fid_.fid_J.zeta, ib);

            if (Jxi.is_allocated())
            {
                Int3 lo = Jxi.inner_lo();
                Int3 hi = Jxi.inner_hi();
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                            Jxi(i, j, k, 0) = 0.0;
            }

            if (Jet.is_allocated())
            {
                Int3 lo = Jet.inner_lo();
                Int3 hi = Jet.inner_hi();
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                            Jet(i, j, k, 0) = 0.0;
            }

            if (Jze.is_allocated())
            {
                Int3 lo = Jze.inner_lo();
                Int3 hi = Jze.inner_hi();
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                            Jze(i, j, k, 0) = 0.0;
            }
        }
    };

    // =============================================================================================
    // RAW_MODE: beta=1, no alpha, no boundary sync.
    // This directly implements the same topological pattern as CurlAdjFaceToEdge.
    // =============================================================================================
    auto compute_Jedge_raw_beta1_noalpha = [&]()
    {
        clear_J_fields();

        for (int ib = 0; ib < fld_->num_blocks(); ++ib)
        {
            auto &Bxi = fld_->field(fid_.fid_B.xi, ib);
            auto &Bet = fld_->field(fid_.fid_B.eta, ib);
            auto &Bze = fld_->field(fid_.fid_B.zeta, ib);

            auto &Jxi = fld_->field(fid_.fid_J.xi, ib);
            auto &Jet = fld_->field(fid_.fid_J.eta, ib);
            auto &Jze = fld_->field(fid_.fid_J.zeta, ib);

            if (!Bxi.is_allocated() || !Bet.is_allocated() || !Bze.is_allocated())
                continue;
            if (!Jxi.is_allocated() || !Jet.is_allocated() || !Jze.is_allocated())
                continue;

            // Edge xi:
            // raw_xi = +Beta(i,j,k) - Beta(i,j,k-1) + Bzeta(i,j-1,k) - Bzeta(i,j,k)
            {
                Int3 lo = Jxi.inner_lo();
                Int3 hi = Jxi.inner_hi();

                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                        {
                            Jxi(i, j, k, 0) =
                                (Bet(i, j, k, 0) - Bet(i, j, k - 1, 0)) + (Bze(i, j - 1, k, 0) - Bze(i, j, k, 0));
                        }
            }

            // Edge eta:
            // raw_eta = +Bxi(i,j,k-1) - Bxi(i,j,k) + Bzeta(i,j,k) - Bzeta(i-1,j,k)
            {
                Int3 lo = Jet.inner_lo();
                Int3 hi = Jet.inner_hi();

                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                        {
                            Jet(i, j, k, 0) =
                                (Bxi(i, j, k - 1, 0) - Bxi(i, j, k, 0)) + (Bze(i, j, k, 0) - Bze(i - 1, j, k, 0));
                        }
            }

            // Edge zeta:
            // raw_zeta = +Bxi(i,j-1,k) - Bxi(i,j,k) + Beta(i,j,k) - Beta(i-1,j,k)
            {
                Int3 lo = Jze.inner_lo();
                Int3 hi = Jze.inner_hi();

                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                        {
                            Jze(i, j, k, 0) =
                                (Bxi(i, j - 1, k, 0) - Bxi(i, j, k, 0)) + (Bet(i, j, k, 0) - Bet(i - 1, j, k, 0));
                        }
            }
        }
    };

    // =============================================================================================
    // PRODUCTION_MODE: call your real Calc_J_Edge().
    // Note: if Calc_J_Edge() internally calls lunar_bound_.Sync("Jedge"), that is included here.
    // =============================================================================================
    auto compute_Jedge_production = [&]()
    {
        clear_J_fields();
        // Production singular-edge current is assembled from the physical
        // cell reconstruction of the manufactured face fluxes.
        calc_Bcell();
        Calc_J_Edge();
    };

    auto compute_Jcell_from_current_Jedge = [&]()
    {
        auto &Jc0 = fld_->field(fid_.fid_Jcell, 0);
        (void)Jc0;

        for (int ib = 0; ib < fld_->num_blocks(); ++ib)
        {
            auto &Jc = fld_->field(fid_.fid_Jcell, ib);
            if (!Jc.is_allocated())
                continue;

            Int3 lo = Jc.inner_lo();
            Int3 hi = Jc.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        for (int m = 0; m < 3; ++m)
                            Jc(i, j, k, m) = 0.0;
        }

        calc_Jcell();
    };

    auto check_Jedge_error = [&](int tid, const std::string &mode_label)
    {
        MaxRec all_xi, all_eta, all_zeta;
        MaxRec safe_xi, safe_eta, safe_zeta;

        double rms_all = 0.0, cnt_all = 0.0;
        double rms_safe = 0.0, cnt_safe = 0.0;

        for (int ib = 0; ib < fld_->num_blocks(); ++ib)
        {
            auto &Jxi = fld_->field(fid_.fid_J.xi, ib);
            auto &Jet = fld_->field(fid_.fid_J.eta, ib);
            auto &Jze = fld_->field(fid_.fid_J.zeta, ib);

            if (!Jxi.is_allocated() || !Jet.is_allocated() || !Jze.is_allocated())
                continue;

            {
                Int3 lo = Jxi.inner_lo();
                Int3 hi = Jxi.inner_hi();
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                        {
                            std::array<double, 3> xm, dr;
                            edge_xi_mid_dr(ib, i, j, k, xm, dr);
                            auto Je = analytic_J(tid, xm[0], xm[1], xm[2]);

                            const double exact = dot3(Je, dr);
                            const double num = Jxi(i, j, k, 0);
                            const double err = std::abs(num - exact);

                            std::array<double, 3> nv{num, 0.0, 0.0};
                            std::array<double, 3> ev{exact, 0.0, 0.0};

                            all_xi.update(err, ib, i, j, k, xm, nv, ev);
                            rms_all += err * err;
                            cnt_all += 1.0;

                            if (is_safe_interior(lo, hi, i, j, k))
                            {
                                safe_xi.update(err, ib, i, j, k, xm, nv, ev);
                                rms_safe += err * err;
                                cnt_safe += 1.0;
                            }
                        }
            }

            {
                Int3 lo = Jet.inner_lo();
                Int3 hi = Jet.inner_hi();
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                        {
                            std::array<double, 3> xm, dr;
                            edge_eta_mid_dr(ib, i, j, k, xm, dr);
                            auto Je = analytic_J(tid, xm[0], xm[1], xm[2]);

                            const double exact = dot3(Je, dr);
                            const double num = Jet(i, j, k, 0);
                            const double err = std::abs(num - exact);

                            std::array<double, 3> nv{num, 0.0, 0.0};
                            std::array<double, 3> ev{exact, 0.0, 0.0};

                            all_eta.update(err, ib, i, j, k, xm, nv, ev);
                            rms_all += err * err;
                            cnt_all += 1.0;

                            if (is_safe_interior(lo, hi, i, j, k))
                            {
                                safe_eta.update(err, ib, i, j, k, xm, nv, ev);
                                rms_safe += err * err;
                                cnt_safe += 1.0;
                            }
                        }
            }

            {
                Int3 lo = Jze.inner_lo();
                Int3 hi = Jze.inner_hi();
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                        {
                            std::array<double, 3> xm, dr;
                            edge_zeta_mid_dr(ib, i, j, k, xm, dr);
                            auto Je = analytic_J(tid, xm[0], xm[1], xm[2]);

                            const double exact = dot3(Je, dr);
                            const double num = Jze(i, j, k, 0);
                            const double err = std::abs(num - exact);

                            std::array<double, 3> nv{num, 0.0, 0.0};
                            std::array<double, 3> ev{exact, 0.0, 0.0};

                            all_zeta.update(err, ib, i, j, k, xm, nv, ev);
                            rms_all += err * err;
                            cnt_all += 1.0;

                            if (is_safe_interior(lo, hi, i, j, k))
                            {
                                safe_zeta.update(err, ib, i, j, k, xm, nv, ev);
                                rms_safe += err * err;
                                cnt_safe += 1.0;
                            }
                        }
            }
        }

        out << "\n[JTEST][" << test_name(tid) << "][" << mode_label << "] Jedge 1-form error\n";
        out << std::scientific << std::setprecision(8);
        out << "  RMS_all_inner     = " << std::sqrt(rms_all / std::max(1.0, cnt_all)) << "\n";
        out << "  RMS_safe_interior = " << std::sqrt(rms_safe / std::max(1.0, cnt_safe)) << "\n";

        print_maxrec("ALL  Jedge_xi   max_abs", all_xi);
        print_maxrec("ALL  Jedge_eta  max_abs", all_eta);
        print_maxrec("ALL  Jedge_zeta max_abs", all_zeta);

        print_maxrec("SAFE Jedge_xi   max_abs", safe_xi);
        print_maxrec("SAFE Jedge_eta  max_abs", safe_eta);
        print_maxrec("SAFE Jedge_zeta max_abs", safe_zeta);
    };

    auto check_Jcell_error = [&](int tid, const std::string &mode_label)
    {
        MaxRec all_abs, safe_abs;

        double rms_all = 0.0, cnt_all = 0.0;
        double rms_safe = 0.0, cnt_safe = 0.0;

        for (int ib = 0; ib < fld_->num_blocks(); ++ib)
        {
            auto &Jc = fld_->field(fid_.fid_Jcell, ib);
            if (!Jc.is_allocated())
                continue;

            Int3 lo = Jc.inner_lo();
            Int3 hi = Jc.inner_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        auto xc = cell_center(ib, i, j, k);
                        auto Je = analytic_J(tid, xc[0], xc[1], xc[2]);

                        std::array<double, 3> Jn{
                            Jc(i, j, k, 0),
                            Jc(i, j, k, 1),
                            Jc(i, j, k, 2)};

                        const double err = norm3(sub3(Jn, Je));

                        all_abs.update(err, ib, i, j, k, xc, Jn, Je);
                        rms_all += err * err;
                        cnt_all += 1.0;

                        if (is_safe_interior(lo, hi, i, j, k))
                        {
                            safe_abs.update(err, ib, i, j, k, xc, Jn, Je);
                            rms_safe += err * err;
                            cnt_safe += 1.0;
                        }
                    }
        }

        out << "\n[JTEST][" << test_name(tid) << "][" << mode_label << "] Jcell Cartesian error\n";
        out << std::scientific << std::setprecision(8);
        out << "  RMS_all_inner     = " << std::sqrt(rms_all / std::max(1.0, cnt_all)) << "\n";
        out << "  RMS_safe_interior = " << std::sqrt(rms_safe / std::max(1.0, cnt_safe)) << "\n";

        print_maxrec("ALL  Jcell max_abs", all_abs);
        print_maxrec("SAFE Jcell max_abs", safe_abs);
    };

    auto diagnose_Jcell_weights = [&](const std::string &mode_label)
    {
        out << "\\n[JTEST][" << mode_label
            << "] implicit-Hall weight diagnostics removed (explicit Hall only).\\n";
    };
    auto run_one_test = [&](int tid)
    {
        out << "\n----------------------------------------------------------------------------------------------------\n";
        out << "[JTEST] Running test " << tid << " : " << test_name(tid) << "\n";
        out << "----------------------------------------------------------------------------------------------------\n";

        fill_analytic_Bface(tid);

        // RAW_MODE
        {
            const std::string mode = "RAW_beta1_noalpha";

            compute_Jedge_raw_beta1_noalpha();
            check_Jedge_error(tid, mode);

            compute_Jcell_from_current_Jedge();
            check_Jcell_error(tid, mode);
            diagnose_Jcell_weights(mode);
        }

        // PRODUCTION_MODE
        {
            const std::string mode = "PRODUCTION_Calc_J_Edge";

            compute_Jedge_production();
            check_Jedge_error(tid, mode);

            compute_Jcell_from_current_Jedge();
            check_Jcell_error(tid, mode);
            diagnose_Jcell_weights(mode);
        }
    };

    if (test_id >= 0 && test_id <= 5)
    {
        run_one_test(test_id);
    }
    else
    {
        for (int tid = 0; tid <= 5; ++tid)
            run_one_test(tid);
    }

    out << "\n====================================================================================================\n";
    out << " Debug_TestJOperator_Manufactured END\n";
    out << " log file = " << fname.str() << "\n";
    out << "====================================================================================================\n";
}
