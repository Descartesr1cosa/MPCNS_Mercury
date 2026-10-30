#include "9_runtime/PartitionProfiler.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sys/stat.h>

namespace PROFILE
{
bool g_enabled = false;
double g_blocked_seconds = 0.0;

void PhaseScope::start_()
{
    b0_ = g_blocked_seconds;
    t0_ = MPI_Wtime();
}

void PhaseScope::stop_()
{
    const double wall = MPI_Wtime() - t0_;
    p_.add_phase(id_, wall, g_blocked_seconds - b0_);
}

void PartitionProfiler::configure(bool enabled, int warmup_steps, int measured_steps)
{
    g_enabled = enabled;
    warmup_ = std::max(0, warmup_steps);
    measured_ = std::max(1, measured_steps);
    step_ = 0;
    measuring_ = enabled && warmup_ == 0;
    for (int i = 0; i < kMaxPhases; ++i)
        wall_[i] = blocked_[i] = 0.0;
    g_blocked_seconds = 0.0;
}

bool PartitionProfiler::end_step()
{
    if (!g_enabled)
        return false;
    ++step_;
    // Phases of step `step_` are recorded iff the step is inside the window;
    // measuring_ is switched for the *next* step here.
    measuring_ = (step_ >= warmup_) && (step_ < warmup_ + measured_);
    return step_ >= warmup_ + measured_;
}

namespace
{
struct Fit
{
    double a_fluid = 0.0, a_solid = 0.0, c = 0.0;
    double r2 = 0.0, rms = 0.0, cond = 0.0;
    bool ok = false;
};

// Solve the 3x3 normal equations restricted to the `active` columns.
bool solve_subset(const std::vector<std::vector<double>> &X, const std::vector<double> &y,
                  const bool active[3], double coef[3])
{
    int idx[3], n = 0;
    for (int j = 0; j < 3; ++j)
        if (active[j])
            idx[n++] = j;
    for (int j = 0; j < 3; ++j)
        coef[j] = 0.0;
    if (n == 0)
        return true;
    double A[3][4] = {};
    for (size_t r = 0; r < X.size(); ++r)
        for (int p = 0; p < n; ++p)
        {
            for (int q = 0; q < n; ++q)
                A[p][q] += X[r][idx[p]] * X[r][idx[q]];
            A[p][3] += X[r][idx[p]] * y[r];
        }
    for (int col = 0; col < n; ++col)
    {
        int piv = col;
        for (int r = col + 1; r < n; ++r)
            if (std::fabs(A[r][col]) > std::fabs(A[piv][col]))
                piv = r;
        if (std::fabs(A[piv][col]) < 1e-300)
            return false;
        for (int k = 0; k < 4; ++k)
            std::swap(A[col][k], A[piv][k]);
        for (int r = 0; r < n; ++r)
        {
            if (r == col)
                continue;
            const double f = A[r][col] / A[col][col];
            for (int k = col; k < 4; ++k)
                A[r][k] -= f * A[col][k];
        }
    }
    for (int p = 0; p < n; ++p)
        coef[idx[p]] = A[p][3] / A[p][p];
    return true;
}

// Non-negative least squares over the 3 columns by enumerating active sets.
Fit fit_nonneg(const std::vector<std::vector<double>> &X, const std::vector<double> &y)
{
    Fit best;
    double best_rss = 1e300;
    for (int mask = 1; mask < 8; ++mask)
    {
        const bool active[3] = {(mask & 1) != 0, (mask & 2) != 0, (mask & 4) != 0};
        double coef[3];
        if (!solve_subset(X, y, active, coef))
            continue;
        if (coef[0] < 0 || coef[1] < 0 || coef[2] < 0)
            continue;
        double rss = 0.0;
        for (size_t r = 0; r < X.size(); ++r)
        {
            const double e = y[r] - (coef[0] * X[r][0] + coef[1] * X[r][1] + coef[2] * X[r][2]);
            rss += e * e;
        }
        if (rss < best_rss)
        {
            best_rss = rss;
            best.a_fluid = coef[0];
            best.a_solid = coef[1];
            best.c = coef[2];
            best.ok = true;
        }
    }
    if (!best.ok)
        return best;
    double mean = 0.0;
    for (double v : y)
        mean += v;
    mean /= std::max<size_t>(1, y.size());
    double tss = 0.0;
    for (double v : y)
        tss += (v - mean) * (v - mean);
    best.rms = std::sqrt(best_rss / std::max<size_t>(1, y.size()));
    best.r2 = tss > 0.0 ? 1.0 - best_rss / tss : 0.0;

    // Condition number of the column-normalised design matrix (3x3 Gram).
    double colnorm[3] = {};
    for (const auto &row : X)
        for (int j = 0; j < 3; ++j)
            colnorm[j] += row[j] * row[j];
    double G[3][3] = {};
    for (const auto &row : X)
        for (int p = 0; p < 3; ++p)
            for (int q = 0; q < 3; ++q)
            {
                const double sp = colnorm[p] > 0 ? 1.0 / std::sqrt(colnorm[p]) : 0.0;
                const double sq = colnorm[q] > 0 ? 1.0 / std::sqrt(colnorm[q]) : 0.0;
                G[p][q] += row[p] * sp * row[q] * sq;
            }
    // power iteration for the largest / smallest eigenvalue of a 3x3 SPD matrix
    auto matvec = [&](const double *v, double *o)
    { for (int p = 0; p < 3; ++p) { o[p] = 0; for (int q = 0; q < 3; ++q) o[p] += G[p][q] * v[q]; } };
    double v[3] = {1, 1, 1}, w[3], lmax = 0.0;
    for (int it = 0; it < 200; ++it)
    {
        matvec(v, w);
        double nrm = std::sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
        if (nrm == 0.0)
            break;
        for (int p = 0; p < 3; ++p)
            v[p] = w[p] / nrm;
        lmax = nrm;
    }
    // smallest eigenvalue via shifted power iteration on (lmax*I - G)
    double u[3] = {1, 0.5, 0.25}, lmin = lmax;
    for (int it = 0; it < 500; ++it)
    {
        matvec(u, w);
        for (int p = 0; p < 3; ++p)
            w[p] = lmax * u[p] - w[p];
        double nrm = std::sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
        if (nrm == 0.0)
            break;
        for (int p = 0; p < 3; ++p)
            u[p] = w[p] / nrm;
        double gu[3];
        matvec(u, gu);
        lmin = u[0] * gu[0] + u[1] * gu[1] + u[2] * gu[2];
    }
    best.cond = lmin > 1e-300 ? lmax / lmin : 1e300;
    return best;
}
} // namespace

void PartitionProfiler::write_json(const std::string &path, int nranks, int rank) const
{
    const int nphase = [&]
    {
        int n = 0;
        for (int i = 0; i < kMaxPhases; ++i)
            if (!names_[i].empty())
                n = i + 1;
        return n;
    }();

    // Per-rank record: [fluid_cells, solid_cells, alloc_fluid, alloc_solid, blocks,
    //                   measured_steps, then (wall, blocked) per phase]
    const int stride = 6 + 2 * nphase;
    std::vector<double> mine(stride, 0.0);
    for (const auto &b : blocks_)
    {
        (b.fluid ? mine[0] : mine[1]) += static_cast<double>(b.cells);
        (b.fluid ? mine[2] : mine[3]) += static_cast<double>(b.alloc_cells);
    }
    mine[4] = static_cast<double>(blocks_.size());
    mine[5] = static_cast<double>(measured_);
    for (int i = 0; i < nphase; ++i)
    {
        mine[6 + 2 * i] = wall_[i];
        mine[7 + 2 * i] = blocked_[i];
    }
    std::vector<double> all(rank == 0 ? static_cast<size_t>(stride) * nranks : 0);
    PMPI_Gather(mine.data(), stride, MPI_DOUBLE, all.data(), stride, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Per-block names/cells gathered as text so the JSON lists every block.
    std::string blocks_txt;
    for (size_t i = 0; i < blocks_.size(); ++i)
        blocks_txt += blocks_[i].name + "," + std::to_string(blocks_[i].cells) + "," +
                      std::to_string(blocks_[i].alloc_cells) + ";";
    int len = static_cast<int>(blocks_txt.size());
    std::vector<int> lens(rank == 0 ? nranks : 0);
    PMPI_Gather(&len, 1, MPI_INT, lens.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    std::vector<int> displs(rank == 0 ? nranks : 0);
    int total = 0;
    if (rank == 0)
        for (int r = 0; r < nranks; ++r)
        {
            displs[r] = total;
            total += lens[r];
        }
    std::vector<char> text(rank == 0 ? total : 0);
    PMPI_Gatherv(blocks_txt.data(), len, MPI_CHAR, text.data(), lens.data(), displs.data(),
                 MPI_CHAR, 0, MPI_COMM_WORLD);
    if (rank != 0)
        return;

    auto steps_of = [&](int r) { return all[static_cast<size_t>(r) * stride + 5]; };
    // regression data: per rank, per-step compute (s/step)
    std::vector<std::vector<double>> X;
    std::vector<double> y_total;
    std::vector<std::vector<double>> y_phase(nphase);
    for (int r = 0; r < nranks; ++r)
    {
        const double *p = &all[static_cast<size_t>(r) * stride];
        X.push_back({p[2], p[3], p[4]}); // ghost-inclusive cells by type, blocks
        double total_compute = 0.0;
        for (int i = 0; i < nphase; ++i)
        {
            const double c = (p[6 + 2 * i] - p[7 + 2 * i]) / steps_of(r);
            y_phase[i].push_back(c);
            total_compute += c;
        }
        y_total.push_back(total_compute);
    }
    const Fit total_fit = fit_nonneg(X, y_total);

    mkdir("DATA", 0755);
    std::ofstream out(path);
    if (!out)
    {
        std::fprintf(stderr, "[PartitionProfiler] cannot write %s\n", path.c_str());
        return;
    }
    out.precision(10);
    const double w = total_fit.a_solid > 0.0 ? total_fit.a_fluid / total_fit.a_solid : 1.0;
    const double oh = total_fit.a_solid > 0.0 ? total_fit.c / total_fit.a_solid : 0.0;
    const bool well = total_fit.ok && total_fit.cond < 1.0e3 && total_fit.r2 > 0.8 &&
                      total_fit.a_solid > 0.0 && total_fit.a_fluid > 0.0;
    out << "{\n";
    out << "  \"schema\": \"mpcns_partition_profile_v1\",\n";
    out << "  \"ranks\": " << nranks << ",\n";
    out << "  \"measured_steps\": " << measured_ << ",\n";
    out << "  \"model\": \"compute_s_per_step = a_fluid*alloc_fluid_cells + a_solid*alloc_solid_cells + c_block*blocks (alloc cells include ghost layers)\",\n";
    out << "  \"recommended\": {\n";
    out << "    \"fluid_cost_weight\": " << w << ",\n";
    out << "    \"block_cost_overhead\": " << oh << ",\n";
    out << "    \"cost_ghost_layers\": " << ghost_layers_ << ",\n";
    out << "    \"well_conditioned\": " << (well ? "true" : "false") << "\n";
    out << "  },\n";
    out << "  \"fit\": {\"a_fluid\": " << total_fit.a_fluid << ", \"a_solid\": " << total_fit.a_solid
        << ", \"c_block\": " << total_fit.c << ", \"r2\": " << total_fit.r2 << ", \"rms\": "
        << total_fit.rms << ", \"condition\": " << total_fit.cond << "},\n";
    out << "  \"phases\": [\n";
    for (int i = 0; i < nphase; ++i)
    {
        const Fit f = fit_nonneg(X, y_phase[i]);
        double mean_c = 0.0, max_c = 0.0;
        for (double v : y_phase[i])
        {
            mean_c += v;
            max_c = std::max(max_c, v);
        }
        mean_c /= nranks;
        out << "    {\"name\": \"" << names_[i] << "\", \"mean_compute\": " << mean_c
            << ", \"max_compute\": " << max_c << ", \"a_fluid\": " << f.a_fluid << ", \"a_solid\": "
            << f.a_solid << ", \"c_block\": " << f.c << ", \"r2\": " << f.r2 << "}"
            << (i + 1 < nphase ? "," : "") << "\n";
    }
    out << "  ],\n";
    out << "  \"per_rank\": [\n";
    size_t off = 0;
    for (int r = 0; r < nranks; ++r)
    {
        const double *p = &all[static_cast<size_t>(r) * stride];
        double wall = 0.0, blocked = 0.0;
        for (int i = 0; i < nphase; ++i)
        {
            wall += p[6 + 2 * i];
            blocked += p[7 + 2 * i];
        }
        out << "    {\"rank\": " << r << ", \"fluid_cells\": " << p[0] << ", \"solid_cells\": " << p[1]
            << ", \"alloc_fluid_cells\": " << p[2] << ", \"alloc_solid_cells\": " << p[3]
            << ", \"blocks\": " << p[4] << ", \"step_wall\": " << wall / steps_of(r)
            << ", \"step_blocked\": " << blocked / steps_of(r) << ", \"step_compute\": "
            << (wall - blocked) / steps_of(r) << ", \"block_list\": \""
            << std::string(text.begin() + off, text.begin() + off + lens[r]) << "\"}"
            << (r + 1 < nranks ? "," : "") << "\n";
        off += lens[r];
    }
    out << "  ]\n}\n";
    std::printf("[PartitionProfiler] wrote %s: fluid_cost_weight=%.4g block_cost_overhead=%.4g "
                "(R2=%.3f, cond=%.3g, %s)\n",
                path.c_str(), w, oh, total_fit.r2, total_fit.cond,
                well ? "well conditioned" : "WEAK FIT - use a more varied partition");
}
} // namespace PROFILE

//==============================================================================
// PMPI wrappers: accumulate time spent blocked in synchronising MPI calls.
// They only measure when profiling is enabled; otherwise they forward directly.
//==============================================================================
#undef MPI_Waitall
#undef MPI_Allreduce
#undef MPI_Allgather
#undef MPI_Allgatherv
#undef MPI_Alltoall
#undef MPI_Alltoallv
#undef MPI_Barrier
#undef MPI_Bcast
#undef MPI_Gather
#undef MPI_Exscan
#undef MPI_Reduce

#define PROFILE_BLOCKED(call)                                         \
    if (!PROFILE::g_enabled)                                          \
        return call;                                                  \
    const double t0_ = PMPI_Wtime();                                  \
    const int err_ = call;                                            \
    PROFILE::g_blocked_seconds += PMPI_Wtime() - t0_;                 \
    return err_

extern "C"
{
int MPI_Waitall(int n, MPI_Request r[], MPI_Status s[])
{
    PROFILE_BLOCKED(PMPI_Waitall(n, r, s));
}
int MPI_Allreduce(const void *a, void *b, int c, MPI_Datatype d, MPI_Op o, MPI_Comm m)
{
    PROFILE_BLOCKED(PMPI_Allreduce(a, b, c, d, o, m));
}
int MPI_Reduce(const void *a, void *b, int c, MPI_Datatype d, MPI_Op o, int root, MPI_Comm m)
{
    PROFILE_BLOCKED(PMPI_Reduce(a, b, c, d, o, root, m));
}
int MPI_Allgather(const void *a, int ac, MPI_Datatype ad, void *b, int bc, MPI_Datatype bd, MPI_Comm m)
{
    PROFILE_BLOCKED(PMPI_Allgather(a, ac, ad, b, bc, bd, m));
}
int MPI_Allgatherv(const void *a, int ac, MPI_Datatype ad, void *b, const int bc[], const int disp[],
                   MPI_Datatype bd, MPI_Comm m)
{
    PROFILE_BLOCKED(PMPI_Allgatherv(a, ac, ad, b, bc, disp, bd, m));
}
int MPI_Alltoall(const void *a, int ac, MPI_Datatype ad, void *b, int bc, MPI_Datatype bd, MPI_Comm m)
{
    PROFILE_BLOCKED(PMPI_Alltoall(a, ac, ad, b, bc, bd, m));
}
int MPI_Alltoallv(const void *a, const int ac[], const int ad[], MPI_Datatype adt, void *b, const int bc[],
                  const int bd[], MPI_Datatype bdt, MPI_Comm m)
{
    PROFILE_BLOCKED(PMPI_Alltoallv(a, ac, ad, adt, b, bc, bd, bdt, m));
}
int MPI_Barrier(MPI_Comm m)
{
    PROFILE_BLOCKED(PMPI_Barrier(m));
}
int MPI_Bcast(void *b, int c, MPI_Datatype d, int root, MPI_Comm m)
{
    PROFILE_BLOCKED(PMPI_Bcast(b, c, d, root, m));
}
int MPI_Gather(const void *a, int ac, MPI_Datatype ad, void *b, int bc, MPI_Datatype bd, int root, MPI_Comm m)
{
    PROFILE_BLOCKED(PMPI_Gather(a, ac, ad, b, bc, bd, root, m));
}
int MPI_Exscan(const void *a, void *b, int c, MPI_Datatype d, MPI_Op o, MPI_Comm m)
{
    PROFILE_BLOCKED(PMPI_Exscan(a, b, c, d, o, m));
}
}
