#include "7_metric/Metric.h"

#include "3_field/Field.h"
#include "mpcns/compat/legacy_geometry_field_registration.hpp"
#include "mpcns/compat/legacy_lumped_metric_publisher.hpp"
#include "mpcns/compat/legacy_primal_metric_publisher.hpp"
#include "mpcns/compat/legacy_reconstruction_metric_publisher.hpp"

#include <limits>

namespace METRIC
{
void register_metric_fields(Field &fields, int geomtry_ghost_)
{
    mpcns::compat::register_legacy_geometry_fields(fields, geomtry_ghost_);
}

void compute_metric_fields(Field &fields, Grid &grid)
{
    // Dependency order is intentional; compatibility publishers retain legacy storage.
    // Do not reorder without an exact field-by-field regression gate.
    mpcns::compat::populate_legacy_primal_face_and_cell_metrics(fields, grid);

    mpcns::compat::populate_legacy_lumped_face_and_edge_metrics(fields, grid);

    mpcns::compat::populate_legacy_reconstruction_cache(fields, grid);

    mpcns::compat::populate_legacy_primal_edge_vectors(fields, grid);

    mpcns::compat::populate_legacy_diagonal_hodge_closure(fields, grid);

}

void build_field_geometry(Field &fields, Grid &grid, int geometry_ghost)
{
    register_metric_fields(fields, geometry_ghost);
    compute_metric_fields(fields, grid);
}

MetricDiagnostics diagnose_metric_fields(Field &fields)
{
    MetricDiagnostics diag;

    auto scan_positive = [](FieldBlock &fb, std::int64_t &counter)
    {
        Int3 lo = fb.get_lo();
        Int3 hi = fb.get_hi();
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                    if (!(fb(i, j, k, 0) > 0.0))
                        ++counter;
    };

    auto scan_finite = [](FieldBlock &fb, std::int64_t &counter)
    {
        Int3 lo = fb.get_lo();
        Int3 hi = fb.get_hi();
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                    if (!std::isfinite(fb(i, j, k, 0)))
                        ++counter;
    };

    auto scan_zero = [](FieldBlock &fb, std::int64_t &counter)
    {
        Int3 lo = fb.get_lo();
        Int3 hi = fb.get_hi();
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                    if (fb(i, j, k, 0) == 0.0)
                        ++counter;
    };

    for (int ib = 0; ib < fields.num_blocks(); ++ib)
    {
        scan_positive(fields.field("Jac", ib), diag.jac_nonpositive);

        scan_positive(fields.field("Area_xi", ib), diag.area_nonpositive);
        scan_positive(fields.field("Area_eta", ib), diag.area_nonpositive);
        scan_positive(fields.field("Area_zeta", ib), diag.area_nonpositive);

        scan_positive(fields.field("dl_xi", ib), diag.dl_nonpositive);
        scan_positive(fields.field("dl_eta", ib), diag.dl_nonpositive);
        scan_positive(fields.field("dl_zeta", ib), diag.dl_nonpositive);

        scan_finite(fields.field("Hodge_star_inverse_2form_to_1form_edge_xi_lumped", ib),
                    diag.Hodge_star_inverse_2form_to_1form_nonfinite);
        scan_finite(fields.field("Hodge_star_inverse_2form_to_1form_edge_eta_lumped", ib),
                    diag.Hodge_star_inverse_2form_to_1form_nonfinite);
        scan_finite(fields.field("Hodge_star_inverse_2form_to_1form_edge_zeta_lumped", ib),
                    diag.Hodge_star_inverse_2form_to_1form_nonfinite);

        scan_finite(fields.field("Hodge_star_2form_to_1form_face_xi_lumped", ib),
                    diag.Hodge_star_2form_to_1form_nonfinite);
        scan_finite(fields.field("Hodge_star_2form_to_1form_face_eta_lumped", ib),
                    diag.Hodge_star_2form_to_1form_nonfinite);
        scan_finite(fields.field("Hodge_star_2form_to_1form_face_zeta_lumped", ib),
                    diag.Hodge_star_2form_to_1form_nonfinite);

        scan_zero(fields.field("Astar_xi", ib), diag.near_axis_singular);
        scan_zero(fields.field("Astar_eta", ib), diag.near_axis_singular);
        scan_zero(fields.field("Astar_zeta", ib), diag.near_axis_singular);

        scan_zero(fields.field("Hodge_star_inverse_2form_to_1form_edge_xi_lumped", ib), diag.near_axis_capped);
        scan_zero(fields.field("Hodge_star_inverse_2form_to_1form_edge_eta_lumped", ib), diag.near_axis_capped);
        scan_zero(fields.field("Hodge_star_inverse_2form_to_1form_edge_zeta_lumped", ib), diag.near_axis_capped);
    }

    return diag;
}

namespace
{
bool metric_in_range(const FieldBlock &fb, int i, int j, int k)
{
    const Int3 lo = fb.get_lo();
    const Int3 hi = fb.get_hi();
    return i >= lo.i && i < hi.i &&
           j >= lo.j && j < hi.j &&
           k >= lo.k && k < hi.k;
}

void metric_zero(double out_xyz[3])
{
    out_xyz[0] = 0.0;
    out_xyz[1] = 0.0;
    out_xyz[2] = 0.0;
}

bool metric_solve3x3(double A[3][3], double b[3], double x[3])
{
    double M[3][4] = {
        {A[0][0], A[0][1], A[0][2], b[0]},
        {A[1][0], A[1][1], A[1][2], b[1]},
        {A[2][0], A[2][1], A[2][2], b[2]}};

    for (int c = 0; c < 3; ++c)
    {
        int piv = c;
        double amax = std::fabs(M[c][c]);
        for (int r = c + 1; r < 3; ++r)
        {
            const double av = std::fabs(M[r][c]);
            if (av > amax)
            {
                amax = av;
                piv = r;
            }
        }
        if (amax < 1.0e-300)
            return false;
        if (piv != c)
            for (int q = c; q < 4; ++q)
                std::swap(M[c][q], M[piv][q]);

        const double inv = 1.0 / M[c][c];
        for (int q = c; q < 4; ++q)
            M[c][q] *= inv;

        for (int r = 0; r < 3; ++r)
        {
            if (r == c)
                continue;
            const double f = M[r][c];
            for (int q = c; q < 4; ++q)
                M[r][q] -= f * M[c][q];
        }
    }

    x[0] = M[0][3];
    x[1] = M[1][3];
    x[2] = M[2][3];
    return true;
}
}

bool reconstruct_face_2form_to_cell(Field &fields,
                                    int fid_xi, int fid_eta, int fid_zeta,
                                    int iblock, int i, int j, int k,
                                    double out_xyz[3])
{
    metric_zero(out_xyz);

    FieldBlock &Fxi = fields.field(fid_xi, iblock);
    FieldBlock &Feta = fields.field(fid_eta, iblock);
    FieldBlock &Fzeta = fields.field(fid_zeta, iblock);

    if (!Fxi.is_allocated() || !Feta.is_allocated() || !Fzeta.is_allocated())
        return false;

    if (!metric_in_range(Fxi, i, j, k) ||
        !metric_in_range(Fxi, i + 1, j, k) ||
        !metric_in_range(Feta, i, j, k) ||
        !metric_in_range(Feta, i, j + 1, k) ||
        !metric_in_range(Fzeta, i, j, k) ||
        !metric_in_range(Fzeta, i, j, k + 1))
        return false;

    const double phi[6] = {
        -Fxi(i, j, k, 0),
        Fxi(i + 1, j, k, 0),
        -Feta(i, j, k, 0),
        Feta(i, j + 1, k, 0),
        -Fzeta(i, j, k, 0),
        Fzeta(i, j, k + 1, 0)};

    if (fields.has_field("Bcell_from_Bface_w"))
    {
        FieldBlock &W = fields.field("Bcell_from_Bface_w", iblock);
        if (W.is_allocated() && metric_in_range(W, i, j, k))
        {
            for (int n = 0; n < 6; ++n)
            {
                out_xyz[0] += W(i, j, k, n) * phi[n];
                out_xyz[1] += W(i, j, k, 6 + n) * phi[n];
                out_xyz[2] += W(i, j, k, 12 + n) * phi[n];
            }
            return true;
        }
    }

    if (!fields.has_field("JDxi") || !fields.has_field("JDet") || !fields.has_field("JDze"))
        return false;

    FieldBlock &JDxi = fields.field("JDxi", iblock);
    FieldBlock &JDet = fields.field("JDet", iblock);
    FieldBlock &JDze = fields.field("JDze", iblock);
    if (!JDxi.is_allocated() || !JDet.is_allocated() || !JDze.is_allocated())
        return false;
    if (!metric_in_range(JDxi, i, j, k) || !metric_in_range(JDxi, i + 1, j, k) ||
        !metric_in_range(JDet, i, j, k) || !metric_in_range(JDet, i, j + 1, k) ||
        !metric_in_range(JDze, i, j, k) || !metric_in_range(JDze, i, j, k + 1))
        return false;

    double A[3][3] = {
        {0.5 * (JDxi(i, j, k, 0) + JDxi(i + 1, j, k, 0)),
         0.5 * (JDxi(i, j, k, 1) + JDxi(i + 1, j, k, 1)),
         0.5 * (JDxi(i, j, k, 2) + JDxi(i + 1, j, k, 2))},
        {0.5 * (JDet(i, j, k, 0) + JDet(i, j + 1, k, 0)),
         0.5 * (JDet(i, j, k, 1) + JDet(i, j + 1, k, 1)),
         0.5 * (JDet(i, j, k, 2) + JDet(i, j + 1, k, 2))},
        {0.5 * (JDze(i, j, k, 0) + JDze(i, j, k + 1, 0)),
         0.5 * (JDze(i, j, k, 1) + JDze(i, j, k + 1, 1)),
         0.5 * (JDze(i, j, k, 2) + JDze(i, j, k + 1, 2))}};

    double rhs[3] = {
        0.5 * (Fxi(i, j, k, 0) + Fxi(i + 1, j, k, 0)),
        0.5 * (Feta(i, j, k, 0) + Feta(i, j + 1, k, 0)),
        0.5 * (Fzeta(i, j, k, 0) + Fzeta(i, j, k + 1, 0))};

    return metric_solve3x3(A, rhs, out_xyz);
}

bool reconstruct_edge_1form_to_cell(Field &fields,
                                    int fid_xi, int fid_eta, int fid_zeta,
                                    int iblock, int i, int j, int k,
                                    double out_xyz[3])
{
    metric_zero(out_xyz);

    FieldBlock &Exi = fields.field(fid_xi, iblock);
    FieldBlock &Eeta = fields.field(fid_eta, iblock);
    FieldBlock &Ezeta = fields.field(fid_zeta, iblock);

    if (!Exi.is_allocated() || !Eeta.is_allocated() || !Ezeta.is_allocated())
        return false;

    if (!metric_in_range(Exi, i, j, k) ||
        !metric_in_range(Exi, i, j + 1, k) ||
        !metric_in_range(Exi, i, j, k + 1) ||
        !metric_in_range(Exi, i, j + 1, k + 1) ||
        !metric_in_range(Eeta, i, j, k) ||
        !metric_in_range(Eeta, i + 1, j, k) ||
        !metric_in_range(Eeta, i, j, k + 1) ||
        !metric_in_range(Eeta, i + 1, j, k + 1) ||
        !metric_in_range(Ezeta, i, j, k) ||
        !metric_in_range(Ezeta, i + 1, j, k) ||
        !metric_in_range(Ezeta, i, j + 1, k) ||
        !metric_in_range(Ezeta, i + 1, j + 1, k))
        return false;

    const double s[12] = {
        Exi(i, j, k, 0),
        Exi(i, j + 1, k, 0),
        Exi(i, j, k + 1, 0),
        Exi(i, j + 1, k + 1, 0),
        Eeta(i, j, k, 0),
        Eeta(i + 1, j, k, 0),
        Eeta(i, j, k + 1, 0),
        Eeta(i + 1, j, k + 1, 0),
        Ezeta(i, j, k, 0),
        Ezeta(i + 1, j, k, 0),
        Ezeta(i, j + 1, k, 0),
        Ezeta(i + 1, j + 1, k, 0)};

    if (fields.has_field("Jcell_from_Jedge_w"))
    {
        FieldBlock &W = fields.field("Jcell_from_Jedge_w", iblock);
        if (W.is_allocated() && metric_in_range(W, i, j, k))
        {
            for (int n = 0; n < 12; ++n)
            {
                out_xyz[0] += W(i, j, k, n) * s[n];
                out_xyz[1] += W(i, j, k, 12 + n) * s[n];
                out_xyz[2] += W(i, j, k, 24 + n) * s[n];
            }
            return true;
        }
    }

    if (!fields.has_field("pinvGT_cell"))
        return false;
    FieldBlock &pinv = fields.field("pinvGT_cell", iblock);
    if (!pinv.is_allocated() || !metric_in_range(pinv, i, j, k))
        return false;

    const double w[3] = {
        0.25 * (s[0] + s[1] + s[2] + s[3]),
        0.25 * (s[4] + s[5] + s[6] + s[7]),
        0.25 * (s[8] + s[9] + s[10] + s[11])};

    for (int r = 0; r < 3; ++r)
        out_xyz[r] = pinv(i, j, k, 3 * r + 0) * w[0] +
                     pinv(i, j, k, 3 * r + 1) * w[1] +
                     pinv(i, j, k, 3 * r + 2) * w[2];

    return true;
}

bool reconstruct_face_2form_to_node(Field &fields,
                                    int fid_xi, int fid_eta, int fid_zeta,
                                    int iblock, int i, int j, int k,
                                    double out_xyz[3])
{
    metric_zero(out_xyz);

    const FieldBlock &ref = fields.field(fid_xi, iblock);
    const Block &blk = ref.get_block();
    const int dim = blk.dimension;

    double sum[3] = {0.0, 0.0, 0.0};
    int count = 0;
    for (int di = -1; di <= 0; ++di)
        for (int dj = -1; dj <= 0; ++dj)
            for (int dk = (dim == 3 ? -1 : 0); dk <= 0; ++dk)
            {
                if (dim == 2 && dk != 0)
                    continue;
                double cell_xyz[3];
                if (!reconstruct_face_2form_to_cell(fields, fid_xi, fid_eta, fid_zeta,
                                                     iblock, i + di, j + dj, k + dk, cell_xyz))
                    continue;
                sum[0] += cell_xyz[0];
                sum[1] += cell_xyz[1];
                sum[2] += cell_xyz[2];
                ++count;
            }

    if (count == 0)
        return false;

    out_xyz[0] = sum[0] / count;
    out_xyz[1] = sum[1] / count;
    out_xyz[2] = sum[2] / count;
    return true;
}

bool reconstruct_edge_1form_to_node(Field &fields,
                                    int fid_xi, int fid_eta, int fid_zeta,
                                    int iblock, int i, int j, int k,
                                    double out_xyz[3])
{
    metric_zero(out_xyz);

    const FieldBlock &ref = fields.field(fid_xi, iblock);
    const Block &blk = ref.get_block();
    const int dim = blk.dimension;

    double sum[3] = {0.0, 0.0, 0.0};
    int count = 0;
    for (int di = -1; di <= 0; ++di)
        for (int dj = -1; dj <= 0; ++dj)
            for (int dk = (dim == 3 ? -1 : 0); dk <= 0; ++dk)
            {
                if (dim == 2 && dk != 0)
                    continue;
                double cell_xyz[3];
                if (!reconstruct_edge_1form_to_cell(fields, fid_xi, fid_eta, fid_zeta,
                                                     iblock, i + di, j + dj, k + dk, cell_xyz))
                    continue;
                sum[0] += cell_xyz[0];
                sum[1] += cell_xyz[1];
                sum[2] += cell_xyz[2];
                ++count;
            }

    if (count == 0)
        return false;

    out_xyz[0] = sum[0] / count;
    out_xyz[1] = sum[1] / count;
    out_xyz[2] = sum[2] / count;
    return true;
}
}
