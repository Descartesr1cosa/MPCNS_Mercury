#include "mpcns/operators/pole_current_constraint.hpp"

#include "mpcns/constraints/least_squares_3.hpp"

#include <cmath>

namespace mpcns::operators
{
namespace
{

constexpr double kMinimumEdgeLength = 1.0e-14;

void add_edge_sample(FieldBlock &edge_current, FieldBlock &edge_vector,
                     int i, int j, int k,
                     double normal_matrix[3][3], double rhs[3], int &count)
{
    const double dx = edge_vector(i, j, k, 0);
    const double dy = edge_vector(i, j, k, 1);
    const double dz = edge_vector(i, j, k, 2);
    const double length = std::sqrt(dx * dx + dy * dy + dz * dz);
    if (length < kMinimumEdgeLength)
        return;

    const double tangent[3] = {dx / length, dy / length, dz / length};
    const double sample = edge_current(i, j, k, 0) / length;
    for (int row = 0; row < 3; ++row)
    {
        rhs[row] += sample * tangent[row];
        for (int column = 0; column < 3; ++column)
            normal_matrix[row][column] += tangent[row] * tangent[column];
    }
    ++count;
}

void write_ring(FieldBlock &cell_current, int i, int j,
                double normal_matrix[3][3], double rhs[3], int count,
                int k_lo, int k_hi)
{
    if (count <= 0)
        return;

    double current[3] = {0.0, 0.0, 0.0};
    if (!constraints::solve_regularized_least_squares_3(
            normal_matrix, rhs, current))
        return;

    for (int k = k_lo; k < k_hi; ++k)
    {
        cell_current(i, j, k, 0) = current[0];
        cell_current(i, j, k, 1) = current[1];
        cell_current(i, j, k, 2) = current[2];
    }
}

}

void PoleCurrentCellConstraint::apply(
    const std::vector<PoleCurrentCellPatch> &patches)
{
    for (const PoleCurrentCellPatch &patch : patches)
    {
        const int normal_axis = std::abs(patch.signed_normal_axis);
        if (normal_axis != 1 && normal_axis != 2)
            continue;

        if (!patch.cell_current || !patch.edge_current_xi ||
            !patch.edge_current_eta || !patch.edge_current_zeta ||
            !patch.edge_vector_xi || !patch.edge_vector_eta ||
            !patch.edge_vector_zeta)
            continue;

        FieldBlock &cell_current = *patch.cell_current;
        FieldBlock &current_xi = *patch.edge_current_xi;
        FieldBlock &current_eta = *patch.edge_current_eta;
        FieldBlock &current_zeta = *patch.edge_current_zeta;
        FieldBlock &vector_xi = *patch.edge_vector_xi;
        FieldBlock &vector_eta = *patch.edge_vector_eta;
        FieldBlock &vector_zeta = *patch.edge_vector_zeta;
        if (!cell_current.is_allocated() || !current_xi.is_allocated() ||
            !current_eta.is_allocated() || !current_zeta.is_allocated() ||
            !vector_xi.is_allocated() || !vector_eta.is_allocated() ||
            !vector_zeta.is_allocated())
            continue;

        const bool high_side = patch.signed_normal_axis > 0;
        if (normal_axis == 1)
        {
            const int cell_i = high_side ? patch.node_lo_i - 1 : patch.node_lo_i;
            const int tangent_i = high_side ? cell_i : cell_i + 1;
            for (int j = patch.node_lo_j; j < patch.node_hi_j - 1; ++j)
            {
                double normal_matrix[3][3] = {};
                double rhs[3] = {};
                int count = 0;
                for (int k = patch.node_lo_k; k < patch.node_hi_k - 1; ++k)
                {
                    add_edge_sample(current_xi, vector_xi, cell_i, j, k, normal_matrix, rhs, count);
                    add_edge_sample(current_xi, vector_xi, cell_i, j + 1, k, normal_matrix, rhs, count);
                    add_edge_sample(current_xi, vector_xi, cell_i, j, k + 1, normal_matrix, rhs, count);
                    add_edge_sample(current_xi, vector_xi, cell_i, j + 1, k + 1, normal_matrix, rhs, count);

                    add_edge_sample(current_eta, vector_eta, tangent_i, j, k, normal_matrix, rhs, count);
                    add_edge_sample(current_eta, vector_eta, tangent_i, j, k + 1, normal_matrix, rhs, count);
                    add_edge_sample(current_zeta, vector_zeta, tangent_i, j, k, normal_matrix, rhs, count);
                    add_edge_sample(current_zeta, vector_zeta, tangent_i, j + 1, k, normal_matrix, rhs, count);
                }
                write_ring(cell_current, cell_i, j, normal_matrix, rhs, count,
                           patch.node_lo_k, patch.node_hi_k - 1);
            }
        }
        else
        {
            const int cell_j = high_side ? patch.node_lo_j - 1 : patch.node_lo_j;
            const int tangent_j = high_side ? cell_j : cell_j + 1;
            for (int i = patch.node_lo_i; i < patch.node_hi_i - 1; ++i)
            {
                double normal_matrix[3][3] = {};
                double rhs[3] = {};
                int count = 0;
                for (int k = patch.node_lo_k; k < patch.node_hi_k - 1; ++k)
                {
                    add_edge_sample(current_eta, vector_eta, i, cell_j, k, normal_matrix, rhs, count);
                    add_edge_sample(current_eta, vector_eta, i + 1, cell_j, k, normal_matrix, rhs, count);
                    add_edge_sample(current_eta, vector_eta, i, cell_j, k + 1, normal_matrix, rhs, count);
                    add_edge_sample(current_eta, vector_eta, i + 1, cell_j, k + 1, normal_matrix, rhs, count);

                    add_edge_sample(current_xi, vector_xi, i, tangent_j, k, normal_matrix, rhs, count);
                    add_edge_sample(current_xi, vector_xi, i, tangent_j, k + 1, normal_matrix, rhs, count);
                    add_edge_sample(current_zeta, vector_zeta, i, tangent_j, k, normal_matrix, rhs, count);
                    add_edge_sample(current_zeta, vector_zeta, i + 1, tangent_j, k, normal_matrix, rhs, count);
                }
                write_ring(cell_current, i, cell_j, normal_matrix, rhs, count,
                           patch.node_lo_k, patch.node_hi_k - 1);
            }
        }
    }
}

} // namespace mpcns::operators
