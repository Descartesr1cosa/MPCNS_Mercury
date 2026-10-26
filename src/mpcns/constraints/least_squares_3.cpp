#include "mpcns/constraints/least_squares_3.hpp"

#include <algorithm>
#include <cmath>

namespace mpcns::constraints
{

bool solve_regularized_least_squares_3(
    double normal_matrix[3][3],
    double rhs[3],
    double solution[3])
{
    const double trace = normal_matrix[0][0] + normal_matrix[1][1] + normal_matrix[2][2];
    const double regularization = 1.0e-12 * std::max(1.0, trace);
    normal_matrix[0][0] += regularization;
    normal_matrix[1][1] += regularization;
    normal_matrix[2][2] += regularization;

    double augmented[3][4] = {
        {normal_matrix[0][0], normal_matrix[0][1], normal_matrix[0][2], rhs[0]},
        {normal_matrix[1][0], normal_matrix[1][1], normal_matrix[1][2], rhs[1]},
        {normal_matrix[2][0], normal_matrix[2][1], normal_matrix[2][2], rhs[2]}};

    for (int column = 0; column < 3; ++column)
    {
        int pivot = column;
        double maximum = std::abs(augmented[column][column]);
        for (int row = column + 1; row < 3; ++row)
        {
            const double candidate = std::abs(augmented[row][column]);
            if (candidate > maximum)
            {
                maximum = candidate;
                pivot = row;
            }
        }

        if (maximum < 1.0e-300)
            return false;

        if (pivot != column)
            for (int q = column; q < 4; ++q)
                std::swap(augmented[column][q], augmented[pivot][q]);

        const double inverse_pivot = 1.0 / augmented[column][column];
        for (int q = column; q < 4; ++q)
            augmented[column][q] *= inverse_pivot;

        for (int row = 0; row < 3; ++row)
        {
            if (row == column)
                continue;
            const double factor = augmented[row][column];
            for (int q = column; q < 4; ++q)
                augmented[row][q] -= factor * augmented[column][q];
        }
    }

    solution[0] = augmented[0][3];
    solution[1] = augmented[1][3];
    solution[2] = augmented[2][3];
    return true;
}

}
