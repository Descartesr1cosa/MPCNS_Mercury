#pragma once

namespace mpcns::constraints
{

bool solve_regularized_least_squares_3(
    double normal_matrix[3][3],
    double rhs[3],
    double solution[3]);

}
