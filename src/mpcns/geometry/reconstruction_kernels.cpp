#include "mpcns/geometry/reconstruction_kernels.hpp"

#include <cmath>

namespace mpcns::geometry
{
namespace
{

CartesianVector subtract(const CartesianVector &lhs,
                         const CartesianVector &rhs)
{
    return {lhs[0] - rhs[0],
            lhs[1] - rhs[1],
            lhs[2] - rhs[2]};
}

CartesianVector add(const CartesianVector &lhs,
                    const CartesianVector &rhs)
{
    return {lhs[0] + rhs[0],
            lhs[1] + rhs[1],
            lhs[2] + rhs[2]};
}

bool invert_matrix3(const Matrix3 &matrix, Matrix3 &inverse)
{
    const double a00 = matrix[0], a01 = matrix[1], a02 = matrix[2];
    const double a10 = matrix[3], a11 = matrix[4], a12 = matrix[5];
    const double a20 = matrix[6], a21 = matrix[7], a22 = matrix[8];

    const double c00 = a11 * a22 - a12 * a21;
    const double c01 = a02 * a21 - a01 * a22;
    const double c02 = a01 * a12 - a02 * a11;

    const double c10 = a12 * a20 - a10 * a22;
    const double c11 = a00 * a22 - a02 * a20;
    const double c12 = a02 * a10 - a00 * a12;

    const double c20 = a10 * a21 - a11 * a20;
    const double c21 = a01 * a20 - a00 * a21;
    const double c22 = a00 * a11 - a01 * a10;

    const double determinant = a00 * c00 + a01 * c10 + a02 * c20;
    if (std::fabs(determinant) < 1e-300)
        return false;

    const double inverse_determinant = 1.0 / determinant;
    inverse[0] = c00 * inverse_determinant;
    inverse[1] = c01 * inverse_determinant;
    inverse[2] = c02 * inverse_determinant;
    inverse[3] = c10 * inverse_determinant;
    inverse[4] = c11 * inverse_determinant;
    inverse[5] = c12 * inverse_determinant;
    inverse[6] = c20 * inverse_determinant;
    inverse[7] = c21 * inverse_determinant;
    inverse[8] = c22 * inverse_determinant;
    return true;
}

} // namespace

CellCovariantBasis cell_covariant_basis_from_corners(
    const std::array<CartesianVector, 8> &corners)
{
    const CartesianVector &r000 = corners[0];
    const CartesianVector &r100 = corners[1];
    const CartesianVector &r010 = corners[2];
    const CartesianVector &r110 = corners[3];
    const CartesianVector &r001 = corners[4];
    const CartesianVector &r101 = corners[5];
    const CartesianVector &r011 = corners[6];
    const CartesianVector &r111 = corners[7];

    const CartesianVector ex0 = subtract(r100, r000);
    const CartesianVector ex1 = subtract(r110, r010);
    const CartesianVector ex2 = subtract(r101, r001);
    const CartesianVector ex3 = subtract(r111, r011);
    const CartesianVector exs = add(add(ex0, ex1), add(ex2, ex3));

    const CartesianVector ey0 = subtract(r010, r000);
    const CartesianVector ey1 = subtract(r110, r100);
    const CartesianVector ey2 = subtract(r011, r001);
    const CartesianVector ey3 = subtract(r111, r101);
    const CartesianVector eys = add(add(ey0, ey1), add(ey2, ey3));

    const CartesianVector ez0 = subtract(r001, r000);
    const CartesianVector ez1 = subtract(r101, r100);
    const CartesianVector ez2 = subtract(r011, r010);
    const CartesianVector ez3 = subtract(r111, r110);
    const CartesianVector ezs = add(add(ez0, ez1), add(ez2, ez3));

    return CellCovariantBasis{
        {0.25 * exs[0], 0.25 * exs[1], 0.25 * exs[2]},
        {0.25 * eys[0], 0.25 * eys[1], 0.25 * eys[2]},
        {0.25 * ezs[0], 0.25 * ezs[1], 0.25 * ezs[2]}};
}

Matrix3 matrix3_from_row_major(const double values[9])
{
    return {values[0], values[1], values[2],
            values[3], values[4], values[5],
            values[6], values[7], values[8]};
}

Matrix3 regularized_pseudoinverse_transpose(const Matrix3 &matrix)
{
    Matrix3 gram{{0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0}};
    for (int row = 0; row < 3; ++row)
        for (int column = 0; column < 3; ++column)
        {
            double sum = 0.0;
            for (int inner = 0; inner < 3; ++inner)
                sum += matrix[row * 3 + inner] *
                       matrix[column * 3 + inner];
            gram[row * 3 + column] = sum;
        }

    const double trace = gram[0] + gram[4] + gram[8];
    double regularization = 1e-14 * (trace + 1e-300);

    Matrix3 inverse;
    bool is_invertible = false;
    for (int attempt = 0; attempt < 8; ++attempt)
    {
        const Matrix3 regularized{{
            gram[0] + regularization, gram[1], gram[2],
            gram[3], gram[4] + regularization, gram[5],
            gram[6], gram[7], gram[8] + regularization}};

        if (invert_matrix3(regularized, inverse))
        {
            is_invertible = true;
            break;
        }
        regularization *= 10.0;
    }

    Matrix3 pseudoinverse{{0.0, 0.0, 0.0,
                           0.0, 0.0, 0.0,
                           0.0, 0.0, 0.0}};
    if (!is_invertible)
        return pseudoinverse;

    for (int row = 0; row < 3; ++row)
        for (int column = 0; column < 3; ++column)
        {
            double sum = 0.0;
            for (int inner = 0; inner < 3; ++inner)
                sum += inverse[row * 3 + inner] *
                       matrix[inner * 3 + column];
            pseudoinverse[row * 3 + column] = sum;
        }
    return pseudoinverse;
}

} // namespace mpcns::geometry
