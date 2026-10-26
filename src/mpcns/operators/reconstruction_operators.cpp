#include "mpcns/operators/reconstruction_operators.hpp"

namespace mpcns::operators
{
namespace
{
bool face_inputs_allocated(const state::PrimalFaceFormView &flux,
                           const FieldBlock &weights,
                           const FieldBlock &cell_vector)
{
    return flux.xi.is_allocated() && flux.eta.is_allocated() &&
           flux.zeta.is_allocated() && weights.is_allocated() &&
           cell_vector.is_allocated();
}

bool edge_inputs_allocated(const state::PrimalEdgeFormView &edge,
                           const FieldBlock &weights,
                           const FieldBlock &cell_vector)
{
    return edge.xi.is_allocated() && edge.eta.is_allocated() &&
           edge.zeta.is_allocated() && weights.is_allocated() &&
           cell_vector.is_allocated();
}

void write_weighted_face_vector(const double phi[6],
                                const FieldBlock &weights,
                                FieldBlock &cell_vector,
                                int i, int j, int k)
{
    double value_x = 0.0;
    double value_y = 0.0;
    double value_z = 0.0;
    for (int n = 0; n < 6; ++n)
    {
        value_x += weights(i, j, k, n) * phi[n];
        value_y += weights(i, j, k, 6 + n) * phi[n];
        value_z += weights(i, j, k, 12 + n) * phi[n];
    }
    cell_vector(i, j, k, 0) = value_x;
    cell_vector(i, j, k, 1) = value_y;
    cell_vector(i, j, k, 2) = value_z;
}

}

bool WeightedCellReconstruction::face_2form_to_cell_vector(
    const state::PrimalFaceFormView &flux,
    FieldBlock &weights,
    FieldBlock &cell_vector)
{
    if (!face_inputs_allocated(flux, weights, cell_vector))
        return false;

    const Int3 lo = cell_vector.inner_lo();
    const Int3 hi = cell_vector.inner_hi();
    for (int i = lo.i; i < hi.i; ++i)
        for (int j = lo.j; j < hi.j; ++j)
            for (int k = lo.k; k < hi.k; ++k)
            {
                const double phi[6] = {
                    -flux.xi(i, j, k, flux.component),
                    flux.xi(i + 1, j, k, flux.component),
                    -flux.eta(i, j, k, flux.component),
                    flux.eta(i, j + 1, k, flux.component),
                    -flux.zeta(i, j, k, flux.component),
                    flux.zeta(i, j, k + 1, flux.component)};
                write_weighted_face_vector(phi, weights, cell_vector, i, j, k);
            }
    return true;
}

bool WeightedCellReconstruction::sum_face_2forms_to_cell_vector(
    const state::PrimalFaceFormView &primary_flux,
    const state::PrimalFaceFormView &additional_flux,
    FieldBlock &weights,
    FieldBlock &cell_vector)
{
    if (!face_inputs_allocated(primary_flux, weights, cell_vector) ||
        !face_inputs_allocated(additional_flux, weights, cell_vector))
        return false;

    const Int3 lo = cell_vector.inner_lo();
    const Int3 hi = cell_vector.inner_hi();
    for (int i = lo.i; i < hi.i; ++i)
        for (int j = lo.j; j < hi.j; ++j)
            for (int k = lo.k; k < hi.k; ++k)
            {
                const double phi[6] = {
                    -primary_flux.xi(i, j, k, primary_flux.component) - additional_flux.xi(i, j, k, additional_flux.component),
                    primary_flux.xi(i + 1, j, k, primary_flux.component) + additional_flux.xi(i + 1, j, k, additional_flux.component),
                    -primary_flux.eta(i, j, k, primary_flux.component) - additional_flux.eta(i, j, k, additional_flux.component),
                    primary_flux.eta(i, j + 1, k, primary_flux.component) + additional_flux.eta(i, j + 1, k, additional_flux.component),
                    -primary_flux.zeta(i, j, k, primary_flux.component) - additional_flux.zeta(i, j, k, additional_flux.component),
                    primary_flux.zeta(i, j, k + 1, primary_flux.component) + additional_flux.zeta(i, j, k + 1, additional_flux.component)};
                write_weighted_face_vector(phi, weights, cell_vector, i, j, k);
            }
    return true;
}

bool WeightedCellReconstruction::edge_1form_to_cell_vector(
    const state::PrimalEdgeFormView &edge,
    FieldBlock &weights,
    FieldBlock &cell_vector)
{
    if (!edge_inputs_allocated(edge, weights, cell_vector))
        return false;

    const Int3 lo = cell_vector.inner_lo();
    const Int3 hi = cell_vector.inner_hi();
    for (int i = lo.i; i < hi.i; ++i)
        for (int j = lo.j; j < hi.j; ++j)
            for (int k = lo.k; k < hi.k; ++k)
            {
                const double samples[12] = {
                    edge.xi(i, j, k, edge.component),
                    edge.xi(i, j + 1, k, edge.component),
                    edge.xi(i, j, k + 1, edge.component),
                    edge.xi(i, j + 1, k + 1, edge.component),
                    edge.eta(i, j, k, edge.component),
                    edge.eta(i + 1, j, k, edge.component),
                    edge.eta(i, j, k + 1, edge.component),
                    edge.eta(i + 1, j, k + 1, edge.component),
                    edge.zeta(i, j, k, edge.component),
                    edge.zeta(i + 1, j, k, edge.component),
                    edge.zeta(i, j + 1, k, edge.component),
                    edge.zeta(i + 1, j + 1, k, edge.component)};

                double value_x = 0.0;
                double value_y = 0.0;
                double value_z = 0.0;
                for (int n = 0; n < 12; ++n)
                {
                    value_x += weights(i, j, k, n) * samples[n];
                    value_y += weights(i, j, k, 12 + n) * samples[n];
                    value_z += weights(i, j, k, 24 + n) * samples[n];
                }
                cell_vector(i, j, k, 0) = value_x;
                cell_vector(i, j, k, 1) = value_y;
                cell_vector(i, j, k, 2) = value_z;
            }
    return true;
}

}
