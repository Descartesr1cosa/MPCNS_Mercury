#pragma once

#include "mpcns/compat/legacy_cochain_view.hpp"
#include "mpcns/state/cochain_view.hpp"

namespace mpcns::operators
{

class WeightedCellReconstruction
{
public:
    static bool face_2form_to_cell_vector(
        const state::PrimalFaceFormView &flux,
        FieldBlock &weights,
        FieldBlock &cell_vector);

    static bool sum_face_2forms_to_cell_vector(
        const state::PrimalFaceFormView &primary_flux,
        const state::PrimalFaceFormView &additional_flux,
        FieldBlock &weights,
        FieldBlock &cell_vector);

    static bool edge_1form_to_cell_vector(
        const state::PrimalEdgeFormView &edge,
        FieldBlock &weights,
        FieldBlock &cell_vector);
};

}
