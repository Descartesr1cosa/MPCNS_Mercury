#pragma once

#include "3_field/FieldBlock.h"

#include <vector>

namespace mpcns::operators
{

// A storage-level description of one Pole patch.  The Mercury application
// adapts its topology records to this type; the constraint itself has no
// Grid, topology, Field catalogue, boundary, or MPI dependency.
struct PoleCurrentCellPatch
{
    int signed_normal_axis = 0;
    int node_lo_i = 0;
    int node_lo_j = 0;
    int node_lo_k = 0;
    int node_hi_i = 0;
    int node_hi_j = 0;
    int node_hi_k = 0;

    FieldBlock *cell_current = nullptr;
    FieldBlock *edge_current_xi = nullptr;
    FieldBlock *edge_current_eta = nullptr;
    FieldBlock *edge_current_zeta = nullptr;
    FieldBlock *edge_vector_xi = nullptr;
    FieldBlock *edge_vector_eta = nullptr;
    FieldBlock *edge_vector_zeta = nullptr;
};

// Enforces the legacy Pole current reconstruction.  Patch order, edge sample
// order, regularization, and cell overwrite order are part of the operation's
// numerical contract.
class PoleCurrentCellConstraint
{
public:
    static void apply(const std::vector<PoleCurrentCellPatch> &patches);
};

} // namespace mpcns::operators
