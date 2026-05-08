#pragma once
#include "1_grid/1_MPCNS_Grid.h"
#include "2_topology/TopologyBuilder.h"
#include "3_field/Field.h"
#include "4_halo/Halo.h"

namespace CTOperators
{
    // edge E (scalar, E·dr) -> face F (scalar,  F·dn) = curl E * multiper
    void CurlEdgeToFace(int iblk,
                        FieldBlock &Edge_xi, FieldBlock &Edge_eta, FieldBlock &Edge_zeta,
                        FieldBlock &Face_xi, FieldBlock &Face_eta, FieldBlock &Face_zeta, double multiper);

    // face F (scalar) -> edge E (scalars per edge) = curl F * multiper
    void CurlAdjFaceToEdge(int iblk,
                           FieldBlock &Face_xi, FieldBlock &Face_eta, FieldBlock &Face_zeta,
                           FieldBlock &beta_xi, FieldBlock &beta_eta, FieldBlock &beta_zeta,
                           FieldBlock &Edge_xi, FieldBlock &Edge_eta, FieldBlock &Edge_zeta, double multiper);
}
