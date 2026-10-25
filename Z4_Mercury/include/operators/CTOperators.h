#pragma once

#include "3_field/FieldBlock.h"

namespace CTOperators
{
    // edge E (scalar, E·dr) -> face F (scalar,  F·dn) = curl E * multiper
    void CurlEdgeToFace(int iblk,
                        FieldBlock &Edge_xi, FieldBlock &Edge_eta, FieldBlock &Edge_zeta,
                        FieldBlock &Face_xi, FieldBlock &Face_eta, FieldBlock &Face_zeta, double multiper);

    // face F (scalar) -> edge E (scalars per edge) = curl F * multiper
    void CurlAdjFaceToEdge(int iblk,
                           FieldBlock &Face_xi, FieldBlock &Face_eta, FieldBlock &Face_zeta,
                           FieldBlock &Hodge_star_2form_to_1form_face_xi, FieldBlock &Hodge_star_2form_to_1form_face_eta, FieldBlock &Hodge_star_2form_to_1form_face_zeta,
                           FieldBlock &Edge_xi, FieldBlock &Edge_eta, FieldBlock &Edge_zeta, double multiper);
}
