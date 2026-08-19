#pragma once

class FieldBlock;

namespace CTOperators
{
// Oriented edge circulation E·dr -> oriented face flux curl(E)·dn.
void CurlEdgeToFace(int block,
                    FieldBlock &edge_xi, FieldBlock &edge_eta, FieldBlock &edge_zeta,
                    FieldBlock &face_xi, FieldBlock &face_eta, FieldBlock &face_zeta,
                    double multiplier);

// Hodge-weighted oriented face flux -> adjoint edge circulation.
void CurlAdjFaceToEdge(int block,
                       FieldBlock &face_xi, FieldBlock &face_eta, FieldBlock &face_zeta,
                       FieldBlock &beta_xi, FieldBlock &beta_eta, FieldBlock &beta_zeta,
                       FieldBlock &edge_xi, FieldBlock &edge_eta, FieldBlock &edge_zeta,
                       double multiplier);
}
