#pragma once

class FieldBlock;

namespace CTOperators
{
// Oriented edge circulation E·dr -> oriented face flux curl(E)·dn.
// include_ghost requires edge buffers covering the complete face boundary,
// including ghost cells (equal edge/face ghost widths satisfy this).
void CurlEdgeToFace(int block,
                    FieldBlock &edge_xi, FieldBlock &edge_eta, FieldBlock &edge_zeta,
                    FieldBlock &face_xi, FieldBlock &face_eta, FieldBlock &face_zeta,
                    double multiplier, bool include_ghost = false);

// Hodge-weighted oriented face flux -> adjoint edge circulation.
void CurlAdjFaceToEdge(int block,
                       FieldBlock &face_xi, FieldBlock &face_eta, FieldBlock &face_zeta,
                       FieldBlock &beta_xi, FieldBlock &beta_eta, FieldBlock &beta_zeta,
                       FieldBlock &edge_xi, FieldBlock &edge_eta, FieldBlock &edge_zeta,
                       double multiplier);
}
