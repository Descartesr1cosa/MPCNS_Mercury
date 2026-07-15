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
                           FieldBlock &Hodge_star_2form_to_1form_face_xi, FieldBlock &Hodge_star_2form_to_1form_face_eta, FieldBlock &Hodge_star_2form_to_1form_face_zeta,
                           FieldBlock &Edge_xi, FieldBlock &Edge_eta, FieldBlock &Edge_zeta, double multiper);

    // D^T applied after a consistent face Hodge has already produced the
    // dual face 1-form.  No diagonal beta is applied in this path.
    void CurlAdjMetricFaceToEdge(int iblk,
                                 FieldBlock &MetricFace_xi,
                                 FieldBlock &MetricFace_eta,
                                 FieldBlock &MetricFace_zeta,
                                 FieldBlock &Edge_xi,
                                 FieldBlock &Edge_eta,
                                 FieldBlock &Edge_zeta,
                                 double multiper);
}
