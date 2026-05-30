#pragma once

#include <string>

#include "3_field/Field.h"
#include "3_field/FieldBlock.h"
#include "8_dec/DecTypes.h"

namespace DEC
{
    void d0_node_to_edge(FieldBlock &node,
                         FieldBlock &edge_xi,
                         FieldBlock &edge_eta,
                         FieldBlock &edge_zeta,
                         int node_comp = 0,
                         int edge_comp = 0);

    void d0_node_to_edge(Field &fields,
                         const std::string &node_name,
                         const EdgeFormNames &edge_names,
                         int node_comp = 0,
                         int edge_comp = 0);

    void d1_edge_to_face(FieldBlock &edge_xi,
                         FieldBlock &edge_eta,
                         FieldBlock &edge_zeta,
                         FieldBlock &face_xi,
                         FieldBlock &face_eta,
                         FieldBlock &face_zeta,
                         int edge_comp = 0,
                         int face_comp = 0);

    void d1_edge_to_face(Field &fields,
                         const EdgeFormNames &edge_names,
                         const FaceFormNames &face_names,
                         int edge_comp = 0,
                         int face_comp = 0);

    void d2_face_to_cell(FieldBlock &face_xi,
                         FieldBlock &face_eta,
                         FieldBlock &face_zeta,
                         FieldBlock &cell,
                         int face_comp = 0,
                         int cell_comp = 0);

    void d2_face_to_cell(Field &fields,
                         const FaceFormNames &face_names,
                         const std::string &cell_name,
                         int face_comp = 0,
                         int cell_comp = 0);

    void faraday_update_face_2form(FieldBlock &edge_xi,
                                   FieldBlock &edge_eta,
                                   FieldBlock &edge_zeta,
                                   FieldBlock &face_xi,
                                   FieldBlock &face_eta,
                                   FieldBlock &face_zeta,
                                   double dt,
                                   int edge_comp = 0,
                                   int face_comp = 0);

    void faraday_update_face_2form(Field &fields,
                                   const EdgeFormNames &edge_names,
                                   const FaceFormNames &face_names,
                                   double dt,
                                   int edge_comp = 0,
                                   int face_comp = 0);
}
