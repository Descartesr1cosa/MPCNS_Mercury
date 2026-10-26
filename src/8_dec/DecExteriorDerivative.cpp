#include "8_dec/DecOps.h"

#include "mpcns/discrete/incidence_operator.hpp"

namespace DEC
{
namespace
{

bool allocated(const FieldBlock &a, const FieldBlock &b, const FieldBlock &c)
{
    return a.is_allocated() && b.is_allocated() && c.is_allocated();
}

} // namespace

void d0_node_to_edge(FieldBlock &node,
                     FieldBlock &edge_xi,
                     FieldBlock &edge_eta,
                     FieldBlock &edge_zeta,
                     int node_comp,
                     int edge_comp)
{
    mpcns::discrete::apply_d0(
        mpcns::state::PrimalZeroFormView{node, node_comp},
        mpcns::state::PrimalEdgeFormView{
            edge_xi, edge_eta, edge_zeta, edge_comp});
}

void d0_node_to_edge(Field &fields,
                     const std::string &node_name,
                     const EdgeFormNames &edge_names,
                     int node_comp,
                     int edge_comp)
{
    for (int ib = 0; ib < fields.num_blocks(); ++ib)
    {
        FieldBlock &node = fields.field(node_name, ib);
        FieldBlock &xi = fields.field(edge_names.xi, ib);
        FieldBlock &eta = fields.field(edge_names.eta, ib);
        FieldBlock &zeta = fields.field(edge_names.zeta, ib);
        if (node.is_allocated() && allocated(xi, eta, zeta))
            d0_node_to_edge(node, xi, eta, zeta, node_comp, edge_comp);
    }
}

void d1_edge_to_face(FieldBlock &edge_xi,
                     FieldBlock &edge_eta,
                     FieldBlock &edge_zeta,
                     FieldBlock &face_xi,
                     FieldBlock &face_eta,
                     FieldBlock &face_zeta,
                     int edge_comp,
                     int face_comp)
{
    mpcns::discrete::apply_d1(
        mpcns::state::PrimalEdgeFormView{
            edge_xi, edge_eta, edge_zeta, edge_comp},
        mpcns::state::PrimalFaceFormView{
            face_xi, face_eta, face_zeta, face_comp});
}

void d1_edge_to_face(Field &fields,
                     const EdgeFormNames &edge_names,
                     const FaceFormNames &face_names,
                     int edge_comp,
                     int face_comp)
{
    for (int ib = 0; ib < fields.num_blocks(); ++ib)
    {
        FieldBlock &ex = fields.field(edge_names.xi, ib);
        FieldBlock &ey = fields.field(edge_names.eta, ib);
        FieldBlock &ez = fields.field(edge_names.zeta, ib);
        FieldBlock &fx = fields.field(face_names.xi, ib);
        FieldBlock &fy = fields.field(face_names.eta, ib);
        FieldBlock &fz = fields.field(face_names.zeta, ib);
        if (allocated(ex, ey, ez) && allocated(fx, fy, fz))
            d1_edge_to_face(ex, ey, ez, fx, fy, fz, edge_comp, face_comp);
    }
}

void d2_face_to_cell(FieldBlock &face_xi,
                     FieldBlock &face_eta,
                     FieldBlock &face_zeta,
                     FieldBlock &cell,
                     int face_comp,
                     int cell_comp)
{
    mpcns::discrete::apply_d2(
        mpcns::state::PrimalFaceFormView{
            face_xi, face_eta, face_zeta, face_comp},
        mpcns::state::PrimalThreeFormView{cell, cell_comp});
}

void d2_face_to_cell(Field &fields,
                     const FaceFormNames &face_names,
                     const std::string &cell_name,
                     int face_comp,
                     int cell_comp)
{
    for (int ib = 0; ib < fields.num_blocks(); ++ib)
    {
        FieldBlock &xi = fields.field(face_names.xi, ib);
        FieldBlock &eta = fields.field(face_names.eta, ib);
        FieldBlock &zeta = fields.field(face_names.zeta, ib);
        FieldBlock &cell = fields.field(cell_name, ib);
        if (allocated(xi, eta, zeta) && cell.is_allocated())
            d2_face_to_cell(xi, eta, zeta, cell, face_comp, cell_comp);
    }
}

void faraday_update_face_2form(FieldBlock &edge_xi,
                               FieldBlock &edge_eta,
                               FieldBlock &edge_zeta,
                               FieldBlock &face_xi,
                               FieldBlock &face_eta,
                               FieldBlock &face_zeta,
                               double dt,
                               int edge_comp,
                               int face_comp)
{
    mpcns::discrete::subtract_d1(
        mpcns::state::PrimalEdgeFormView{
            edge_xi, edge_eta, edge_zeta, edge_comp},
        mpcns::state::PrimalFaceFormView{
            face_xi, face_eta, face_zeta, face_comp},
        dt);
}

void faraday_update_face_2form(Field &fields,
                               const EdgeFormNames &edge_names,
                               const FaceFormNames &face_names,
                               double dt,
                               int edge_comp,
                               int face_comp)
{
    for (int ib = 0; ib < fields.num_blocks(); ++ib)
    {
        FieldBlock &ex = fields.field(edge_names.xi, ib);
        FieldBlock &ey = fields.field(edge_names.eta, ib);
        FieldBlock &ez = fields.field(edge_names.zeta, ib);
        FieldBlock &fx = fields.field(face_names.xi, ib);
        FieldBlock &fy = fields.field(face_names.eta, ib);
        FieldBlock &fz = fields.field(face_names.zeta, ib);
        if (allocated(ex, ey, ez) && allocated(fx, fy, fz))
            faraday_update_face_2form(ex, ey, ez, fx, fy, fz, dt,
                                      edge_comp, face_comp);
    }
}

} // namespace DEC
