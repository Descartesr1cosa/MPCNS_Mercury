#include "8_dec/DecOps.h"

namespace DEC
{
    namespace
    {
        double d1_xi_at(FieldBlock &edge_eta,
                        FieldBlock &edge_zeta,
                        int i,
                        int j,
                        int k,
                        int edge_comp)
        {
            return (edge_eta(i, j, k, edge_comp) - edge_eta(i, j, k + 1, edge_comp)) +
                   (edge_zeta(i, j + 1, k, edge_comp) - edge_zeta(i, j, k, edge_comp));
        }

        double d1_eta_at(FieldBlock &edge_xi,
                         FieldBlock &edge_zeta,
                         int i,
                         int j,
                         int k,
                         int edge_comp)
        {
            return (edge_xi(i, j, k + 1, edge_comp) - edge_xi(i, j, k, edge_comp)) +
                   (edge_zeta(i, j, k, edge_comp) - edge_zeta(i + 1, j, k, edge_comp));
        }

        double d1_zeta_at(FieldBlock &edge_xi,
                          FieldBlock &edge_eta,
                          int i,
                          int j,
                          int k,
                          int edge_comp)
        {
            return (edge_xi(i, j, k, edge_comp) - edge_xi(i, j + 1, k, edge_comp)) +
                   (edge_eta(i + 1, j, k, edge_comp) - edge_eta(i, j, k, edge_comp));
        }
    }

    void d0_node_to_edge(FieldBlock &node,
                         FieldBlock &edge_xi,
                         FieldBlock &edge_eta,
                         FieldBlock &edge_zeta,
                         int node_comp,
                         int edge_comp)
    {
        {
            const Int3 lo = edge_xi.inner_lo();
            const Int3 hi = edge_xi.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        edge_xi(i, j, k, edge_comp) =
                            node(i + 1, j, k, node_comp) - node(i, j, k, node_comp);
        }

        {
            const Int3 lo = edge_eta.inner_lo();
            const Int3 hi = edge_eta.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        edge_eta(i, j, k, edge_comp) =
                            node(i, j + 1, k, node_comp) - node(i, j, k, node_comp);
        }

        {
            const Int3 lo = edge_zeta.inner_lo();
            const Int3 hi = edge_zeta.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        edge_zeta(i, j, k, edge_comp) =
                            node(i, j, k + 1, node_comp) - node(i, j, k, node_comp);
        }
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
            FieldBlock &edge_xi = fields.field(edge_names.xi, ib);
            FieldBlock &edge_eta = fields.field(edge_names.eta, ib);
            FieldBlock &edge_zeta = fields.field(edge_names.zeta, ib);
            if (!node.is_allocated() || !edge_xi.is_allocated() ||
                !edge_eta.is_allocated() || !edge_zeta.is_allocated())
                continue;

            d0_node_to_edge(node, edge_xi, edge_eta, edge_zeta, node_comp, edge_comp);
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
        {
            const Int3 lo = face_xi.inner_lo();
            const Int3 hi = face_xi.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        face_xi(i, j, k, face_comp) =
                            d1_xi_at(edge_eta, edge_zeta, i, j, k, edge_comp);
        }

        {
            const Int3 lo = face_eta.inner_lo();
            const Int3 hi = face_eta.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        face_eta(i, j, k, face_comp) =
                            d1_eta_at(edge_xi, edge_zeta, i, j, k, edge_comp);
        }

        {
            const Int3 lo = face_zeta.inner_lo();
            const Int3 hi = face_zeta.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        face_zeta(i, j, k, face_comp) =
                            d1_zeta_at(edge_xi, edge_eta, i, j, k, edge_comp);
        }
    }

    void d1_edge_to_face(Field &fields,
                         const EdgeFormNames &edge_names,
                         const FaceFormNames &face_names,
                         int edge_comp,
                         int face_comp)
    {
        for (int ib = 0; ib < fields.num_blocks(); ++ib)
        {
            FieldBlock &edge_xi = fields.field(edge_names.xi, ib);
            FieldBlock &edge_eta = fields.field(edge_names.eta, ib);
            FieldBlock &edge_zeta = fields.field(edge_names.zeta, ib);
            FieldBlock &face_xi = fields.field(face_names.xi, ib);
            FieldBlock &face_eta = fields.field(face_names.eta, ib);
            FieldBlock &face_zeta = fields.field(face_names.zeta, ib);
            if (!edge_xi.is_allocated() || !edge_eta.is_allocated() ||
                !edge_zeta.is_allocated() || !face_xi.is_allocated() ||
                !face_eta.is_allocated() || !face_zeta.is_allocated())
                continue;

            d1_edge_to_face(edge_xi, edge_eta, edge_zeta,
                            face_xi, face_eta, face_zeta,
                            edge_comp, face_comp);
        }
    }

    void d2_face_to_cell(FieldBlock &face_xi,
                         FieldBlock &face_eta,
                         FieldBlock &face_zeta,
                         FieldBlock &cell,
                         int face_comp,
                         int cell_comp)
    {
        const Int3 lo = cell.inner_lo();
        const Int3 hi = cell.inner_hi();
        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                    cell(i, j, k, cell_comp) =
                        (face_xi(i + 1, j, k, face_comp) - face_xi(i, j, k, face_comp)) +
                        (face_eta(i, j + 1, k, face_comp) - face_eta(i, j, k, face_comp)) +
                        (face_zeta(i, j, k + 1, face_comp) - face_zeta(i, j, k, face_comp));
    }

    void d2_face_to_cell(Field &fields,
                         const FaceFormNames &face_names,
                         const std::string &cell_name,
                         int face_comp,
                         int cell_comp)
    {
        for (int ib = 0; ib < fields.num_blocks(); ++ib)
        {
            FieldBlock &face_xi = fields.field(face_names.xi, ib);
            FieldBlock &face_eta = fields.field(face_names.eta, ib);
            FieldBlock &face_zeta = fields.field(face_names.zeta, ib);
            FieldBlock &cell = fields.field(cell_name, ib);
            if (!face_xi.is_allocated() || !face_eta.is_allocated() ||
                !face_zeta.is_allocated() || !cell.is_allocated())
                continue;

            d2_face_to_cell(face_xi, face_eta, face_zeta, cell, face_comp, cell_comp);
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
        {
            const Int3 lo = face_xi.inner_lo();
            const Int3 hi = face_xi.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        face_xi(i, j, k, face_comp) -=
                            dt * d1_xi_at(edge_eta, edge_zeta, i, j, k, edge_comp);
        }

        {
            const Int3 lo = face_eta.inner_lo();
            const Int3 hi = face_eta.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        face_eta(i, j, k, face_comp) -=
                            dt * d1_eta_at(edge_xi, edge_zeta, i, j, k, edge_comp);
        }

        {
            const Int3 lo = face_zeta.inner_lo();
            const Int3 hi = face_zeta.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        face_zeta(i, j, k, face_comp) -=
                            dt * d1_zeta_at(edge_xi, edge_eta, i, j, k, edge_comp);
        }
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
            FieldBlock &edge_xi = fields.field(edge_names.xi, ib);
            FieldBlock &edge_eta = fields.field(edge_names.eta, ib);
            FieldBlock &edge_zeta = fields.field(edge_names.zeta, ib);
            FieldBlock &face_xi = fields.field(face_names.xi, ib);
            FieldBlock &face_eta = fields.field(face_names.eta, ib);
            FieldBlock &face_zeta = fields.field(face_names.zeta, ib);
            if (!edge_xi.is_allocated() || !edge_eta.is_allocated() ||
                !edge_zeta.is_allocated() || !face_xi.is_allocated() ||
                !face_eta.is_allocated() || !face_zeta.is_allocated())
                continue;

            faraday_update_face_2form(edge_xi, edge_eta, edge_zeta,
                                      face_xi, face_eta, face_zeta,
                                      dt, edge_comp, face_comp);
        }
    }
}
