#include "8_dec/DecDiagnostics.h"

#include <algorithm>
#include <cmath>

namespace DEC
{
    namespace
    {
        double abs_curl_grad_xi(FieldBlock &node, int i, int j, int k, int comp)
        {
            const double edge_eta_0 = node(i, j + 1, k, comp) - node(i, j, k, comp);
            const double edge_eta_1 = node(i, j + 1, k + 1, comp) - node(i, j, k + 1, comp);
            const double edge_zeta_0 = node(i, j, k + 1, comp) - node(i, j, k, comp);
            const double edge_zeta_1 = node(i, j + 1, k + 1, comp) - node(i, j + 1, k, comp);
            return (edge_eta_0 - edge_eta_1) + (edge_zeta_1 - edge_zeta_0);
        }

        double abs_curl_grad_eta(FieldBlock &node, int i, int j, int k, int comp)
        {
            const double edge_xi_0 = node(i + 1, j, k, comp) - node(i, j, k, comp);
            const double edge_xi_1 = node(i + 1, j, k + 1, comp) - node(i, j, k + 1, comp);
            const double edge_zeta_0 = node(i, j, k + 1, comp) - node(i, j, k, comp);
            const double edge_zeta_1 = node(i + 1, j, k + 1, comp) - node(i + 1, j, k, comp);
            return (edge_xi_1 - edge_xi_0) + (edge_zeta_0 - edge_zeta_1);
        }

        double abs_curl_grad_zeta(FieldBlock &node, int i, int j, int k, int comp)
        {
            const double edge_xi_0 = node(i + 1, j, k, comp) - node(i, j, k, comp);
            const double edge_xi_1 = node(i + 1, j + 1, k, comp) - node(i, j + 1, k, comp);
            const double edge_eta_0 = node(i, j + 1, k, comp) - node(i, j, k, comp);
            const double edge_eta_1 = node(i + 1, j + 1, k, comp) - node(i + 1, j, k, comp);
            return (edge_xi_0 - edge_xi_1) + (edge_eta_1 - edge_eta_0);
        }

        double curl_xi_at(FieldBlock &edge_eta,
                          FieldBlock &edge_zeta,
                          int i,
                          int j,
                          int k,
                          int comp)
        {
            return (edge_eta(i, j, k, comp) - edge_eta(i, j, k + 1, comp)) +
                   (edge_zeta(i, j + 1, k, comp) - edge_zeta(i, j, k, comp));
        }

        double curl_eta_at(FieldBlock &edge_xi,
                           FieldBlock &edge_zeta,
                           int i,
                           int j,
                           int k,
                           int comp)
        {
            return (edge_xi(i, j, k + 1, comp) - edge_xi(i, j, k, comp)) +
                   (edge_zeta(i, j, k, comp) - edge_zeta(i + 1, j, k, comp));
        }

        double curl_zeta_at(FieldBlock &edge_xi,
                            FieldBlock &edge_eta,
                            int i,
                            int j,
                            int k,
                            int comp)
        {
            return (edge_xi(i, j, k, comp) - edge_xi(i, j + 1, k, comp)) +
                   (edge_eta(i + 1, j, k, comp) - edge_eta(i, j, k, comp));
        }
    }

    double max_abs_curl_grad(FieldBlock &node, int node_comp)
    {
        double max_abs = 0.0;
        const Int3 lo = node.inner_lo();
        const Int3 hi = node.inner_hi();

        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j - 1; ++j)
                for (int k = lo.k; k < hi.k - 1; ++k)
                    max_abs = std::max(max_abs, std::abs(abs_curl_grad_xi(node, i, j, k, node_comp)));

        for (int i = lo.i; i < hi.i - 1; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k - 1; ++k)
                    max_abs = std::max(max_abs, std::abs(abs_curl_grad_eta(node, i, j, k, node_comp)));

        for (int i = lo.i; i < hi.i - 1; ++i)
            for (int j = lo.j; j < hi.j - 1; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                    max_abs = std::max(max_abs, std::abs(abs_curl_grad_zeta(node, i, j, k, node_comp)));

        return max_abs;
    }

    double max_abs_curl_grad(Field &fields,
                             const std::string &node_name,
                             int node_comp)
    {
        double max_abs = 0.0;
        for (int ib = 0; ib < fields.num_blocks(); ++ib)
        {
            FieldBlock &node = fields.field(node_name, ib);
            if (!node.is_allocated())
                continue;
            max_abs = std::max(max_abs, max_abs_curl_grad(node, node_comp));
        }
        return max_abs;
    }

    double max_abs_div_curl(FieldBlock &edge_xi,
                            FieldBlock &edge_eta,
                            FieldBlock &edge_zeta,
                            int edge_comp)
    {
        double max_abs = 0.0;
        const Int3 lo = edge_xi.inner_lo();
        const Int3 hi = edge_xi.inner_hi();

        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j - 1; ++j)
                for (int k = lo.k; k < hi.k - 1; ++k)
                {
                    const double div_curl =
                        (curl_xi_at(edge_eta, edge_zeta, i + 1, j, k, edge_comp) -
                         curl_xi_at(edge_eta, edge_zeta, i, j, k, edge_comp)) +
                        (curl_eta_at(edge_xi, edge_zeta, i, j + 1, k, edge_comp) -
                         curl_eta_at(edge_xi, edge_zeta, i, j, k, edge_comp)) +
                        (curl_zeta_at(edge_xi, edge_eta, i, j, k + 1, edge_comp) -
                         curl_zeta_at(edge_xi, edge_eta, i, j, k, edge_comp));

                    max_abs = std::max(max_abs, std::abs(div_curl));
                }

        return max_abs;
    }

    double max_abs_div_curl(Field &fields,
                            const EdgeFormNames &edge_names,
                            int edge_comp)
    {
        double max_abs = 0.0;
        for (int ib = 0; ib < fields.num_blocks(); ++ib)
        {
            FieldBlock &edge_xi = fields.field(edge_names.xi, ib);
            FieldBlock &edge_eta = fields.field(edge_names.eta, ib);
            FieldBlock &edge_zeta = fields.field(edge_names.zeta, ib);
            if (!edge_xi.is_allocated() || !edge_eta.is_allocated() || !edge_zeta.is_allocated())
                continue;
            max_abs = std::max(max_abs, max_abs_div_curl(edge_xi, edge_eta, edge_zeta, edge_comp));
        }
        return max_abs;
    }
}
