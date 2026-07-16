#include "operators/CTOperators.h"

namespace CTOperators
{
    // edge E (scalar, E·dr) -> face F (scalar,  F·dn) = curl E * multiper
    void CurlEdgeToFace(int iblk,
                        FieldBlock &Edge_xi, FieldBlock &Edge_eta, FieldBlock &Edge_zeta,
                        FieldBlock &Face_xi, FieldBlock &Face_eta, FieldBlock &Face_zeta, double multiper)
    {
        // (void)iblk; // 当前实现不依赖 block id（保留接口一致性）

        // -------- Face_xi = multiper * [ (d/dzeta)Edge_eta - (d/deta)Edge_zeta ] 的离散形式
        {
            Int3 sub = Face_xi.inner_lo();
            Int3 sup = Face_xi.inner_hi();
            for (int i = sub.i; i < sup.i; ++i)
                for (int j = sub.j; j < sup.j; ++j)
                    for (int k = sub.k; k < sup.k; ++k)
                    {
                        const double curl_xi =
                            (Edge_eta(i, j, k, 0) - Edge_eta(i, j, k + 1, 0)) +
                            (Edge_zeta(i, j + 1, k, 0) - Edge_zeta(i, j, k, 0));
                        Face_xi(i, j, k, 0) = multiper * curl_xi;
                    }
        }

        // -------- Face_eta = multiper * [ (d/dxi)Edge_zeta - (d/dzeta)Edge_xi ] 的离散形式
        {
            Int3 sub = Face_eta.inner_lo();
            Int3 sup = Face_eta.inner_hi();
            for (int i = sub.i; i < sup.i; ++i)
                for (int j = sub.j; j < sup.j; ++j)
                    for (int k = sub.k; k < sup.k; ++k)
                    {
                        const double curl_eta =
                            (Edge_xi(i, j, k + 1, 0) - Edge_xi(i, j, k, 0)) +
                            (Edge_zeta(i, j, k, 0) - Edge_zeta(i + 1, j, k, 0));
                        Face_eta(i, j, k, 0) = multiper * curl_eta;
                    }
        }

        // -------- Face_zeta = multiper * [ (d/deta)Edge_xi - (d/dxi)Edge_eta ] 的离散形式
        {
            Int3 sub = Face_zeta.inner_lo();
            Int3 sup = Face_zeta.inner_hi();
            for (int i = sub.i; i < sup.i; ++i)
                for (int j = sub.j; j < sup.j; ++j)
                    for (int k = sub.k; k < sup.k; ++k)
                    {
                        const double curl_zeta =
                            (Edge_xi(i, j, k, 0) - Edge_xi(i, j + 1, k, 0)) +
                            (Edge_eta(i + 1, j, k, 0) - Edge_eta(i, j, k, 0));
                        Face_zeta(i, j, k, 0) = multiper * curl_zeta;
                    }
        }
    }

    // face F (scalar) -> edge E (scalars per edge) = curl F * multiper, in fact circulation , and Hodge *:2form -> 1form is required
    void CurlAdjFaceToEdge(int iblk,
                           FieldBlock &Face_xi, FieldBlock &Face_eta, FieldBlock &Face_zeta,
                           FieldBlock &beta_xi, FieldBlock &beta_eta, FieldBlock &beta_zeta,
                           FieldBlock &Edge_xi, FieldBlock &Edge_eta, FieldBlock &Edge_zeta, double multiper)
    {
        // (void)iblk;

        // Edge_xi
        {
            Int3 sub = Edge_xi.inner_lo();
            Int3 sup = Edge_xi.inner_hi();
            for (int i = sub.i; i < sup.i; ++i)
                for (int j = sub.j; j < sup.j; ++j)
                    for (int k = sub.k; k < sup.k; ++k)
                    {
                        const double val =
                            (beta_eta(i, j, k - 1, 0) * Face_eta(i, j, k - 1, 0) - beta_eta(i, j, k, 0) * Face_eta(i, j, k, 0)) +
                            (beta_zeta(i, j, k, 0) * Face_zeta(i, j, k, 0) - beta_zeta(i, j - 1, k, 0) * Face_zeta(i, j - 1, k, 0));
                        Edge_xi(i, j, k, 0) = multiper * val;
                    }
        }

        // Edge_eta
        {
            Int3 sub = Edge_eta.inner_lo();
            Int3 sup = Edge_eta.inner_hi();
            for (int i = sub.i; i < sup.i; ++i)
                for (int j = sub.j; j < sup.j; ++j)
                    for (int k = sub.k; k < sup.k; ++k)
                    {
                        const double val =
                            (beta_xi(i, j, k, 0) * Face_xi(i, j, k, 0) - beta_xi(i, j, k - 1, 0) * Face_xi(i, j, k - 1, 0)) +
                            (beta_zeta(i - 1, j, k, 0) * Face_zeta(i - 1, j, k, 0) - beta_zeta(i, j, k, 0) * Face_zeta(i, j, k, 0));
                        Edge_eta(i, j, k, 0) = multiper * val;
                    }
        }

        // Edge_zeta
        {
            Int3 sub = Edge_zeta.inner_lo();
            Int3 sup = Edge_zeta.inner_hi();
            for (int i = sub.i; i < sup.i; ++i)
                for (int j = sub.j; j < sup.j; ++j)
                    for (int k = sub.k; k < sup.k; ++k)
                    {
                        const double val =
                            (beta_xi(i, j - 1, k, 0) * Face_xi(i, j - 1, k, 0) - beta_xi(i, j, k, 0) * Face_xi(i, j, k, 0)) +
                            (beta_eta(i, j, k, 0) * Face_eta(i, j, k, 0) - beta_eta(i - 1, j, k, 0) * Face_eta(i - 1, j, k, 0));
                        Edge_zeta(i, j, k, 0) = multiper * val;
                    }
        }
    }
}
