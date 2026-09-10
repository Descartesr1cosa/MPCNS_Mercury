#include "8_dec/CTOperators.h"

#include "3_field/FieldBlock.h"

namespace CTOperators
{
void CurlEdgeToFace(int block,
                    FieldBlock &edge_xi, FieldBlock &edge_eta, FieldBlock &edge_zeta,
                    FieldBlock &face_xi, FieldBlock &face_eta, FieldBlock &face_zeta,
                    double multiplier, bool include_ghost)
{
    (void)block;
    {
        const Int3 lo=include_ghost?face_xi.get_lo():face_xi.inner_lo();
        const Int3 hi=include_ghost?face_xi.get_hi():face_xi.inner_hi();
        for(int i=lo.i;i<hi.i;++i) for(int j=lo.j;j<hi.j;++j) for(int k=lo.k;k<hi.k;++k)
        {
            const double curl=(edge_eta(i,j,k,0)-edge_eta(i,j,k+1,0))+
                              (edge_zeta(i,j+1,k,0)-edge_zeta(i,j,k,0));
            face_xi(i,j,k,0)=multiplier*curl;
        }
    }
    {
        const Int3 lo=include_ghost?face_eta.get_lo():face_eta.inner_lo();
        const Int3 hi=include_ghost?face_eta.get_hi():face_eta.inner_hi();
        for(int i=lo.i;i<hi.i;++i) for(int j=lo.j;j<hi.j;++j) for(int k=lo.k;k<hi.k;++k)
        {
            const double curl=(edge_xi(i,j,k+1,0)-edge_xi(i,j,k,0))+
                              (edge_zeta(i,j,k,0)-edge_zeta(i+1,j,k,0));
            face_eta(i,j,k,0)=multiplier*curl;
        }
    }
    {
        const Int3 lo=include_ghost?face_zeta.get_lo():face_zeta.inner_lo();
        const Int3 hi=include_ghost?face_zeta.get_hi():face_zeta.inner_hi();
        for(int i=lo.i;i<hi.i;++i) for(int j=lo.j;j<hi.j;++j) for(int k=lo.k;k<hi.k;++k)
        {
            const double curl=(edge_xi(i,j,k,0)-edge_xi(i,j+1,k,0))+
                              (edge_eta(i+1,j,k,0)-edge_eta(i,j,k,0));
            face_zeta(i,j,k,0)=multiplier*curl;
        }
    }
}

void CurlAdjFaceToEdge(int block,
                       FieldBlock &face_xi, FieldBlock &face_eta, FieldBlock &face_zeta,
                       FieldBlock &beta_xi, FieldBlock &beta_eta, FieldBlock &beta_zeta,
                       FieldBlock &edge_xi, FieldBlock &edge_eta, FieldBlock &edge_zeta,
                       double multiplier)
{
    (void)block;
    {
        const Int3 lo=edge_xi.inner_lo(),hi=edge_xi.inner_hi();
        for(int i=lo.i;i<hi.i;++i) for(int j=lo.j;j<hi.j;++j) for(int k=lo.k;k<hi.k;++k)
        {
            const double value=
                (beta_eta(i,j,k-1,0)*face_eta(i,j,k-1,0)-beta_eta(i,j,k,0)*face_eta(i,j,k,0))+
                (beta_zeta(i,j,k,0)*face_zeta(i,j,k,0)-beta_zeta(i,j-1,k,0)*face_zeta(i,j-1,k,0));
            edge_xi(i,j,k,0)=multiplier*value;
        }
    }
    {
        const Int3 lo=edge_eta.inner_lo(),hi=edge_eta.inner_hi();
        for(int i=lo.i;i<hi.i;++i) for(int j=lo.j;j<hi.j;++j) for(int k=lo.k;k<hi.k;++k)
        {
            const double value=
                (beta_xi(i,j,k,0)*face_xi(i,j,k,0)-beta_xi(i,j,k-1,0)*face_xi(i,j,k-1,0))+
                (beta_zeta(i-1,j,k,0)*face_zeta(i-1,j,k,0)-beta_zeta(i,j,k,0)*face_zeta(i,j,k,0));
            edge_eta(i,j,k,0)=multiplier*value;
        }
    }
    {
        const Int3 lo=edge_zeta.inner_lo(),hi=edge_zeta.inner_hi();
        for(int i=lo.i;i<hi.i;++i) for(int j=lo.j;j<hi.j;++j) for(int k=lo.k;k<hi.k;++k)
        {
            const double value=
                (beta_xi(i,j-1,k,0)*face_xi(i,j-1,k,0)-beta_xi(i,j,k,0)*face_xi(i,j,k,0))+
                (beta_eta(i,j,k,0)*face_eta(i,j,k,0)-beta_eta(i-1,j,k,0)*face_eta(i-1,j,k,0));
            edge_zeta(i,j,k,0)=multiplier*value;
        }
    }
}
}
