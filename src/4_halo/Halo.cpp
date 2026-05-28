#include "4_halo/Halo.h"

Halo::Halo(Field *field, TOPO::Topology *topo)
    : fld_(field), topo_(topo)
{
}

void Halo::set_topology_equiv(const TOPO::Topology *equiv)
{
    equiv_ = equiv;
}

const TOPO::Topology *Halo::topology_equiv() const
{
    return equiv_;
}

void Halo::data_trans_1DCorner(std::string &field_name)
{
    exchange_inner(field_name);
    exchange_parallel(field_name);
}

void Halo::data_trans_2DCorner(std::string &field_name)
{
    exchange_inner_edge(field_name);
    exchange_parallel_edge(field_name);
}

void Halo::data_trans_3DCorner(std::string &field_name)
{
    exchange_inner_vertex(field_name);
    exchange_parallel_vertex(field_name);
}

void Halo::coupling_trans_1DCorner(std::string &src, std::string &dst)
{
    coupling_inner_face(src, dst);
    coupling_parallel_face(src, dst);
}

void Halo::coupling_trans_1DCorner(std::string &src,
                                   std::string &dst,
                                   std::vector<int32_t> &field_cids)
{
    coupling_inner_face(src, dst, field_cids);
    coupling_parallel_face(src, dst, field_cids);
}

void Halo::coupling_trans_2DCorner(std::string &src, std::string &dst)
{
    coupling_inner_edge(src, dst);
    coupling_parallel_edge(src, dst);
}

void Halo::coupling_trans_2DCorner(std::string &src,
                                   std::string &dst,
                                   std::vector<int32_t> &field_cids)
{
    coupling_inner_edge(src, dst, field_cids);
    coupling_parallel_edge(src, dst, field_cids);
}

void Halo::coupling_trans_3DCorner(std::string &src, std::string &dst)
{
    coupling_inner_vertex(src, dst);
    coupling_parallel_vertex(src, dst);
}

void Halo::coupling_trans_3DCorner(std::string &src,
                                   std::string &dst,
                                   std::vector<int32_t> &field_cids)
{
    coupling_inner_vertex(src, dst, field_cids);
    coupling_parallel_vertex(src, dst, field_cids);
}
