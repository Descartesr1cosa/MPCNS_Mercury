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

std::vector<HaloRegion> Halo::debug_halo_regions(StaggerLocation location,
                                                 int nghost,
                                                 HaloLevel stage) const
{
    const PatternKey key = {location, nghost};
    std::vector<HaloRegion> regions;

    auto append = [&regions](const std::map<PatternKey, HaloPattern> &patterns,
                             const PatternKey &pattern_key)
    {
        auto it = patterns.find(pattern_key);
        if (it == patterns.end())
            return;
        regions.insert(regions.end(), it->second.regions.begin(), it->second.regions.end());
    };

    if (stage == HaloLevel::FaceOnly)
    {
        append(inner_patterns_, key);
        append(parallel_patterns_, key);
    }
    else if (stage == HaloLevel::Edge)
    {
        append(inner_edge_patterns_, key);
        append(parallel_edge_patterns_recv, key);
    }
    else if (stage == HaloLevel::Vertex)
    {
        append(inner_vertex_patterns_, key);
        append(parallel_vertex_patterns_recv, key);
    }

    return regions;
}

std::vector<HaloRegion> Halo::debug_halo_send_regions(StaggerLocation location,
                                                       int nghost,
                                                       HaloLevel stage) const
{
    const PatternKey key = {location, nghost};
    std::vector<HaloRegion> regions;

    auto append = [&regions](const std::map<PatternKey, HaloPattern> &patterns,
                             const PatternKey &pattern_key)
    {
        auto it = patterns.find(pattern_key);
        if (it == patterns.end())
            return;
        regions.insert(regions.end(), it->second.regions.begin(), it->second.regions.end());
    };

    if (stage == HaloLevel::FaceOnly)
    {
        append(parallel_patterns_, key);
    }
    else if (stage == HaloLevel::Edge)
    {
        append(parallel_edge_patterns_send, key);
    }
    else if (stage == HaloLevel::Vertex)
    {
        append(parallel_vertex_patterns_send, key);
    }

    return regions;
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
