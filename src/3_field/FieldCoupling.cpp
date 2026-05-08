#include "3_field/Field.h"
#include "2_topology/TopologyView.h"
#include "0_basic/Error.h"
#include "0_basic/LayoutTraits.h"
#include "1_grid/BlockTraits.h"

void Field::register_coupling_channel(const std::string &src,
                                      const std::string &dst,
                                      const std::string &tag,
                                      StaggerLocation location,
                                      FieldValueKind value_kind,
                                      int ncomp,
                                      int nghost,
                                      bool orientation_aware)
{
    if (has_field(tag))
    {
        const auto &desc = descriptor(field_id(tag));
        const bool same = desc.location == location &&
                          desc.value_kind == value_kind &&
                          desc.ncomp == ncomp &&
                          desc.nghost == nghost &&
                          desc.sync.orientation_aware == orientation_aware;
        if (!same)
        {
            ERROR::Abort("Field::register_coupling_channel: channel spec does not match field descriptor: " + tag);
        }
    }

    auto add_one = [&](const std::string &a, const std::string &b)
    {
        PairKey key{a, b};
        auto &pd = coupling_pairs_[key]; // 不存在会默认构造

        // 第一次创建时补 pair
        if (pd.pair.src.empty() && pd.pair.dst.empty())
        {
            pd.pair.src = a;
            pd.pair.dst = b;
        }

        // enforce: 同一 (a->b) 下 tag 唯一；若重复则必须完全一致
        for (const auto &ch : pd.channels)
        {
            if (ch.tag != tag)
                continue;

            const bool same =
                (ch.location == location) &&
                (ch.value_kind == value_kind) &&
                (ch.ncomp == ncomp) &&
                (ch.nghost == nghost) &&
                (ch.orientation_aware == orientation_aware);

            if (same)
                return; // 幂等：重复注册同样的 channel，直接忽略

            std::cout << "Fatal: coupling channel tag duplicated but spec differs: "
                      << a << " -> " << b << " tag=" << tag << "\n";
            std::exit(-1);
        }

        CouplingChannelSpec spec;
        spec.tag = tag;
        spec.location = location;
        spec.value_kind = value_kind;
        spec.ncomp = ncomp;
        spec.nghost = nghost;
        spec.orientation_aware = orientation_aware;

        pd.channels.push_back(spec);
    };

    add_one(src, dst);
}

void Field::register_coupling_channel(const std::string &src,
                                      const std::string &dst,
                                      const std::string &tag,
                                      StaggerLocation location,
                                      int ncomp,
                                      int nghost)
{
    register_coupling_channel(src, dst, tag, location, FieldValueKind::Scalar, ncomp, nghost, false);
}

void Field::register_coupling_channel(const std::string &src,
                                      const std::string &dst,
                                      const std::string &field_name)
{
    const auto &desc = descriptor(field_name);
    register_coupling_channel(src, dst, field_name,
                              desc.location,
                              desc.value_kind,
                              desc.ncomp,
                              desc.nghost,
                              desc.sync.orientation_aware);
}

void Field::register_declared_coupling_channels(const std::vector<PairKey> &directed_pairs)
{
    const std::vector<std::string> names = coupled_field_names();
    for (const auto &pair : directed_pairs)
    {
        for (const auto &field_name : names)
            register_coupling_channel(pair.first, pair.second, field_name);
    }
}

bool Field::has_coupling_pair(const std::string &src, const std::string &dst) const
{
    return coupling_pairs_.count(PairKey{src, dst}) > 0;
}

const CouplingPairDesc &Field::coupling_pair(const std::string &src, const std::string &dst) const
{
    PairKey key{src, dst};
    auto it = coupling_pairs_.find(key);
    if (it == coupling_pairs_.end())
    {
        std::cout << "Coupling pair not found: (" << src << " -> " << dst << ")\n";
        exit(-1);
    }
    return it->second;
}

CouplingBuffersForPair &Field::coupling_buffers(const std::string &src, const std::string &dst)
{
    PairKey k{src, dst};
    auto it = coupling_buffers_.find(k);
    if (it == coupling_buffers_.end())
        ERROR::Abort("Field::coupling_buffers: buffers not built for coupling pair");
    return it->second;
}

const std::map<Field::PairKey, CouplingPairDesc> &Field::coupling_pairs() const
{
    return coupling_pairs_;
}

void Field::build_coupling_buffers(const TOPO::Topology &topo, int dimension)
{
    // 0) 清空并为每个已注册 (src,dst) 建空壳
    coupling_buffers_.clear();

    const std::vector<TOPO_VIEW::FacePatchView> inner_faces = TOPO_VIEW::inner_faces(topo);
    const std::vector<TOPO_VIEW::FacePatchView> parallel_faces = TOPO_VIEW::parallel_faces(topo);
    const auto &inner_edges = TOPO_VIEW::edge_patches(topo, TOPO::PatchKind::Inner);
    const auto &parallel_edges = TOPO_VIEW::edge_patches(topo, TOPO::PatchKind::Parallel);
    const auto &inner_vertices = TOPO_VIEW::vertex_patches(topo, TOPO::PatchKind::Inner);
    const auto &parallel_vertices = TOPO_VIEW::vertex_patches(topo, TOPO::PatchKind::Parallel);

    // 1) 根据coupling_pairs_以及topo中patch的数量开辟[cid][ipatch]数组的空间
    for (const auto &kv : coupling_pairs_)
    {
        const PairKey &key = kv.first;
        const CouplingPairDesc &pd = kv.second;

        CouplingBuffersForPair bs;
        bs.desc = pd;

        const int nc = (int)pd.channels.size();

        bs.inner_face.assign(nc, std::vector<CouplingBufferBlock>(inner_faces.size()));
        bs.parallel_face.assign(nc, std::vector<CouplingBufferBlock>(parallel_faces.size()));

        bs.inner_edge.assign(nc, std::vector<CouplingBufferBlock>(inner_edges.size()));
        bs.parallel_edge.assign(nc, std::vector<CouplingBufferBlock>(parallel_edges.size()));

        bs.inner_vertex.assign(nc, std::vector<CouplingBufferBlock>(inner_vertices.size()));
        bs.parallel_vertex.assign(nc, std::vector<CouplingBufferBlock>(parallel_vertices.size()));

        coupling_buffers_[key] = std::move(bs);
    }

    // ---------- allocate one buffer block ----------
    auto alloc_block = [&](CouplingBufferBlock &cb,
                           TOPO::PatchKind kind,
                           int this_block,
                           const CouplingChannelSpec &ch,
                           const Box3 &box)
    {
        int Ni = box.hi.i - box.lo.i;
        int Nj = box.hi.j - box.lo.j;
        int Nk = box.hi.k - box.lo.k;
        if (Ni <= 0 || Nj <= 0 || Nk <= 0)
            return;

        cb.allocated = true;
        cb.kind = kind;
        cb.this_block = this_block;

        cb.tag = ch.tag;
        cb.location = ch.location;
        cb.value_kind = ch.value_kind;
        cb.ncomp = ch.ncomp;
        cb.nghost = ch.nghost;
        cb.orientation_aware = ch.orientation_aware;

        cb.box = box;
        cb.shift = Int3{-box.lo.i, -box.lo.j, -box.lo.k};

        cb.data.SetSize(Ni, Nj, Nk, 0, ch.ncomp);
    };

    // ---------- build for Interface (face) ----------
    auto build_face_list = [&](const auto &plist, TOPO::PatchKind kind,
                               auto CouplingBuffersForPair::*storage_member)
    {
        for (int ip = 0; ip < (int)plist.size(); ++ip)
        {
            const auto &p = plist[ip];
            if (!p.is_coupling)
                continue;

            PairKey key{p.nb_block_name, p.this_block_name};
            auto it = coupling_buffers_.find(key);
            if (it == coupling_buffers_.end())
                continue;

            CouplingBuffersForPair &bs = it->second;
            auto &storage = bs.*storage_member; // [cid][ip]

            const Block &blk = storage_.block(p.this_block);

            // 对于patch上的每一个channel开辟合适的空间
            for (int cid = 0; cid < (int)bs.desc.channels.size(); ++cid)
            {
                const auto &ch = bs.desc.channels[cid];
                Box3 box = LAYOUT::coupling_face_ghost_slab_from_cells(
                    GRID_TRAITS::cell_counts(blk),
                    ch.location,
                    p.this_box_node,
                    p.direction,
                    ch.nghost);
                alloc_block(storage[cid][ip], kind, p.this_block, ch, box);
            }
        }
    };

    build_face_list(inner_faces, TOPO::PatchKind::Inner, &CouplingBuffersForPair::inner_face);
    build_face_list(parallel_faces, TOPO::PatchKind::Parallel, &CouplingBuffersForPair::parallel_face);

    // ---------- build for Edge ----------
    if (dimension >= 2)
    {
        auto build_edge_list = [&](const auto &elist, TOPO::PatchKind kind,
                                   auto CouplingBuffersForPair::*storage_member)
        {
            for (int ie = 0; ie < (int)elist.size(); ++ie)
            {
                const auto &e = elist[ie];
                if (!e.is_coupling)
                    continue;

                PairKey key{e.nb_block_name, e.this_block_name};
                auto it = coupling_buffers_.find(key);
                if (it == coupling_buffers_.end())
                    continue;

                CouplingBuffersForPair &bs = it->second;
                auto &storage = bs.*storage_member;

                const Block &blk = storage_.block(e.this_block);

                for (int cid = 0; cid < (int)bs.desc.channels.size(); ++cid)
                {
                    const auto &ch = bs.desc.channels[cid];
                    Box3 box = LAYOUT::coupling_edge_ghost_slab_from_cells(
                        GRID_TRAITS::cell_counts(blk),
                        ch.location,
                        e.this_box_node,
                        e.dir1,
                        e.dir2,
                        ch.nghost);
                    alloc_block(storage[cid][ie], kind, e.this_block, ch, box);
                }
            }
        };

        build_edge_list(inner_edges, TOPO::PatchKind::Inner, &CouplingBuffersForPair::inner_edge);
        build_edge_list(parallel_edges, TOPO::PatchKind::Parallel, &CouplingBuffersForPair::parallel_edge);
    }

    // ---------- build for Vertex ----------
    if (dimension >= 3)
    {
        auto build_vertex_list = [&](const auto &vlist, TOPO::PatchKind kind,
                                     auto CouplingBuffersForPair::*storage_member)
        {
            for (int iv = 0; iv < (int)vlist.size(); ++iv)
            {
                const auto &v = vlist[iv];
                if (!v.is_coupling)
                    continue;

                PairKey key{v.nb_block_name, v.this_block_name};
                auto it = coupling_buffers_.find(key);
                if (it == coupling_buffers_.end())
                    continue;

                CouplingBuffersForPair &bs = it->second;
                auto &storage = bs.*storage_member;

                const Block &blk = storage_.block(v.this_block);

                for (int cid = 0; cid < (int)bs.desc.channels.size(); ++cid)
                {
                    const auto &ch = bs.desc.channels[cid];
                    Box3 box = LAYOUT::coupling_vertex_ghost_slab_from_cells(
                        GRID_TRAITS::cell_counts(blk),
                        ch.location,
                        v.dir1,
                        v.dir2,
                        v.dir3,
                        ch.nghost);
                    alloc_block(storage[cid][iv], kind, v.this_block, ch, box);
                }
            }
        };

        build_vertex_list(inner_vertices, TOPO::PatchKind::Inner, &CouplingBuffersForPair::inner_vertex);
        build_vertex_list(parallel_vertices, TOPO::PatchKind::Parallel, &CouplingBuffersForPair::parallel_vertex);
    }
}
