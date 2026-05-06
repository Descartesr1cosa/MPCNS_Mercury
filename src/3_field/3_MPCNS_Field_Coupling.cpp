#include "3_field/2_MPCNS_Field.h"
#include "0_basic/Error.h"

void Field::register_coupling_channel(const std::string &src,
                                      const std::string &dst,
                                      const std::string &tag,
                                      StaggerLocation location,
                                      int ncomp,
                                      int nghost)
{
    if (has_field(tag))
    {
        const auto &desc = descriptor(field_id(tag));
        const bool same = desc.location == location &&
                          desc.ncomp == ncomp &&
                          desc.nghost == nghost;
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
                (ch.ncomp == ncomp) &&
                (ch.nghost == nghost);

            if (same)
                return; // 幂等：重复注册同样的 channel，直接忽略

            std::cout << "Fatal: coupling channel tag duplicated but spec differs: "
                      << a << " -> " << b << " tag=" << tag << "\n";
            std::exit(-1);
        }

        CouplingChannelSpec spec;
        spec.tag = tag;
        spec.location = location;
        spec.ncomp = ncomp;
        spec.nghost = nghost;

        pd.channels.push_back(spec);
    };

    add_one(src, dst);
}

void Field::register_coupling_channel(const std::string &src,
                                      const std::string &dst,
                                      const std::string &field_name)
{
    const auto &desc = descriptor(field_name);
    register_coupling_channel(src, dst, field_name, desc.location, desc.ncomp, desc.nghost);
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

    // 1) 根据coupling_pairs_以及topo中patch的数量开辟[cid][ipatch]数组的空间
    for (const auto &kv : coupling_pairs_)
    {
        const PairKey &key = kv.first;
        const CouplingPairDesc &pd = kv.second;

        CouplingBuffersForPair bs;
        bs.desc = pd;

        const int nc = (int)pd.channels.size();

        bs.inner_face.assign(nc, std::vector<CouplingBufferBlock>(topo.inner_patches.size()));
        bs.parallel_face.assign(nc, std::vector<CouplingBufferBlock>(topo.parallel_patches.size()));

        bs.inner_edge.assign(nc, std::vector<CouplingBufferBlock>(topo.inner_edge_patches.size()));
        bs.parallel_edge.assign(nc, std::vector<CouplingBufferBlock>(topo.parallel_edge_patches.size()));

        bs.inner_vertex.assign(nc, std::vector<CouplingBufferBlock>(topo.inner_vertex_patches.size()));
        bs.parallel_vertex.assign(nc, std::vector<CouplingBufferBlock>(topo.parallel_vertex_patches.size()));

        coupling_buffers_[key] = std::move(bs);
    }

    // ---------- lambdas: location geometry ----------
    auto loc_delta = [](StaggerLocation loc) -> Int3
    {
        // delta=1 表示该轴是 cell-like（比 node 少1）；delta=0 表示 node-like
        switch (loc)
        {
        case StaggerLocation::Cell:
            return {1, 1, 1};
        case StaggerLocation::Node:
            return {0, 0, 0};
        case StaggerLocation::FaceXi:
            return {0, 1, 1};
        case StaggerLocation::FaceEt:
            return {1, 0, 1};
        case StaggerLocation::FaceZe:
            return {1, 1, 0};
        case StaggerLocation::EdgeXi:
            return {1, 0, 0};
        case StaggerLocation::EdgeEt:
            return {0, 1, 0};
        case StaggerLocation::EdgeZe:
            return {0, 0, 1};
        }
        return {0, 0, 0};
    };

    auto loc_inner_hi = [&](const Block &blk, StaggerLocation loc) -> Int3
    {
        const Int3 nodes = {blk.mx + 1, blk.my + 1, blk.mz + 1}; // node counts
        const Int3 d = loc_delta(loc);
        return {nodes.i - d.i, nodes.j - d.j, nodes.k - d.k}; // half-open [0,hi)
    };

    auto convert_tangent = [](int lo_n, int hi_n, int delta, int &lo, int &hi)
    {
        lo = lo_n;
        hi = (delta == 0) ? hi_n : (hi_n - 1);
    };

    // ---------- slab box makers (face/edge/vertex) ----------
    auto make_face_slab_box = [&](const Block &blk, StaggerLocation loc,
                                  const Box3 &face_node_box, int dir_code, int nghost) -> Box3
    {
        int ax = std::abs(dir_code); // 1,2,3
        int sgn = (dir_code > 0) ? +1 : -1;

        const Int3 hi_in = loc_inner_hi(blk, loc);
        const Int3 d = loc_delta(loc);

        int t1, t2;
        if (ax == 1)
        {
            t1 = 2;
            t2 = 3;
        }
        else if (ax == 2)
        {
            t1 = 1;
            t2 = 3;
        }
        else
        {
            t1 = 1;
            t2 = 2;
        }

        Box3 b{};
        // normal ghost
        if (ax == 1)
        {
            b.lo.i = (sgn < 0) ? -nghost : hi_in.i;
            b.hi.i = (sgn < 0) ? 0 : hi_in.i + nghost;
        }
        if (ax == 2)
        {
            b.lo.j = (sgn < 0) ? -nghost : hi_in.j;
            b.hi.j = (sgn < 0) ? 0 : hi_in.j + nghost;
        }
        if (ax == 3)
        {
            b.lo.k = (sgn < 0) ? -nghost : hi_in.k;
            b.hi.k = (sgn < 0) ? 0 : hi_in.k + nghost;
        }

        // tangential from node box -> loc box
        auto set_tangent = [&](int t)
        {
            int lo, hi;
            if (t == 1)
                convert_tangent(face_node_box.lo.i, face_node_box.hi.i, d.i, lo, hi);
            if (t == 2)
                convert_tangent(face_node_box.lo.j, face_node_box.hi.j, d.j, lo, hi);
            if (t == 3)
                convert_tangent(face_node_box.lo.k, face_node_box.hi.k, d.k, lo, hi);

            if (t == 1)
            {
                b.lo.i = lo;
                b.hi.i = hi;
            }
            if (t == 2)
            {
                b.lo.j = lo;
                b.hi.j = hi;
            }
            if (t == 3)
            {
                b.lo.k = lo;
                b.hi.k = hi;
            }
        };
        set_tangent(t1);
        set_tangent(t2);

        return b;
    };

    auto make_edge_slab_box = [&](const Block &blk, StaggerLocation loc,
                                  const Box3 &edge_node_box, int dir1, int dir2, int nghost) -> Box3
    {
        int a1 = std::abs(dir1), a2 = std::abs(dir2);
        int s1 = (dir1 > 0) ? +1 : -1;
        int s2 = (dir2 > 0) ? +1 : -1;

        int ae = 6 - a1 - a2; // remaining axis (1+2+3=6)

        const Int3 hi_in = loc_inner_hi(blk, loc);
        const Int3 d = loc_delta(loc);

        Box3 b{};

        auto set_ghost = [&](int ax, int sgn)
        {
            if (ax == 1)
            {
                b.lo.i = (sgn < 0) ? -nghost : hi_in.i;
                b.hi.i = (sgn < 0) ? 0 : hi_in.i + nghost;
            }
            if (ax == 2)
            {
                b.lo.j = (sgn < 0) ? -nghost : hi_in.j;
                b.hi.j = (sgn < 0) ? 0 : hi_in.j + nghost;
            }
            if (ax == 3)
            {
                b.lo.k = (sgn < 0) ? -nghost : hi_in.k;
                b.hi.k = (sgn < 0) ? 0 : hi_in.k + nghost;
            }
        };
        set_ghost(a1, s1);
        set_ghost(a2, s2);

        // along-edge axis range from node box -> loc box
        {
            int lo, hi;
            if (ae == 1)
                convert_tangent(edge_node_box.lo.i, edge_node_box.hi.i, d.i, lo, hi);
            if (ae == 2)
                convert_tangent(edge_node_box.lo.j, edge_node_box.hi.j, d.j, lo, hi);
            if (ae == 3)
                convert_tangent(edge_node_box.lo.k, edge_node_box.hi.k, d.k, lo, hi);

            if (ae == 1)
            {
                b.lo.i = lo;
                b.hi.i = hi;
            }
            if (ae == 2)
            {
                b.lo.j = lo;
                b.hi.j = hi;
            }
            if (ae == 3)
            {
                b.lo.k = lo;
                b.hi.k = hi;
            }
        }

        return b;
    };

    auto make_vertex_slab_box = [&](const Block &blk, StaggerLocation loc,
                                    int dir1, int dir2, int dir3, int nghost) -> Box3
    {
        const Int3 hi_in = loc_inner_hi(blk, loc);

        auto ghost_rng = [&](int dir, int hi_axis, int &lo, int &hi)
        {
            if (dir < 0)
            {
                lo = -nghost;
                hi = 0;
            }
            else
            {
                lo = hi_axis;
                hi = hi_axis + nghost;
            }
        };

        Box3 b{};
        ghost_rng(dir1, hi_in.i, b.lo.i, b.hi.i);
        ghost_rng(dir2, hi_in.j, b.lo.j, b.hi.j);
        ghost_rng(dir3, hi_in.k, b.lo.k, b.hi.k);
        return b;
    };

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
        cb.ncomp = ch.ncomp;
        cb.nghost = ch.nghost;

        cb.box = box;
        cb.shift = Int3{-box.lo.i, -box.lo.j, -box.lo.k};

        cb.data.SetSize(Ni, Nj, Nk, 0, ch.ncomp);
    };

    auto detect_dir_code_from_node_box = [&](const Box3 &face_node_box, const Block &blk) -> int
    {
        // 假设 topo 的 node index 是本块局部索引：node i范围 [0, mx+1)
        // 若你不是 0-based，把 0 和 (mx+1) 替换成 blk 的 node_lo / node_hi 即可
        const int nx = blk.mx + 1;
        const int ny = blk.my + 1;
        const int nz = blk.mz + 1;

        // 正常情况下 interface 的 node box 在法向方向厚度为 1
        if (face_node_box.hi.i - face_node_box.lo.i == 1)
        {
            if (face_node_box.lo.i == 0)
                return -1; // X-
            if (face_node_box.hi.i == nx)
                return +1; // X+
        }
        if (face_node_box.hi.j - face_node_box.lo.j == 1)
        {
            if (face_node_box.lo.j == 0)
                return -2; // Y-
            if (face_node_box.hi.j == ny)
                return +2; // Y+
        }
        if (face_node_box.hi.k - face_node_box.lo.k == 1)
        {
            if (face_node_box.lo.k == 0)
                return -3; // Z-
            if (face_node_box.hi.k == nz)
                return +3; // Z+
        }

        return 0; // 推不出来说明 this_box_node 不是贴边的“纯面”
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

            const Block &blk = *blocks_[p.this_block];

            // 这里要求dir_code =（±1/±2/±3）
            int dir_code = detect_dir_code_from_node_box(p.this_box_node, blk);
            if (dir_code == 0)
            {
                std::cout << "Error: cannot detect dir_code for coupling interface patch\n";
                exit(-1);
            }

            // 对于patch上的每一个channel开辟合适的空间
            for (int cid = 0; cid < (int)bs.desc.channels.size(); ++cid)
            {
                const auto &ch = bs.desc.channels[cid];
                Box3 box = make_face_slab_box(blk, ch.location, p.this_box_node, dir_code, ch.nghost);
                alloc_block(storage[cid][ip], kind, p.this_block, ch, box);
            }
        }
    };

    build_face_list(topo.inner_patches, TOPO::PatchKind::Inner, &CouplingBuffersForPair::inner_face);
    build_face_list(topo.parallel_patches, TOPO::PatchKind::Parallel, &CouplingBuffersForPair::parallel_face);

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

                const Block &blk = *blocks_[e.this_block];

                for (int cid = 0; cid < (int)bs.desc.channels.size(); ++cid)
                {
                    const auto &ch = bs.desc.channels[cid];
                    Box3 box = make_edge_slab_box(blk, ch.location, e.this_box_node, e.dir1, e.dir2, ch.nghost);
                    alloc_block(storage[cid][ie], kind, e.this_block, ch, box);
                }
            }
        };

        build_edge_list(topo.inner_edge_patches, TOPO::PatchKind::Inner, &CouplingBuffersForPair::inner_edge);
        build_edge_list(topo.parallel_edge_patches, TOPO::PatchKind::Parallel, &CouplingBuffersForPair::parallel_edge);
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

                const Block &blk = *blocks_[v.this_block];

                for (int cid = 0; cid < (int)bs.desc.channels.size(); ++cid)
                {
                    const auto &ch = bs.desc.channels[cid];
                    Box3 box = make_vertex_slab_box(blk, ch.location, v.dir1, v.dir2, v.dir3, ch.nghost);
                    alloc_block(storage[cid][iv], kind, v.this_block, ch, box);
                }
            }
        };

        build_vertex_list(topo.inner_vertex_patches, TOPO::PatchKind::Inner, &CouplingBuffersForPair::inner_vertex);
        build_vertex_list(topo.parallel_vertex_patches, TOPO::PatchKind::Parallel, &CouplingBuffersForPair::parallel_vertex);
    }
}
