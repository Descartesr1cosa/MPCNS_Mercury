#include "4_halo/Halo.h"
#include "2_topology/TopologyView.h"
#include "4_halo/detail/HaloBuildTools.h"

void Halo::coupling_inner_face(const std::string &src, const std::string &dst)
{
    // 0) 必须先有 coupling 定义 + 已经 build_coupling_buffers
    if (!fld_->has_coupling_pair(src, dst))
        return;

    CouplingBuffersForPair &cb = fld_->coupling_buffers(src, dst); // 缓冲区
    const CouplingPairDesc &desc = fld_->coupling_pair(src, dst);  // 描述

    // inner_face 为空则无需做
    if (cb.inner_face.empty())
        return;

    // 1) 遍历 inner face patches：用 ip 作为 buffer 的第二维下标
    const auto inner_faces = TOPO_VIEW::inner_faces(*topo_);
    for (size_t ip = 0; ip < inner_faces.size(); ++ip)
    {
        const auto &p = inner_faces[ip];

        if (!p.is_coupling)
            continue;

        // 仅处理 dst 接收侧方向：src -> dst
        // InterfacePatch 语义：this = 接收侧（dst），nb = 发送侧（src）
        if (p.this_block_name != dst || p.nb_block_name != src)
            continue;

        const int src_block = p.nb_block;        // 提供数据的块（src）
        const TOPO::IndexTransform &T = p.trans; // dst(this) -> src(nb)

        // 2) 遍历所有 channel：用 cid 作为 buffer 第一维下标
        const size_t nch = desc.channels.size();
        for (size_t cid = 0; cid < nch; ++cid)
        {
            // 防御式：保证 inner_face[cid][ip] 存在
            if (cid >= cb.inner_face.size())
            {
                std::cout << "Error! Not Enough Buffers for channel id: " << cid << std::endl;
                exit(-1);
            }
            if (ip >= cb.inner_face[cid].size())
            {
                std::cout << "Error! Not Enough Buffers for patch id: " << ip << std::endl;
                exit(-1);
            }

            CouplingBufferBlock &buf = cb.inner_face[cid][ip];
            if (!buf.allocated)
                continue; // 允许“只分配某些patch/某些channel”

            const CouplingChannelSpec &ch = desc.channels[cid];
            const std::string dst_tag = (!buf.tag.empty()) ? buf.tag : ch.tag;
            std::string src_tag = dst_tag;
            double factor = 1.0;
            if (coupling_channel_needs_form_transfer_(ch))
            {
                const int dst_axis = coupling_form_axis_from_location_(ch.location);
                const int src_axis = coupling_src_axis_from_dst_to_src_transform_(T, dst_axis);
                src_tag = find_triplet_field_name_(dst_tag, ch.value_kind, src_axis);
                factor = static_cast<double>(
                    coupling_form_orientation_sign_dst_to_src_(ch, T, dst_axis));
            }
            const int fid = fld_->field_id(src_tag);

            FieldBlock &fb_src = fld_->field(fid, src_block);

            const Box3 &b = buf.box;
            const int ncomp = buf.ncomp;

            // 3) buf.box 是 dst 侧接收区域（通常是 ghost slab）
            //    对 box 内每个点 (i,j,k)，通过 trans 映射到 src 侧 (is,js,ks) 读值
            for (int i = b.lo.i; i < b.hi.i; ++i)
            {
                for (int j = b.lo.j; j < b.hi.j; ++j)
                {
                    for (int k = b.lo.k; k < b.hi.k; ++k)
                    {
                        int is, js, ks;
                        HALO_TOOLS::apply_transform(T, i, j, k, is, js, ks);

                        for (int m = 0; m < ncomp; ++m)
                        {
                            buf(i, j, k, m) = factor * fb_src(is, js, ks, m);
                        }
                    }
                }
            }
        }
    }
}

void Halo::coupling_inner_face(const std::string &src, const std::string &dst, std::vector<int32_t> &field_cids)
{
    // 0) 必须先有 coupling 定义 + 已经 build_coupling_buffers
    if (!fld_->has_coupling_pair(src, dst))
        return;

    CouplingBuffersForPair &cb = fld_->coupling_buffers(src, dst); // 缓冲区
    const CouplingPairDesc &desc = fld_->coupling_pair(src, dst);  // 描述

    // inner_face 为空则无需做
    if (cb.inner_face.empty())
        return;

    // 1) 遍历 inner face patches：用 ip 作为 buffer 的第二维下标
    const auto inner_faces = TOPO_VIEW::inner_faces(*topo_);
    for (size_t ip = 0; ip < inner_faces.size(); ++ip)
    {
        const auto &p = inner_faces[ip];

        if (!p.is_coupling)
            continue;

        // 仅处理 dst 接收侧方向：src -> dst
        // InterfacePatch 语义：this = 接收侧（dst），nb = 发送侧（src）
        if (p.this_block_name != dst || p.nb_block_name != src)
            continue;

        const int src_block = p.nb_block;        // 提供数据的块（src）
        const TOPO::IndexTransform &T = p.trans; // dst(this) -> src(nb)

        // 2) 遍历所有 channel：用 cid 作为 buffer 第一维下标
        const size_t nch = desc.channels.size();

        for (auto cid_for_string : field_cids)
        {
            // 防御式：保证 inner_face[cid][ip] 存在
            if (cid_for_string >= cb.inner_face.size())
            {
                std::cout << "Error! Not Enough Buffers for channel id: " << cid_for_string << std::endl;
                exit(-1);
            }
            if (ip >= cb.inner_face[cid_for_string].size())
            {
                std::cout << "Error! Not Enough Buffers for patch id: " << ip << std::endl;
                exit(-1);
            }

            CouplingBufferBlock &buf = cb.inner_face[cid_for_string][ip];
            if (!buf.allocated)
                continue; // 允许“只分配某些patch/某些channel”

            const CouplingChannelSpec &ch = desc.channels[cid_for_string];
            const std::string dst_tag = (!buf.tag.empty()) ? buf.tag : ch.tag;
            std::string src_tag = dst_tag;
            double factor = 1.0;
            if (coupling_channel_needs_form_transfer_(ch))
            {
                const int dst_axis = coupling_form_axis_from_location_(ch.location);
                const int src_axis = coupling_src_axis_from_dst_to_src_transform_(T, dst_axis);
                src_tag = find_triplet_field_name_(dst_tag, ch.value_kind, src_axis);
                factor = static_cast<double>(
                    coupling_form_orientation_sign_dst_to_src_(ch, T, dst_axis));
            }
            const int fid = fld_->field_id(src_tag);

            FieldBlock &fb_src = fld_->field(fid, src_block);

            const Box3 &b = buf.box;
            const int ncomp = buf.ncomp;

            // 3) buf.box 是 dst 侧接收区域（通常是 ghost slab）
            //    对 box 内每个点 (i,j,k)，通过 trans 映射到 src 侧 (is,js,ks) 读值
            for (int i = b.lo.i; i < b.hi.i; ++i)
            {
                for (int j = b.lo.j; j < b.hi.j; ++j)
                {
                    for (int k = b.lo.k; k < b.hi.k; ++k)
                    {
                        int is, js, ks;
                        HALO_TOOLS::apply_transform(T, i, j, k, is, js, ks);

                        for (int m = 0; m < ncomp; ++m)
                        {
                            buf(i, j, k, m) = factor * fb_src(is, js, ks, m);
                        }
                    }
                }
            }
        }
    }
}

void Halo::coupling_inner_edge(const std::string &src, const std::string &dst)
{
    if (!fld_->has_coupling_pair(src, dst))
        return;

    CouplingBuffersForPair &cb = fld_->coupling_buffers(src, dst);
    const CouplingPairDesc &desc = fld_->coupling_pair(src, dst);

    if (cb.inner_edge.empty())
        return;

    const auto &inner_edges = TOPO_VIEW::edge_patches(*topo_, TOPO::PatchKind::Inner);
    for (size_t ie = 0; ie < inner_edges.size(); ++ie)
    {
        const TOPO::EdgePatch &e = inner_edges[ie];

        if (!e.is_coupling)
            continue;

        // 只处理 src -> dst（nb=src, this=dst）
        if (e.this_block_name != dst || e.nb_block_name != src)
            continue;

        const int src_block = e.nb_block;        // src block (send)
        const TOPO::IndexTransform &T = e.trans; // dst(this) -> src(nb)

        const size_t nch = desc.channels.size();
        for (size_t cid = 0; cid < nch; ++cid)
        {
            if (cid >= cb.inner_edge.size())
            {
                std::cout << "Error! Not Enough inner_edge buffers for channel id: " << cid << std::endl;
                exit(-1);
            }
            if (ie >= cb.inner_edge[cid].size())
            {
                std::cout << "Error! Not Enough inner_edge buffers for patch id: " << ie << std::endl;
                exit(-1);
            }

            CouplingBufferBlock &buf = cb.inner_edge[cid][ie];
            if (!buf.allocated)
                continue;

            const CouplingChannelSpec &ch = desc.channels[cid];
            const std::string dst_tag = (!buf.tag.empty()) ? buf.tag : ch.tag;
            std::string src_tag = dst_tag;
            double factor = 1.0;
            if (coupling_channel_needs_form_transfer_(ch))
            {
                const int dst_axis = coupling_form_axis_from_location_(ch.location);
                const int src_axis = coupling_src_axis_from_dst_to_src_transform_(T, dst_axis);
                src_tag = find_triplet_field_name_(dst_tag, ch.value_kind, src_axis);
                factor = static_cast<double>(
                    coupling_form_orientation_sign_dst_to_src_(ch, T, dst_axis));
            }
            const int fid = fld_->field_id(src_tag);

            FieldBlock &fb_src = fld_->field(fid, src_block);

            const Box3 &b = buf.box;
            const int ncomp = buf.ncomp;

            for (int i = b.lo.i; i < b.hi.i; ++i)
                for (int j = b.lo.j; j < b.hi.j; ++j)
                    for (int k = b.lo.k; k < b.hi.k; ++k)
                    {
                        int is, js, ks;
                        HALO_TOOLS::apply_transform(T, i, j, k, is, js, ks);

                        for (int m = 0; m < ncomp; ++m)
                            buf(i, j, k, m) = factor * fb_src(is, js, ks, m);
                    }
        }
    }
}

void Halo::coupling_inner_edge(const std::string &src, const std::string &dst, std::vector<int32_t> &field_cids)
{
    if (!fld_->has_coupling_pair(src, dst))
        return;

    CouplingBuffersForPair &cb = fld_->coupling_buffers(src, dst);
    const CouplingPairDesc &desc = fld_->coupling_pair(src, dst);

    if (cb.inner_edge.empty())
        return;

    const auto &inner_edges = TOPO_VIEW::edge_patches(*topo_, TOPO::PatchKind::Inner);
    for (size_t ie = 0; ie < inner_edges.size(); ++ie)
    {
        const TOPO::EdgePatch &e = inner_edges[ie];

        if (!e.is_coupling)
            continue;

        // 只处理 src -> dst（nb=src, this=dst）
        if (e.this_block_name != dst || e.nb_block_name != src)
            continue;

        const int src_block = e.nb_block;        // src block (send)
        const TOPO::IndexTransform &T = e.trans; // dst(this) -> src(nb)

        const size_t nch = desc.channels.size();
        for (auto cid : field_cids)
        {
            if (cid >= cb.inner_edge.size())
            {
                std::cout << "Error! Not Enough inner_edge buffers for channel id: " << cid << std::endl;
                exit(-1);
            }
            if (ie >= cb.inner_edge[cid].size())
            {
                std::cout << "Error! Not Enough inner_edge buffers for patch id: " << ie << std::endl;
                exit(-1);
            }

            CouplingBufferBlock &buf = cb.inner_edge[cid][ie];
            if (!buf.allocated)
                continue;

            const CouplingChannelSpec &ch = desc.channels[cid];
            const std::string dst_tag = (!buf.tag.empty()) ? buf.tag : ch.tag;
            std::string src_tag = dst_tag;
            double factor = 1.0;
            if (coupling_channel_needs_form_transfer_(ch))
            {
                const int dst_axis = coupling_form_axis_from_location_(ch.location);
                const int src_axis = coupling_src_axis_from_dst_to_src_transform_(T, dst_axis);
                src_tag = find_triplet_field_name_(dst_tag, ch.value_kind, src_axis);
                factor = static_cast<double>(
                    coupling_form_orientation_sign_dst_to_src_(ch, T, dst_axis));
            }
            const int fid = fld_->field_id(src_tag);

            FieldBlock &fb_src = fld_->field(fid, src_block);

            const Box3 &b = buf.box;
            const int ncomp = buf.ncomp;

            for (int i = b.lo.i; i < b.hi.i; ++i)
                for (int j = b.lo.j; j < b.hi.j; ++j)
                    for (int k = b.lo.k; k < b.hi.k; ++k)
                    {
                        int is, js, ks;
                        HALO_TOOLS::apply_transform(T, i, j, k, is, js, ks);

                        for (int m = 0; m < ncomp; ++m)
                            buf(i, j, k, m) = factor * fb_src(is, js, ks, m);
                    }
        }
    }
}

void Halo::coupling_inner_vertex(const std::string &src, const std::string &dst)
{
    if (!fld_->has_coupling_pair(src, dst))
        return;

    CouplingBuffersForPair &cb = fld_->coupling_buffers(src, dst);
    const CouplingPairDesc &desc = fld_->coupling_pair(src, dst);

    if (cb.inner_vertex.empty())
        return;

    const auto &inner_vertices = TOPO_VIEW::vertex_patches(*topo_, TOPO::PatchKind::Inner);
    for (size_t iv = 0; iv < inner_vertices.size(); ++iv)
    {
        const TOPO::VertexPatch &v = inner_vertices[iv];

        if (!v.is_coupling)
            continue;

        // 只处理 src -> dst（nb=src, this=dst）
        if (v.this_block_name != dst || v.nb_block_name != src)
            continue;

        const int src_block = v.nb_block;        // src block (send)
        const TOPO::IndexTransform &T = v.trans; // dst(this) -> src(nb)

        const size_t nch = desc.channels.size();
        for (size_t cid = 0; cid < nch; ++cid)
        {
            if (cid >= cb.inner_vertex.size())
            {
                std::cout << "Error! Not Enough inner_vertex buffers for channel id: " << cid << std::endl;
                exit(-1);
            }
            if (iv >= cb.inner_vertex[cid].size())
            {
                std::cout << "Error! Not Enough inner_vertex buffers for patch id: " << iv << std::endl;
                exit(-1);
            }

            CouplingBufferBlock &buf = cb.inner_vertex[cid][iv];
            if (!buf.allocated)
                continue;

            const CouplingChannelSpec &ch = desc.channels[cid];
            const std::string dst_tag = (!buf.tag.empty()) ? buf.tag : ch.tag;
            std::string src_tag = dst_tag;
            double factor = 1.0;
            if (coupling_channel_needs_form_transfer_(ch))
            {
                const int dst_axis = coupling_form_axis_from_location_(ch.location);
                const int src_axis = coupling_src_axis_from_dst_to_src_transform_(T, dst_axis);
                src_tag = find_triplet_field_name_(dst_tag, ch.value_kind, src_axis);
                factor = static_cast<double>(
                    coupling_form_orientation_sign_dst_to_src_(ch, T, dst_axis));
            }
            const int fid = fld_->field_id(src_tag);

            FieldBlock &fb_src = fld_->field(fid, src_block);

            const Box3 &b = buf.box;
            const int ncomp = buf.ncomp;

            for (int i = b.lo.i; i < b.hi.i; ++i)
                for (int j = b.lo.j; j < b.hi.j; ++j)
                    for (int k = b.lo.k; k < b.hi.k; ++k)
                    {
                        int is, js, ks;
                        HALO_TOOLS::apply_transform(T, i, j, k, is, js, ks);

                        for (int m = 0; m < ncomp; ++m)
                            buf(i, j, k, m) = factor * fb_src(is, js, ks, m);
                    }
        }
    }
}

void Halo::coupling_inner_vertex(const std::string &src, const std::string &dst, std::vector<int32_t> &field_cids)
{
    if (!fld_->has_coupling_pair(src, dst))
        return;

    CouplingBuffersForPair &cb = fld_->coupling_buffers(src, dst);
    const CouplingPairDesc &desc = fld_->coupling_pair(src, dst);

    if (cb.inner_vertex.empty())
        return;

    const auto &inner_vertices = TOPO_VIEW::vertex_patches(*topo_, TOPO::PatchKind::Inner);
    for (size_t iv = 0; iv < inner_vertices.size(); ++iv)
    {
        const TOPO::VertexPatch &v = inner_vertices[iv];

        if (!v.is_coupling)
            continue;

        // 只处理 src -> dst（nb=src, this=dst）
        if (v.this_block_name != dst || v.nb_block_name != src)
            continue;

        const int src_block = v.nb_block;        // src block (send)
        const TOPO::IndexTransform &T = v.trans; // dst(this) -> src(nb)

        const size_t nch = desc.channels.size();
        for (auto cid : field_cids)
        {
            if (cid >= cb.inner_vertex.size())
            {
                std::cout << "Error! Not Enough inner_vertex buffers for channel id: " << cid << std::endl;
                exit(-1);
            }
            if (iv >= cb.inner_vertex[cid].size())
            {
                std::cout << "Error! Not Enough inner_vertex buffers for patch id: " << iv << std::endl;
                exit(-1);
            }

            CouplingBufferBlock &buf = cb.inner_vertex[cid][iv];
            if (!buf.allocated)
                continue;

            const CouplingChannelSpec &ch = desc.channels[cid];
            const std::string dst_tag = (!buf.tag.empty()) ? buf.tag : ch.tag;
            std::string src_tag = dst_tag;
            double factor = 1.0;
            if (coupling_channel_needs_form_transfer_(ch))
            {
                const int dst_axis = coupling_form_axis_from_location_(ch.location);
                const int src_axis = coupling_src_axis_from_dst_to_src_transform_(T, dst_axis);
                src_tag = find_triplet_field_name_(dst_tag, ch.value_kind, src_axis);
                factor = static_cast<double>(
                    coupling_form_orientation_sign_dst_to_src_(ch, T, dst_axis));
            }
            const int fid = fld_->field_id(src_tag);

            FieldBlock &fb_src = fld_->field(fid, src_block);

            const Box3 &b = buf.box;
            const int ncomp = buf.ncomp;

            for (int i = b.lo.i; i < b.hi.i; ++i)
                for (int j = b.lo.j; j < b.hi.j; ++j)
                    for (int k = b.lo.k; k < b.hi.k; ++k)
                    {
                        int is, js, ks;
                        HALO_TOOLS::apply_transform(T, i, j, k, is, js, ks);

                        for (int m = 0; m < ncomp; ++m)
                            buf(i, j, k, m) = factor * fb_src(is, js, ks, m);
                    }
        }
    }
}
