#include "3_field/CouplingTypes.h"
#include "2_topology/TopologyBuilder.h"
#include "5_boundary/Boundary.h"
#include "0_basic/Error.h"

// ------------------------------------------------------------
// Registry: Coupling
// ------------------------------------------------------------
void BoundaryCore::RegisterCoupling(const std::string &src,
                                    const std::string &dst,
                                    StaggerLocation loc,
                                    const std::string &channel_tag,
                                    const std::string &dst_field_name,
                                    BOUND::CouplingHandler h)
{
    // 做一致性检查
    // 例如 dst_field_name 的 descriptor.location 必须等于 loc，否则后面会写错场/错位置
    const int fid = fld_->field_id(dst_field_name);
    const auto &desc = fld_->descriptor(fid);
    if (desc.location != loc)
        ERROR::Abort("[BoundaryCore] RegisterCoupling: dst_field location mismatch: " + dst_field_name);

    BOUND::CouplingKey k{src, dst, loc, channel_tag, dst_field_name};
    cpl_reg_[k] = std::move(h);
}

void BoundaryCore::RegisterCoupling(const std::string &src,
                                    const std::string &dst,
                                    StaggerLocation loc,
                                    BOUND::CouplingHandler h)
{
    BOUND::CouplingKey k{src, dst, loc, "", ""};
    cpl_reg_[k] = std::move(h);
}

namespace
{
    void apply_coupling_buffer_list(Field *fld,
                                    int fid_dst,
                                    BOUND::CouplingHandler handler,
                                    std::vector<CouplingBufferBlock> &buffers,
                                    const std::string &src,
                                    const std::string &dst,
                                    const std::string &tag)
    {
        for (auto &buf : buffers)
        {
            if (!buf.allocated)
                continue;

            FieldBlock &Udst = fld->field(fid_dst, buf.this_block);
            if (handler)
                handler(Udst, fld, buf, src, dst, tag);
            else
                BoundaryCore::DefaultCouplingCopy(Udst, fld, buf, src, dst, tag);
        }
    }
}

void BoundaryCore::ApplyCouplingPair(const std::string &src, const std::string &dst, HaloLevel stage)
{
    ApplyCouplingBuffers_(src, dst, stage, nullptr);
}

void BoundaryCore::ApplyCouplingPair(const std::string &src,
                                     const std::string &dst,
                                     HaloLevel stage,
                                     const std::vector<int32_t> &cids_fields)
{
    ApplyCouplingBuffers_(src, dst, stage, &cids_fields);
}

void BoundaryCore::ApplyCouplingPair_1DCorner(const std::string &src, const std::string &dst)
{
    ApplyCouplingBuffers_(src, dst, HaloLevel::Corner1D, nullptr);
}

void BoundaryCore::ApplyCouplingPair_1DCorner(const std::string &src, const std::string &dst, const std::vector<int32_t> &cids_fields)
{
    ApplyCouplingBuffers_(src, dst, HaloLevel::Corner1D, &cids_fields);
}

void BoundaryCore::ApplyCouplingPair_2DCorner(const std::string &src, const std::string &dst)
{
    ApplyCouplingBuffers_(src, dst, HaloLevel::Corner2D, nullptr);
}

void BoundaryCore::ApplyCouplingPair_2DCorner(const std::string &src, const std::string &dst, const std::vector<int32_t> &cids_fields)
{
    ApplyCouplingBuffers_(src, dst, HaloLevel::Corner2D, &cids_fields);
}

void BoundaryCore::ApplyCouplingPair_3DCorner(const std::string &src, const std::string &dst)
{
    ApplyCouplingBuffers_(src, dst, HaloLevel::Corner3D, nullptr);
}

void BoundaryCore::ApplyCouplingPair_3DCorner(const std::string &src, const std::string &dst, const std::vector<int32_t> &cids_fields)
{
    ApplyCouplingBuffers_(src, dst, HaloLevel::Corner3D, &cids_fields);
}

void BoundaryCore::ApplyCouplingBuffers_(const std::string &src,
                                         const std::string &dst,
                                         HaloLevel stage,
                                         const std::vector<int32_t> *cids_fields)
{
    auto &bs = fld_->coupling_buffers(src, dst);
    const auto &channels = bs.desc.channels;

    auto apply_channel = [&](int cid)
    {
        if (cid < 0 || cid >= static_cast<int>(channels.size()))
            ERROR::Abort("[BoundaryCore] ApplyCouplingPair: invalid coupling channel id.");

        const auto &ch = channels[cid];
        const std::string &tag = ch.tag;
        const StaggerLocation loc = ch.location;

        const std::string dst_field = tag;
        const int fid_dst = fld_->field_id(dst_field);
        auto handler = ResolveCoupling(src, dst, loc, tag, dst_field);

        if (stage == HaloLevel::Corner1D)
        {
            apply_coupling_buffer_list(fld_, fid_dst, handler, bs.inner_face[cid], src, dst, tag);
            apply_coupling_buffer_list(fld_, fid_dst, handler, bs.parallel_face[cid], src, dst, tag);
        }
        else if (stage == HaloLevel::Corner2D)
        {
            apply_coupling_buffer_list(fld_, fid_dst, handler, bs.inner_edge[cid], src, dst, tag);
            apply_coupling_buffer_list(fld_, fid_dst, handler, bs.parallel_edge[cid], src, dst, tag);
        }
        else if (stage == HaloLevel::Corner3D)
        {
            apply_coupling_buffer_list(fld_, fid_dst, handler, bs.inner_vertex[cid], src, dst, tag);
            apply_coupling_buffer_list(fld_, fid_dst, handler, bs.parallel_vertex[cid], src, dst, tag);
        }
        else
        {
            ERROR::Abort("[BoundaryCore] ApplyCouplingPair: unsupported halo stage.");
        }
    };

    if (cids_fields)
    {
        for (int32_t cid : *cids_fields)
            apply_channel(static_cast<int>(cid));
    }
    else
    {
        for (int cid = 0; cid < static_cast<int>(channels.size()); ++cid)
            apply_channel(cid);
    }
}

BOUND::CouplingHandler BoundaryCore::ResolveCoupling(const std::string &src,
                                                     const std::string &dst,
                                                     StaggerLocation loc,
                                                     const std::string &channel_tag,
                                                     const std::string &dst_field_name) const
{
    // 优先级（建议）：
    // 1) (src,dst,loc,channel,dst_field)
    // 2) (src,dst,loc,channel,"")
    // 3) (src,dst,loc,"","")
    //
    // 可选：pair-level default（如果你想用 RegisterCoupling(src,dst,h)）
    // 这里为了不引入 AnyLocation，我们给一个弱约定 fallback：
    // 4) (src,dst,Cell,"","") 作为 pair default
    auto find_one = [&](StaggerLocation L,
                        const std::string &ch,
                        const std::string &df) -> BOUND::CouplingHandler
    {
        BOUND::CouplingKey k;
        k.src = src;
        k.dst = dst;
        k.location = L;
        k.channel_tag = ch;
        k.dst_field_name = df;

        auto it = cpl_reg_.find(k);
        if (it != cpl_reg_.end())
            return it->second;
        return nullptr;
    };

    if (auto h = find_one(loc, channel_tag, dst_field_name))
        return h;
    if (auto h = find_one(loc, channel_tag, ""))
        return h;
    if (auto h = find_one(loc, "", ""))
        return h;
    return nullptr;
}

void BoundaryCore::DefaultCouplingCopy(FieldBlock &Udst, Field * /*fld*/,
                                       CouplingBufferBlock &buf,
                                       const std::string & /*src*/,
                                       const std::string & /*dst*/,
                                       const std::string & /*channel_tag*/)
{
    // 把 buf.box 区域的数据写入 Udst 同样的索引范围
    const Box3 &b = buf.box;
    const int ncomp = buf.ncomp;

    for (int i = b.lo.i; i < b.hi.i; ++i)
        for (int j = b.lo.j; j < b.hi.j; ++j)
            for (int k = b.lo.k; k < b.hi.k; ++k)
            {
                for (int m = 0; m < ncomp; ++m)
                {
                    Udst(i, j, k, m) = buf(i, j, k, m);
                }
            }
}
