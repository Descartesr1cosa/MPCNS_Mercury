#include "4_halo/Halo.h"
#include "0_basic/Error.h"

namespace
{
    void apply_transform(const TOPO::IndexTransform &transform,
                         int i,
                         int j,
                         int k,
                         int &out_i,
                         int &out_j,
                         int &out_k)
    {
        const int src[3] = {i, j, k};
        int dst[3] = {0, 0, 0};
        const int offset[3] = {
            transform.offset.i,
            transform.offset.j,
            transform.offset.k};

        for (int axis = 0; axis < 3; ++axis)
            dst[transform.perm[axis]] = transform.sign[axis] * src[axis] + offset[axis];

        out_i = dst[0];
        out_j = dst[1];
        out_k = dst[2];
    }

    void copy_send_box_to_transformed_recv(FieldBlock &send,
                                           FieldBlock &recv,
                                           const HaloRegion &region,
                                           int ncomp)
    {
        const Box3 &send_box = region.send_box;

        for (int i = send_box.lo.i; i < send_box.hi.i; ++i)
            for (int j = send_box.lo.j; j < send_box.hi.j; ++j)
                for (int k = send_box.lo.k; k < send_box.hi.k; ++k)
                {
                    int recv_i, recv_j, recv_k;
                    apply_transform(region.trans, i, j, k, recv_i, recv_j, recv_k);

                    for (int m = 0; m < ncomp; ++m)
                        recv(recv_i, recv_j, recv_k, m) = send(i, j, k, m);
                }
    }

    void copy_transformed_send_to_recv_box(FieldBlock &send,
                                           FieldBlock &recv,
                                           const HaloRegion &region,
                                           int ncomp)
    {
        const Box3 &recv_box = region.recv_box;

        for (int i = recv_box.lo.i; i < recv_box.hi.i; ++i)
            for (int j = recv_box.lo.j; j < recv_box.hi.j; ++j)
                for (int k = recv_box.lo.k; k < recv_box.hi.k; ++k)
                {
                    int send_i, send_j, send_k;
                    apply_transform(region.trans, i, j, k, send_i, send_j, send_k);

                    for (int m = 0; m < ncomp; ++m)
                        recv(i, j, k, m) = send(send_i, send_j, send_k, m);
                }
    }
}

void Halo::exchange_inner(std::string field_name)
{
    const int fid = fld_->field_id(field_name);
    const FieldDescriptor &desc = fld_->descriptor(fid);

    PatternKey key = {desc.location, desc.nghost};
    auto it = inner_patterns_.find(key);
    if (it == inner_patterns_.end())
        ERROR::Abort("[Halo] exchange_inner: missing inner face pattern for field: " + field_name);

    const HaloPattern &pat = it->second;
    const int ncomp = desc.ncomp;

    for (const HaloRegion &r : pat.regions)
    {
        FieldBlock &fb_send = fld_->field(fid, r.this_block);
        FieldBlock &fb_recv = fld_->field(fid, r.neighbor_block);
        if (!fb_send.is_allocated() || !fb_recv.is_allocated())
            continue;

        copy_send_box_to_transformed_recv(fb_send, fb_recv, r, ncomp);
    }
}

void Halo::exchange_inner_edge(std::string field_name)
{
    const int fid = fld_->field_id(field_name);
    const FieldDescriptor &desc = fld_->descriptor(fid);

    PatternKey key = {desc.location, desc.nghost};
    auto it = inner_edge_patterns_.find(key);
    if (it == inner_edge_patterns_.end())
        ERROR::Abort("[Halo] exchange_inner_edge: missing inner edge pattern for field: " + field_name);

    const HaloPattern &pat = it->second;
    const int ncomp = desc.ncomp;

    for (const HaloRegion &r : pat.regions)
    {
        FieldBlock &fb_recv = fld_->field(fid, r.this_block);
        FieldBlock &fb_send = fld_->field(fid, r.neighbor_block);
        if (!fb_recv.is_allocated() || !fb_send.is_allocated())
            continue;

        copy_transformed_send_to_recv_box(fb_send, fb_recv, r, ncomp);
    }
}

void Halo::exchange_inner_vertex(std::string field_name)
{
    const int fid = fld_->field_id(field_name);
    const FieldDescriptor &desc = fld_->descriptor(fid);

    PatternKey key = {desc.location, desc.nghost};
    auto it = inner_vertex_patterns_.find(key);
    if (it == inner_vertex_patterns_.end())
        ERROR::Abort("[Halo] exchange_inner_vertex: missing inner vertex pattern for field: " + field_name);

    const HaloPattern &pat = it->second;
    const int ncomp = desc.ncomp;

    for (const HaloRegion &r : pat.regions)
    {
        FieldBlock &fb_recv = fld_->field(fid, r.this_block);
        FieldBlock &fb_send = fld_->field(fid, r.neighbor_block);
        if (!fb_recv.is_allocated() || !fb_send.is_allocated())
            continue;

        copy_transformed_send_to_recv_box(fb_send, fb_recv, r, ncomp);
    }
}
