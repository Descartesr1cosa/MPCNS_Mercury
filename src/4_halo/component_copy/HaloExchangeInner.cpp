#include "4_halo/Halo.h"
#include "0_basic/Error.h"

#include <vector>

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

    void apply_inverse_transform(const TOPO::IndexTransform &transform,
                                 int i,
                                 int j,
                                 int k,
                                 int &out_i,
                                 int &out_j,
                                 int &out_k)
    {
        const int dst[3] = {i, j, k};
        int src[3] = {0, 0, 0};
        const int offset[3] = {
            transform.offset.i,
            transform.offset.j,
            transform.offset.k};

        for (int axis = 0; axis < 3; ++axis)
            src[axis] = transform.sign[axis] * (dst[transform.perm[axis]] - offset[axis]);

        out_i = src[0];
        out_j = src[1];
        out_k = src[2];
    }

    struct PendingCopy
    {
        int recv_block = -1;
        int i = 0;
        int j = 0;
        int k = 0;
        std::size_t value_offset = 0;
    };

    struct PendingCopyScratch
    {
        std::vector<PendingCopy> ops;
        std::vector<double> values;

        void clear()
        {
            ops.clear();
            values.clear();
        }
    };

    PendingCopyScratch &pending_copy_scratch()
    {
        // Halo exchange is executed serially by each MPI rank.  Retaining
        // capacity removes both the per-call arrays and the former per-point
        // heap allocation without changing pack/apply ordering.
        static thread_local PendingCopyScratch scratch;
        scratch.clear();
        return scratch;
    }

    bool inside_field_block(const FieldBlock &fb, int i, int j, int k)
    {
        const Int3 lo = fb.get_lo();
        const Int3 hi = fb.get_hi();
        return i >= lo.i && i < hi.i &&
               j >= lo.j && j < hi.j &&
               k >= lo.k && k < hi.k;
    }

    void pack_send_box_to_transformed_recv(FieldBlock &send,
                                           int recv_block,
                                           const HaloRegion &region,
                                           int ncomp,
                                           PendingCopyScratch &pending)
    {
        const Box3 &recv_box = region.recv_box;

        for (int i = recv_box.lo.i; i < recv_box.hi.i; ++i)
            for (int j = recv_box.lo.j; j < recv_box.hi.j; ++j)
                for (int k = recv_box.lo.k; k < recv_box.hi.k; ++k)
                {
                    int send_i, send_j, send_k;
                    apply_inverse_transform(region.trans, i, j, k, send_i, send_j, send_k);
                    if (!inside_field_block(send, send_i, send_j, send_k))
                        continue;

                    PendingCopy op;
                    op.recv_block = recv_block;
                    op.i = i;
                    op.j = j;
                    op.k = k;
                    op.value_offset = pending.values.size();
                    for (int m = 0; m < ncomp; ++m)
                        pending.values.push_back(send(send_i, send_j, send_k, m));
                    pending.ops.push_back(op);
                }
    }

    void pack_transformed_send_to_recv_box(FieldBlock &send,
                                           int recv_block,
                                           const HaloRegion &region,
                                           int ncomp,
                                           PendingCopyScratch &pending)
    {
        const Box3 &recv_box = region.recv_box;

        for (int i = recv_box.lo.i; i < recv_box.hi.i; ++i)
            for (int j = recv_box.lo.j; j < recv_box.hi.j; ++j)
                for (int k = recv_box.lo.k; k < recv_box.hi.k; ++k)
                {
                    int send_i, send_j, send_k;
                    apply_transform(region.trans, i, j, k, send_i, send_j, send_k);
                    if (!inside_field_block(send, send_i, send_j, send_k))
                        continue;

                    PendingCopy op;
                    op.recv_block = recv_block;
                    op.i = i;
                    op.j = j;
                    op.k = k;
                    op.value_offset = pending.values.size();
                    for (int m = 0; m < ncomp; ++m)
                        pending.values.push_back(send(send_i, send_j, send_k, m));
                    pending.ops.push_back(op);
                }
    }

    void apply_pending_copies(Field *field,
                              int fid,
                              int ncomp,
                              const PendingCopyScratch &pending)
    {
        for (const PendingCopy &op : pending.ops)
        {
            if (op.recv_block < 0 || op.recv_block >= field->num_blocks())
                continue;
            FieldBlock &recv = field->field(fid, op.recv_block);
            if (!recv.is_allocated() ||
                !inside_field_block(recv, op.i, op.j, op.k))
                continue;
            for (int m = 0; m < ncomp; ++m)
                recv(op.i, op.j, op.k, m) = pending.values[op.value_offset + m];
        }
    }
}

void Halo::exchange_inner(const std::string &field_name)
{
    const int fid = fld_->field_id(field_name);
    const FieldDescriptor &desc = fld_->descriptor(fid);

    PatternKey key = {desc.location, desc.nghost};
    auto it = inner_patterns_.find(key);
    if (it == inner_patterns_.end())
        ERROR::Abort("[Halo] exchange_inner: missing inner face pattern for field: " + field_name);

    const HaloPattern &pat = it->second;
    const int ncomp = desc.ncomp;
    PendingCopyScratch &pending = pending_copy_scratch();

    for (const HaloRegion &r : pat.regions)
    {
        FieldBlock &fb_send = fld_->field(fid, r.this_block);
        if (!fb_send.is_allocated())
            continue;

        pack_send_box_to_transformed_recv(fb_send, r.neighbor_block, r, ncomp, pending);
    }

    apply_pending_copies(fld_, fid, ncomp, pending);
}

void Halo::exchange_inner_edge(const std::string &field_name)
{
    const int fid = fld_->field_id(field_name);
    const FieldDescriptor &desc = fld_->descriptor(fid);

    PatternKey key = {desc.location, desc.nghost};
    auto it = inner_edge_patterns_.find(key);
    if (it == inner_edge_patterns_.end())
        ERROR::Abort("[Halo] exchange_inner_edge: missing inner edge pattern for field: " + field_name);

    const HaloPattern &pat = it->second;
    const int ncomp = desc.ncomp;
    PendingCopyScratch &pending = pending_copy_scratch();

    for (const HaloRegion &r : pat.regions)
    {
        FieldBlock &fb_send = fld_->field(fid, r.neighbor_block);
        if (!fb_send.is_allocated())
            continue;

        pack_transformed_send_to_recv_box(fb_send, r.this_block, r, ncomp, pending);
    }

    apply_pending_copies(fld_, fid, ncomp, pending);
}

void Halo::exchange_inner_vertex(const std::string &field_name)
{
    const int fid = fld_->field_id(field_name);
    const FieldDescriptor &desc = fld_->descriptor(fid);

    PatternKey key = {desc.location, desc.nghost};
    auto it = inner_vertex_patterns_.find(key);
    if (it == inner_vertex_patterns_.end())
        ERROR::Abort("[Halo] exchange_inner_vertex: missing inner vertex pattern for field: " + field_name);

    const HaloPattern &pat = it->second;
    const int ncomp = desc.ncomp;
    PendingCopyScratch &pending = pending_copy_scratch();

    for (const HaloRegion &r : pat.regions)
    {
        FieldBlock &fb_send = fld_->field(fid, r.neighbor_block);
        if (!fb_send.is_allocated())
            continue;

        pack_transformed_send_to_recv_box(fb_send, r.this_block, r, ncomp, pending);
    }

    apply_pending_copies(fld_, fid, ncomp, pending);
}
