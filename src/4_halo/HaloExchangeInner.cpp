#include "4_halo/Halo.h"

void Halo::exchange_inner(std::string field_name)
{
    const FieldDescriptor *desc = nullptr;
    int fid = -1;
    // 遍历所有物理场
    for (int id = 0; id < fld_->num_fields(); ++id)
    {
        if (field_name == fld_->descriptor(id).name)
        {
            desc = &(fld_->descriptor(id));
            fid = id;
        }
    }
    if (desc == nullptr)
    {
        std::cout << "Fatal Error! ! ! Can not find the field:\t" << field_name << std::endl;
        exit(-1);
    }

    PatternKey key = {desc->location, desc->nghost};

    // 找到该 (location, nghost) 的 halo pattern
    auto it = inner_patterns_.find(key);
    if (it == inner_patterns_.end())
    {
        std::cout << "Can not find the  Halo pattern of field:\t" << field_name << std::endl;
        exit(-1);
    }

    const HaloPattern &pat = it->second;
    const int ncomp = desc->ncomp;

    // 遍历所有 HaloRegion：每个 region 表示一次“块 A 的 inner -> 块 B 的 ghost”
    for (const HaloRegion &r : pat.regions)
    {
        // 发送方和接收方的 FieldBlock
        FieldBlock &fb_send = fld_->field(fid, r.this_block);
        FieldBlock &fb_recv = fld_->field(fid, r.neighbor_block);
        if (!fb_send.is_allocated() || !fb_recv.is_allocated())
            continue;

        const Box3 &sb = r.send_box; // send box（在 this_block 上）
        const Box3 &rb = r.recv_box; // recv box（在 neighbor_block 上）

        // point_local -> point_nb 的索引变换：
        // nb[ perm[a] ] = sign[a] * local[ a ] + offset[a]
        const TOPO::IndexTransform &transform = r.trans;
        int offset[3] = {r.trans.offset.i,
                         r.trans.offset.j,
                         r.trans.offset.k};

        int ijk[3], tar_ijk[3];
        for (ijk[0] = sb.lo.i; ijk[0] < sb.hi.i; ijk[0]++)
            for (ijk[1] = sb.lo.j; ijk[1] < sb.hi.j; ijk[1]++)
                for (ijk[2] = sb.lo.k; ijk[2] < sb.hi.k; ijk[2]++)
                {
                    for (int ii = 0; ii < 3; ii++)
                        tar_ijk[transform.perm[ii]] = transform.sign[ii] * ijk[ii] + offset[ii];
                    // 真正的数据拷贝：inner -> ghost
                    for (int m = 0; m < ncomp; ++m)
                        fb_recv(tar_ijk[0], tar_ijk[1], tar_ijk[2], m) = fb_send(ijk[0], ijk[1], ijk[2], m);
                }
    }
}

// edge的发送与face不完全一致，由对面发送，本块接受
void Halo::exchange_inner_edge(std::string field_name)
{
    const FieldDescriptor *desc = nullptr;
    int fid = -1;
    // 遍历所有物理场
    for (int id = 0; id < fld_->num_fields(); ++id)
    {
        if (field_name == fld_->descriptor(id).name)
        {
            desc = &(fld_->descriptor(id));
            fid = id;
        }
    }
    if (desc == nullptr)
    {
        std::cout << "Fatal Error! ! ! Can not find the field:\t" << field_name << std::endl;
        exit(-1);
    }

    PatternKey key = {desc->location, desc->nghost};

    // 找到该 (location, nghost) 的 halo pattern
    auto it = inner_edge_patterns_.find(key);
    if (it == inner_edge_patterns_.end())
    {
        std::cout << "Can not find the  Halo pattern of field:\t" << field_name << std::endl;
        exit(-1);
    }

    const HaloPattern &pat = it->second;
    const int ncomp = desc->ncomp;

    // 遍历所有 HaloRegion：每个 region 表示一次“块 A 的 ghost -> 块 B 的inner ”
    for (const HaloRegion &r : pat.regions)
    {
        // 发送方和接收方的 FieldBlock
        FieldBlock &fb_recv = fld_->field(fid, r.this_block);
        FieldBlock &fb_send = fld_->field(fid, r.neighbor_block);

        if (!fb_recv.is_allocated() || !fb_send.is_allocated())
            continue;

        const Box3 &sb = r.send_box; // send box（在 neighbor_block 上）
        const Box3 &rb = r.recv_box; // recv box（在 this_block 上）

        // point_local -> point_nb 的索引变换：
        // nb[ perm[a] ] = sign[a] * local[ a ] + offset[a]
        const TOPO::IndexTransform &transform = r.trans;
        int offset[3] = {r.trans.offset.i,
                         r.trans.offset.j,
                         r.trans.offset.k};

        int ijk[3], tar_ijk[3];
        // 站在本块上循环
        for (ijk[0] = rb.lo.i; ijk[0] < rb.hi.i; ijk[0]++)
            for (ijk[1] = rb.lo.j; ijk[1] < rb.hi.j; ijk[1]++)
                for (ijk[2] = rb.lo.k; ijk[2] < rb.hi.k; ijk[2]++)
                {
                    for (int ii = 0; ii < 3; ii++)
                        tar_ijk[transform.perm[ii]] = transform.sign[ii] * ijk[ii] + offset[ii];
                    // 真正的数据拷贝：inner -> ghost
                    for (int m = 0; m < ncomp; ++m)
                        fb_recv(ijk[0], ijk[1], ijk[2], m) = fb_send(tar_ijk[0], tar_ijk[1], tar_ijk[2], m);
                }
    }
}

// edge的发送与face不完全一致，由对面发送，本块接受
void Halo::exchange_inner_vertex(std::string field_name)
{
    const FieldDescriptor *desc = nullptr;
    int fid = -1;
    // 遍历所有物理场
    for (int id = 0; id < fld_->num_fields(); ++id)
    {
        if (field_name == fld_->descriptor(id).name)
        {
            desc = &(fld_->descriptor(id));
            fid = id;
        }
    }
    if (desc == nullptr)
    {
        std::cout << "Fatal Error! ! ! Can not find the field:\t" << field_name << std::endl;
        exit(-1);
    }

    PatternKey key = {desc->location, desc->nghost};

    // 找到该 (location, nghost) 的 halo pattern
    auto it = inner_vertex_patterns_.find(key);
    if (it == inner_vertex_patterns_.end())
    {
        std::cout << "Can not find the  Halo pattern of field:\t" << field_name << std::endl;
        exit(-1);
    }

    const HaloPattern &pat = it->second;
    const int ncomp = desc->ncomp;

    // 遍历所有 HaloRegion：每个 region 表示一次“块 A 的 ghost -> 块 B 的inner ”
    for (const HaloRegion &r : pat.regions)
    {
        // 发送方和接收方的 FieldBlock
        FieldBlock &fb_recv = fld_->field(fid, r.this_block);
        FieldBlock &fb_send = fld_->field(fid, r.neighbor_block);

        if (!fb_recv.is_allocated() || !fb_send.is_allocated())
            continue;

        const Box3 &sb = r.send_box; // send box（在 neighbor_block 上）
        const Box3 &rb = r.recv_box; // recv box（在 this_block 上）

        // point_local -> point_nb 的索引变换：
        // nb[ perm[a] ] = sign[a] * local[ a ] + offset[a]
        const TOPO::IndexTransform &transform = r.trans;
        int offset[3] = {r.trans.offset.i,
                         r.trans.offset.j,
                         r.trans.offset.k};

        int ijk[3], tar_ijk[3];
        // 站在本块上循环
        for (ijk[0] = rb.lo.i; ijk[0] < rb.hi.i; ijk[0]++)
            for (ijk[1] = rb.lo.j; ijk[1] < rb.hi.j; ijk[1]++)
                for (ijk[2] = rb.lo.k; ijk[2] < rb.hi.k; ijk[2]++)
                {
                    for (int ii = 0; ii < 3; ii++)
                        tar_ijk[transform.perm[ii]] = transform.sign[ii] * ijk[ii] + offset[ii];
                    // 真正的数据拷贝：inner -> ghost
                    for (int m = 0; m < ncomp; ++m)
                        fb_recv(ijk[0], ijk[1], ijk[2], m) = fb_send(tar_ijk[0], tar_ijk[1], tar_ijk[2], m);
                }
    }
}