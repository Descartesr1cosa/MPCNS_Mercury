#pragma once

#include <vector>
#include <string>

#include "2_topology/2_MPCNS_Topology.h"
#include "3_field/Field_Type.h"
#include "3_field/Field_Array.h"

// 一个“物理块对”的唯一标识：src -> dst
struct PhysPair
{
    std::string src; // e.g. "SOLID"
    std::string dst; // e.g. "FLUID"

    bool operator==(const PhysPair &o) const
    {
        return src == o.src && dst == o.dst;
    }
};

// 该物理域对下的一条“耦合通道”(payload)
// 说明：同一个 (src->dst) 下，tag 必须唯一
struct CouplingChannelSpec
{
    std::string tag;          // e.g. "T", "B", "J", "Ehall", "u_wall"
    StaggerLocation location; // Cell / FaceXi / FaceEt / FaceZe / EdgeXi / ...

    int ncomp = 0;  // 分量数：T=1, u=3, tensor=6, ...
    int nghost = 1; // 法向 slab 厚度（虚网格层数）
};

// (src->dst) 的耦合定义：唯一一条关系，但允许多个 channel
// 注意src->dst的有向性
struct CouplingPairDesc
{
    PhysPair pair;                             // 唯一：src->dst
    std::vector<CouplingChannelSpec> channels; // 多通道 例如Solid给Fluid可以传输T (Cell) 也可传输 B (Face) 甚至还有E (Edge)
};

struct CouplingBufferBlock
{
    bool allocated = false;

    int this_block = -1;  // 接收侧 block id
    TOPO::PatchKind kind; // Inner / Parallel（Physical 通常不需要 coupling buffer）

    std::string src_name; // 方便调试
    std::string dst_name; // 方便调试
    std::string tag;      // channel tag

    StaggerLocation location;
    int ncomp = 0;
    int nghost = 1;

    Box3 box;   // 该 location 逻辑索引空间的 slab 范围（允许 lo 为负）
    Int3 shift; // shift = -box.lo，用于把逻辑索引映射到 data 的 0-based

    Vector data; // (Ni, Nj, Nk, ncomp)

    inline double &operator()(int i, int j, int k, int m)
    {
        // 不做 allocated 检查，按你的习惯；上层需保证只访问 allocated==true 的 buffer
        return data(i + shift.i, j + shift.j, k + shift.k, m);
    }
};

// 每个 pair：channels 数量 = desc.channels.size()
// face/edge/vertex：分别对齐 topo 的 inner/parallel patch vector 下标
struct CouplingBuffersForPair
{
    CouplingPairDesc desc;

    // [cid][ip]  cid=channels index; ip=patch index in topo.xxx_patches
    std::vector<std::vector<CouplingBufferBlock>> inner_face, parallel_face;
    std::vector<std::vector<CouplingBufferBlock>> inner_edge, parallel_edge;
    std::vector<std::vector<CouplingBufferBlock>> inner_vertex, parallel_vertex;
};