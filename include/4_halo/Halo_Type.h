#pragma once
#include "0_basic/TYPES.h"
#include "2_topology/2_MPCNS_Topology.h"
#include "3_field/Field_Type.h"
#include <vector>

struct HaloRegion
{
    //--------------------------------------------------------------------------
    // 对于 parallel（跨 rank）情况，还可以加：
    int this_rank;
    int neighbor_rank;
    //--------------------------------------------------------------------------
    // Block_ID
    int this_block;     // 我是哪个 block（接收方）
    int neighbor_block; // 对方 block（发送方）
    //--------------------------------------------------------------------------
    // IJK_range
    Box3 recv_box; // 我这边要填的 ghost 区 (在 this_block 的 FieldBlock 索引空间里)
    Box3 send_box; // 对方的 interior 区 (在 neighbor_block 的 FieldBlock 索引空间里)
    //--------------------------------------------------------------------------
    // Mapping patern
    TOPO::IndexTransform trans; // 从 this_block 逻辑坐标 -> neighbor_block 逻辑坐标

    int send_flag, recv_flag;
};

struct HaloPattern
{
    StaggerLocation location; // Cell / FaceXi / ...
    int nghost;               // ghost 层数

    std::vector<HaloRegion> regions;
};

enum class Direction
{
    XMinus,
    XPlus,
    YMinus,
    YPlus,
    ZMinus,
    ZPlus
};

struct EdgeMeta
{
    using PatternKey = std::pair<StaggerLocation, int>;
    PatternKey key; // {location, nghost}，方便多种 pattern 共用
    int recv_rank;  // 拿 ghost 的那边
    int send_rank;  // 提供 inner 的那边

    int recv_block; // 在 recv_rank 上的 block id (this_block on recv side)
    int send_block; // 在 send_rank 上的 block id (nb_block on recv side)

    Box3 edge_node_on_send; // 这条 edge 在 send_block 的 node 范围 (ep.nb_box_node)

    int dir1_send; // 在 send_block 上 edge 的两个方向（±1/±2/±3）
    int dir2_send;

    // 可选：如果你在 send_pattern 也需要 trans，可以携带 recv->send 的 IndexTransform
    TOPO::IndexTransform trans_recv_to_send;

    int tag; // 这一条 edge 通信的唯一标识，在 recv_rank 本地编号后写进去
};

struct VertexMeta
{
    using PatternKey = std::pair<StaggerLocation, int>;
    PatternKey key; // {location, nghost}，方便多种 pattern 共用
    int recv_rank;  // 拿 ghost 的那边
    int send_rank;  // 提供 inner 的那边

    int recv_block; // 在 recv_rank 上的 block id (this_block on recv side)
    int send_block; // 在 send_rank 上的 block id (nb_block on recv side)

    Box3 vertex_node_on_send; // 这vertex在 send_block 的 node 范围 (ep.nb_box_node)

    int dir1_send; // 在 send_block 上 edge 的两个方向（±1/±2/±3）
    int dir2_send;

    int dir3_send;

    // 可选：如果你在 send_pattern 也需要 trans，可以携带 recv->send 的 IndexTransform
    TOPO::IndexTransform trans_recv_to_send;

    int tag; // 这一条 edge 通信的唯一标识，在 recv_rank 本地编号后写进去
};
