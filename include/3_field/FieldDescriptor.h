#pragma once
#include <string>
#include "0_basic/StaggerLocation.h"
#include "3_field/FieldValueKind.h"

enum class OwnerSyncPolicy
{
    None,
    NodeOwner,
    EdgeOwner,
    FaceOwner
};

inline const char *owner_sync_policy_name(OwnerSyncPolicy p)
{
    switch (p)
    {
    case OwnerSyncPolicy::None:
        return "None";
    case OwnerSyncPolicy::NodeOwner:
        return "NodeOwner";
    case OwnerSyncPolicy::EdgeOwner:
        return "EdgeOwner";
    case OwnerSyncPolicy::FaceOwner:
        return "FaceOwner";
    }

    return "Unknown";
}

struct FieldHaloRequest
{
    std::string field_name;

    // 用于把同一个几何对象的多个分量组成同步组。
    // 普通 field 可以使用自己的 field_name 作为 group。
    std::string sync_group;

    StaggerLocation location = StaggerLocation::Cell;
    FieldValueKind value_kind = FieldValueKind::Scalar;
    int ncomp = 1;
    int nghost = 0;
    bool do_halo = true;
    HaloLevel level = HaloLevel::Corner3D;
    OwnerSyncPolicy owner_sync = OwnerSyncPolicy::None;
    bool orientation_aware = false;
};

struct FieldSyncContract
{
    // Runtime sync group. Empty means this field is allocated but does not
    // participate in the boundary/coupling/halo driver by default.
    std::string group;

    bool do_coupling = false; // register as a coupling channel
    bool do_physical = false; // require physical boundary handlers
    bool do_halo = false;     // exchange same-field halo data

    HaloLevel halo_level = HaloLevel::Corner3D;
    OwnerSyncPolicy owner_sync = OwnerSyncPolicy::None;
    bool orientation_aware = false;
};

struct FieldDescriptor
{
    std::string name;         // "U", "B", "J", ...
    StaggerLocation location; // Cell / FaceX / ...
    int ncomp;                // 分量个数：Euler=5, MHD=8 之类
    int nghost;               // ghost 层数（新 halo 只看这个）

    // 新增：所属物理域（= block_name）。空串表示所有 block 都分配
    std::string physics = "";

    // 该 field 在运行时同步系统中的契约。物理算子临时量保持默认空契约。
    FieldSyncContract sync;

    FieldValueKind value_kind = FieldValueKind::Scalar;
};
