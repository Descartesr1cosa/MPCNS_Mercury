#pragma once
#include <string>

enum class StaggerLocation
{
    Cell,   // cell center, e.g. ρ, ρu, E
    Node,   // vertex/node (可选)
    FaceXi, // 法向 xi 的面，比如 B xi
    FaceEt, // 法向 eta 的面，比如 B eta
    FaceZe, // 法向 zeta 的面，比如 B zeta
    EdgeXi, // 沿 xi 的棱，比如 E xi
    EdgeEt, // 沿 eta 的棱，比如 E eta
    EdgeZe  // 沿 zeta 的棱，比如 E zeta
};

enum class HaloLevel : int
{
    FaceOnly = 1, // 只需要 1D（面）halo
    Edge = 2,     // 需要到 2D corner（棱）
    Vertex = 3    // 需要到 3D corner（角点）
};

struct HaloFieldRequest
{
    std::string field_name;
    HaloLevel level = HaloLevel::Vertex;
};

struct FieldSyncContract
{
    // Runtime sync group. Empty means this field is allocated but does not
    // participate in the boundary/coupling/halo driver by default.
    std::string group;

    bool do_coupling = false; // register as a coupling channel
    bool do_physical = false; // require physical boundary handlers
    bool do_halo = false;     // exchange same-field halo data

    HaloLevel halo_level = HaloLevel::Vertex;
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
};
