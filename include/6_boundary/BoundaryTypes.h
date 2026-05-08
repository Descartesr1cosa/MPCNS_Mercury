#pragma once

#include <cstdint>
#include <functional>

#include "3_field/1_Field_Block.h"
#include "3_field/2_MPCNS_Field.h"
#include "3_field/Coupling_Type.h"
namespace BOUND
{
    // =========================================================================
    // 1) Physical Boundary (来自 topo_->physical_patches)
    // =========================================================================

    // PhysicalRegion：把 topo patch 转成“循环友好”的 region 描述。
    // 注意：这里的 base_box 不包含 ghost 扩展；ghost 区域应在 Apply 时按 desc.nghost 计算。
    struct PhysicalRegion
    {
        // ---- 该 region 属于哪个 block
        int this_block = -1;
        std::string this_block_name;

        // ---- topo patch 标识信息
        int bc_id = -1;      // topo 中的 bc id
        std::string bc_name; // topo 中的 bc_name：例如 "Solid_Surface", "Farfield", ...
        int direction = 0;   // topo patch direction: +/-1,+/-2,+/-3
        Int3 cycle{};        // raw->cycle：指向 ghost 的方向(通常只有一个分量非 0)

        // 缓存：域内贴边 1 层 slab（已经是该 location 的索引体系）
        Box3 inner_slab{};

        // 运行时临时：要写入的 ghost slab（Apply 时由 inner_slab + nghost 推出）
        // 注意：inner_slab 是缓存数据；box 是临时工作数据（可不写入缓存容器中）
        Box3 box{};

        // ---- 回指 topo 原始结构（可选，用于 debug 或获取更多信息）
        const Physical_Boundary *raw = nullptr;
    };

    struct PhysicalPattern
    {
        StaggerLocation location{};
        std::vector<PhysicalRegion> regions;
    };

    // PhysicalHandler：对某个 field 的某个 block (FieldBlock& U) 的 ghost 区施加边界。
    //
    // 输入参数解释：
    // - U:     当前 field 在某个 block 上的 FieldBlock（Cell/Face/Edge 都用它）
    // - fld:   允许 handler 访问其他场（比如 B_cell, PV 等）或 descriptor
    // - r:     region 描述（含 bc_name/cycle/base_box）
    // - nghost:当前 field 的 ghost 层数（desc.nghost）
    //
    // handler 的责任：只修改 ghost 区，不要动 base_box 内侧区域（除非你明确想这么做）。
    using PhysicalHandler =
        std::function<void(FieldBlock &U, Field *fld, const PhysicalRegion &r, int nghost)>;

    // ---- Physical boundary 的 Key
    //
    // 你需要支持：
    // 1) 同一个 bc_name，在不同 location 上不同处理（Cell vs Face vs Edge）
    // 2) 同一个 bc_name，在不同 field 上不同处理（例如 U 与 B_cell）
    //
    // 所以 Key 至少包含：location + field_name + bc_name
    //
    // 约定：
    // - field_name 允许为空字符串 ""：表示“对所有 field 通用”
    // - bc_name    允许为空字符串 ""：表示“默认 handler（fallback）”
    struct PhysicalKey
    {
        StaggerLocation location{};
        std::string field_name; // "" => general for all fields
        std::string bc_name;    // "" => default for this (location,field)

        // std::map 需要严格弱序比较（不用 hash），便于你理解与调试
        bool operator<(const PhysicalKey &o) const noexcept
        {
            if (static_cast<int>(location) != static_cast<int>(o.location))
                return static_cast<int>(location) < static_cast<int>(o.location);
            if (field_name != o.field_name)
                return field_name < o.field_name;
            return bc_name < o.bc_name;
        }
    };

    // ---- Physical registry：Key -> Handler
    using PhysicalRegistry = std::map<PhysicalKey, PhysicalHandler>;

    // ---- Physical 查找优先级建议（写在这里作为“约定”，实现放在 BoundaryCore）
    // ResolvePhysicalHandler(location, field_name, bc_name) 时按顺序尝试：
    // 1) (location, field_name, bc_name)
    // 2) (location, "",        bc_name)
    // 3) (location, field_name, "")
    // 4) (location, "",        "")
    //
    // 其中 (location,"","") 作为最终 default（例如 copy）。

    // =========================================================================
    // 2) Coupling Boundary (来自 coupling buffers，而非 topo)
    // =========================================================================

    // 你给出的 coupling handler 签名：完全保留
    //
    // Udst: dst 侧 field 的 FieldBlock（写 ghost 的目标）
    // fld:  允许 handler 访问其他场/参数/几何
    // buf:  coupling buffer（Halo 已把 src 数据搬到这里）
    // src/dst: 这对耦合的名字（例如 "Solid","Fluid"）
    // channel_tag: channel 标识（例如 "U_b"、"T"、"Ehall_xi"...）
    //
    // handler 的责任：在 buf.box 范围把数据写到 Udst 的 ghost。
    using CouplingHandler = std::function<void(
        FieldBlock &Udst,
        Field *fld,
        CouplingBufferBlock &buf,
        const std::string &src,
        const std::string &dst,
        const std::string &channel_tag)>;

    // ---- Coupling 的 Key
    //
    // coupling 需要区分的维度比 physical 多：
    // - src/dst（耦合对）
    // - location（Cell/Face/Edge...）
    // - channel_tag（哪个耦合通道）
    // - dst_field_name（可选：如果 channel_tag 与 dst field 不同，需要这项）
    //
    // 为了保持最清晰，我们都放进去：
    // - dst_field_name 默认等于 channel_tag（通用情况下相同）
    //
    // 约定：
    // - channel_tag 允许 ""：表示“对该 src/dst/location 下所有 channel 通用”
    // - dst_field_name 允许 ""：表示同上；也可用来做更细覆盖
    struct CouplingKey
    {
        std::string src; // 例如 "Solid"
        std::string dst; // 例如 "Fluid"
        StaggerLocation location{};
        std::string channel_tag;    // "" => wildcard for all channels
        std::string dst_field_name; // "" => wildcard；常用时可不填（由实现决定用哪个）

        bool operator<(const CouplingKey &o) const noexcept
        {
            if (src != o.src)
                return src < o.src;
            if (dst != o.dst)
                return dst < o.dst;

            const int a = static_cast<int>(location);
            const int b = static_cast<int>(o.location);
            if (a != b)
                return a < b;

            if (channel_tag != o.channel_tag)
                return channel_tag < o.channel_tag;
            return dst_field_name < o.dst_field_name;
        }
    };

    // ---- Coupling registry：Key -> Handler
    using CouplingRegistry = std::map<CouplingKey, CouplingHandler>;

    // ---- Coupling 查找优先级建议（实现放在 CouplingBoundaryCore/BoundaryCore）
    //
    // ResolveCouplingHandler(src,dst,location,channel_tag,dst_field) 可按顺序尝试：
    //
    // 1) (src,dst,location,channel_tag,dst_field)
    // 2) (src,dst,location,channel_tag,"")
    // 3) (src,dst,location,"","")              // location default for this pair
    // 4) (src,dst,StaggerLocation::Cell? or "Any", "","")  // pair default（若你愿意支持）
    //
    // 最终没有匹配到：使用 DefaultCopyCouplingHandler（由实现层提供）。
    //
    // 说明：
    // - 第 4 条是否需要取决于你是否想要“pair 全局默认”；
    //   你也可以仅做到第 3 条，然后每个 location 设置一个默认。

    // =========================================================================
    // 3) 可选：统一的默认行为名称（仅作为约定，非必需）
    // =========================================================================

    // 对 physical：bc_name=="" 代表 default handler（fallback）
    // 对 coupling：channel_tag=="" 代表 wildcard handler（fallback）

} // namespace BOUND