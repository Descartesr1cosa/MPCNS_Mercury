#pragma once

#include <string>
#include <vector>
#include <map>
#include <set>
#include <stdexcept>
#include <cstdio>

#include "6_boundary/Boundary_Type.h"

#include "1_grid/1_MPCNS_Grid.h"

namespace TOPO
{
    struct Topology;
}

// 维护两个 registry：
// PhysicalRegistry phy_reg_：(location, field_name, bc_name) → PhysicalHandler
// CouplingRegistry cpl_reg_：(src, dst, location, channel_tag, dst_field) → CouplingHandler
// 提供 Apply：
// ApplyPhysical(field_name)：对这个 field 的所有物理 patch 施加 BC
// ApplyCouplingPair(src,dst)：把 coupling buffer 写入 dst ghost（并允许处理）
class BoundaryCore
{
public:
    BoundaryCore() = default;
    ~BoundaryCore() = default;

    //===================================================================================
    // ------------------------------------------------------------
    // Setup & Initialization
    // ------------------------------------------------------------
    // 根据输入的物理场编号，推导所需要的Location的边界条件Pattern
    void SetUp(Grid *grd, Field *fld, TOPO::Topology *topo, Param *par, const std::vector<std::string> &field_names);

private:
    // Build cached physical patch list
    void BuildPhysicalPatterns();
    //===================================================================================

    //===================================================================================
    // Physical Boundary
public:
    // 注册：某个 field + 某个 bc_name 的处理器
    void RegisterPhysical(const std::string &field_name,
                          const std::string &bc_name,
                          BOUND::PhysicalHandler h);

    // 检测是否所有需要添加边界条件的物理场都已经设置了Handlers
    void CheckPhysicalHandlers(const std::vector<std::string> &field_names) const;

    // 施加边界条件
    void ApplyPhysical(const std::string &field_name);
    void ApplyPhysical(const std::vector<std::string> &field_names);
    // 物理角区（edge/vertex）默认拷贝处理
    void ApplyPhysicalCornerDefault(const std::string &field_name);
    void ApplyPhysicalCornerDefault(const std::vector<std::string> &field_names);

    // 提供默认的拷贝边界条件
    static void DefaultPhysicalCopy(FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int nghost);

private:
    // 注册：对某个 location + 某个 field + 某个 bc_name 的处理器
    void RegisterPhysical(StaggerLocation loc,
                          const std::string &field_name,
                          const std::string &bc_name,
                          BOUND::PhysicalHandler h);

    BOUND::PhysicalHandler ResolvePhysical(StaggerLocation loc,
                                           const std::string &field_name,
                                           const std::string &bc_name) const;

    //===================================================================================

    //===================================================================================
    // Coupling Boundary
public:
    // Registry: Coupling BC
    void RegisterCoupling(const std::string &src,
                          const std::string &dst,
                          StaggerLocation loc,
                          const std::string &channel_tag,
                          const std::string &dst_field_name,
                          BOUND::CouplingHandler h);

    // src/dst/location 下所有 channel 的默认 handler
    void RegisterCoupling(const std::string &src,
                          const std::string &dst,
                          StaggerLocation loc,
                          BOUND::CouplingHandler h);

    // 对一个 coupling pair 执行：把 buffer 写入 dst ghost
    void ApplyCouplingPair_1DCorner(const std::string &src, const std::string &dst);
    void ApplyCouplingPair_1DCorner(const std::string &src, const std::string &dst, const std::vector<int32_t> &cids_fields);
    void ApplyCouplingPair_2DCorner(const std::string &src, const std::string &dst);
    void ApplyCouplingPair_2DCorner(const std::string &src, const std::string &dst, const std::vector<int32_t> &cids_fields);
    void ApplyCouplingPair_3DCorner(const std::string &src, const std::string &dst);
    void ApplyCouplingPair_3DCorner(const std::string &src, const std::string &dst, const std::vector<int32_t> &cids_fields);
    static void DefaultCouplingCopy(FieldBlock &Udst, Field *fld,
                                    CouplingBufferBlock &buf,
                                    const std::string &src,
                                    const std::string &dst,
                                    const std::string &channel_tag);
    //===================================================================================
private:
    // ------------------------------------------------------------
    // Pointers
    // ------------------------------------------------------------
    Grid *grd_ = nullptr;
    Field *fld_ = nullptr;
    TOPO::Topology *topo_ = nullptr;
    Param *par_ = nullptr;

    // ------------------------------------------------------------
    // 存储哪些Location需要添加边界Pattern
    // ------------------------------------------------------------
    std::set<StaggerLocation> enabled_locs_;

    // ------------------------------------------------------------
    // Cached physical patch list (from topo_->phy_patterns_)
    //
    // 每个 location 一份 regions：regions 内已缓存 inner_slab（法向1层）
    // Apply 时按 inner_slab/nghost 动态算 box
    // ------------------------------------------------------------
    std::map<StaggerLocation, BOUND::PhysicalPattern> phy_patterns_;

    // ------------------------------------------------------------
    // Registries
    // ------------------------------------------------------------
    BOUND::PhysicalRegistry phy_reg_;
    BOUND::CouplingRegistry cpl_reg_;

    // ------------------------------------------------------------
    // Internal helpers: resolve handlers
    // ------------------------------------------------------------

    BOUND::CouplingHandler ResolveCoupling(const std::string &src,
                                           const std::string &dst,
                                           StaggerLocation loc,
                                           const std::string &channel_tag,
                                           const std::string &dst_field_name) const;

    //===================================================================================
    // ------------------------------------------------------------
    // Internal helpers: geometry for ghost slab on a boundary face
    // This copies the logic style used in your coupling buffer builder.
    // ------------------------------------------------------------
    static Int3 LocDelta(StaggerLocation loc);
    static Int3 LocInnerHi(const Block &blk, StaggerLocation loc);
    static void ConvertTangent(int lo_n, int hi_n, int delta, int &lo, int &hi);

    // 根据 (block, location, topo提供的node面patch盒子, direction) 推 inner_slab（法向1层）
    static Box3 MakeInnerSlabBox_OneLayer(const Block &blk,
                                          StaggerLocation loc,
                                          const Box3 &face_node_box,
                                          int direction);

public:
    // 运行时：由 inner_slab + direction + nghost 得 ghost slab（仅 O(1)）
    static Box3 MakeGhostSlabFromInner(const Box3 &inner_slab,
                                       int direction,
                                       int nghost);
    //===================================================================================
};