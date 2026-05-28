#pragma once

#include <string>
#include <vector>
#include <map>
#include <set>
#include <stdexcept>
#include <cstdio>

#include "5_boundary/BoundaryTypes.h"

#include "1_grid/1_MPCNS_Grid.h"

namespace TOPO
{
    struct Topology;
}

class BoundaryCore
{
public:
    BoundaryCore() = default;
    ~BoundaryCore() = default;

    // Setup
    void SetUp(Grid *grd, Field *fld, TOPO::Topology *topo, Param *par, const std::vector<std::string> &field_names);

    // Physical boundary: Face stage uses registered handlers; Edge/Vertex
    // stages use default corner fill unless explicitly replaced later.
    void RegisterPhysical(const std::string &field_name,
                          const std::string &bc_name,
                          BOUND::PhysicalHandler h);
    void CheckPhysicalHandlers(const std::vector<std::string> &field_names) const;
    void ApplyPhysical(const std::string &field_name);
    void ApplyPhysical(const std::vector<std::string> &field_names);
    void ApplyPhysical(const std::string &field_name, HaloLevel stage);
    void ApplyPhysical(const std::vector<std::string> &field_names, HaloLevel stage);
    void ApplyPhysicalCornerDefault(const std::string &field_name);
    void ApplyPhysicalCornerDefault(const std::vector<std::string> &field_names);
    void ApplyPhysicalEdgeDefault(const std::string &field_name);
    void ApplyPhysicalEdgeDefault(const std::vector<std::string> &field_names);
    void ApplyPhysicalVertexDefault(const std::string &field_name);
    void ApplyPhysicalVertexDefault(const std::vector<std::string> &field_names);
    static void DefaultPhysicalCopy(FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int nghost);

    // Coupling boundary: callers pick the current Face/Edge/Vertex stage.
    void RegisterCoupling(const std::string &src,
                          const std::string &dst,
                          StaggerLocation loc,
                          const std::string &channel_tag,
                          const std::string &dst_field_name,
                          BOUND::CouplingHandler h);
    void RegisterCoupling(const std::string &src,
                          const std::string &dst,
                          StaggerLocation loc,
                          BOUND::CouplingHandler h);
    void ApplyCouplingPair(const std::string &src, const std::string &dst, HaloLevel stage);
    void ApplyCouplingPair(const std::string &src, const std::string &dst, HaloLevel stage, const std::vector<int32_t> &cids_fields);
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

private:
    void BuildPhysicalPatterns();

    void RegisterPhysical(StaggerLocation loc,
                          const std::string &field_name,
                          const std::string &bc_name,
                          BOUND::PhysicalHandler h);

    BOUND::PhysicalHandler ResolvePhysical(StaggerLocation loc,
                                           const std::string &field_name,
                                           const std::string &bc_name) const;

    void ApplyCouplingBuffers_(const std::string &src,
                               const std::string &dst,
                               HaloLevel stage,
                               const std::vector<int32_t> *cids_fields);

    BOUND::CouplingHandler ResolveCoupling(const std::string &src,
                                           const std::string &dst,
                                           StaggerLocation loc,
                                           const std::string &channel_tag,
                                           const std::string &dst_field_name) const;

private:
    Grid *grd_ = nullptr;
    Field *fld_ = nullptr;
    TOPO::Topology *topo_ = nullptr;
    Param *par_ = nullptr;

    std::set<StaggerLocation> enabled_locs_;

    // Cached physical patch list. Each region stores its one-cell inner
    // slab; runtime apply expands it to the requested ghost width.
    std::map<StaggerLocation, BOUND::PhysicalPattern> phy_patterns_;

    BOUND::PhysicalRegistry phy_reg_;
    BOUND::CouplingRegistry cpl_reg_;

    static Int3 LocDelta(StaggerLocation loc);
    static Int3 LocInnerHi(const Block &blk, StaggerLocation loc);
    static void ConvertTangent(int lo_n, int hi_n, int delta, int &lo, int &hi);

    static Box3 MakeInnerSlabBox_OneLayer(const Block &blk,
                                          StaggerLocation loc,
                                          const Box3 &face_node_box,
                                          int direction);

public:
    static Box3 MakeGhostSlabFromInner(const Box3 &inner_slab,
                                       int direction,
                                       int nghost);
};
