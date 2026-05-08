#include "1_Boundary.h"
#include "0_basic/Error.h"
#include "4_halo/HaloEdgeOwner.h"

namespace
{
enum class CornerStage
{
    FaceOnly,
    Edge,
    Vertex
};

CornerStage ToCornerStage(HaloLevel level)
{
    if (level == HaloLevel::FaceOnly)
        return CornerStage::FaceOnly;
    if (level == HaloLevel::Edge)
        return CornerStage::Edge;
    return CornerStage::Vertex;
}
} // namespace

void MercuryBoundary::Sync_(const BoundGroup &g)
{
    // Z4 synchronization path:
    // 1. Physical boundary handlers still live in MercuryBoundary/BoundaryCore.
    // 2. Same-field halo, orientation-aware triplet sync, and owner-alias sync are
    //    delegated to Halo::sync_group(g.name).
    // 3. Coupling buffers are still handled by the existing staged coupling path.

    if (g.do_physical)
    {
        bound_.ApplyPhysical(g.fields);
        bound_.ApplyPhysicalCornerDefault(g.fields);
    }

    if (g.do_halo)
    {
        halo_->sync_group(g.name);
    }

    if (g.do_coupling)
    {
        ApplyCouplingStage_(g, HaloLevel::FaceOnly);

        if (g.halo_level == HaloLevel::Edge || g.halo_level == HaloLevel::Vertex)
            ApplyCouplingStage_(g, HaloLevel::Edge);

        if (g.halo_level == HaloLevel::Vertex)
            ApplyCouplingStage_(g, HaloLevel::Vertex);
    }
}

void MercuryBoundary::ApplyOwnerEdgeSync_(const BoundGroup &g)
{
    if (!g.do_owner_edge_sync || !edge_owner_pat_)
        return;

    IdTriplet fid{
        fld_->field_id(g.fields[0]),
        fld_->field_id(g.fields[1]),
        fld_->field_id(g.fields[2])};

    if (g.owner_edge_is_1form)
        HALO_OWNER::sync_edge_1form(*fld_, fid, *edge_owner_pat_);
    else
        HALO_OWNER::sync_edge_vec(*fld_, fid, *edge_owner_pat_);
}

void MercuryBoundary::ApplySameFieldHaloStage_(const BoundGroup &g, HaloLevel stage)
{
    if (!g.do_halo)
        return;

    if (g.do_owner_edge_sync && g.owner_edge_is_1form)
    {
        halo_->data_trans_edge_1form_triplet(g.fields, stage);
        return;
    }

    for (const auto &fn : g.fields)
    {
        std::string field_name = fn;
        switch (ToCornerStage(stage))
        {
        case CornerStage::FaceOnly:
            halo_->data_trans_1DCorner(field_name);
            break;
        case CornerStage::Edge:
            halo_->data_trans_2DCorner(field_name);
            break;
        case CornerStage::Vertex:
            halo_->data_trans_3DCorner(field_name);
            break;
        }
    }
}

void MercuryBoundary::ApplyCouplingStage_(const BoundGroup &g, HaloLevel stage)
{
    if (!g.do_coupling)
        return;

    for (const auto &pr : g.coupling_pairs)
    {
        std::string src = pr.first;
        std::string dst = pr.second;
        std::vector<int32_t> tmp_cids = g.fields_cids.at(pr);

        switch (ToCornerStage(stage))
        {
        case CornerStage::FaceOnly:
            halo_->coupling_trans_1DCorner(src, dst, tmp_cids);
            bound_.ApplyCouplingPair_1DCorner(src, dst, tmp_cids);
            break;
        case CornerStage::Edge:
            halo_->coupling_trans_2DCorner(src, dst, tmp_cids);
            bound_.ApplyCouplingPair_2DCorner(src, dst, tmp_cids);
            break;
        case CornerStage::Vertex:
            halo_->coupling_trans_3DCorner(src, dst, tmp_cids);
            bound_.ApplyCouplingPair_3DCorner(src, dst, tmp_cids);
            break;
        }
    }
}
