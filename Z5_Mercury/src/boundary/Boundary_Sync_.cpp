#include "1_Boundary.h"
#include "0_basic/Error.h"

void MercuryBoundary::Sync_(const BoundGroup &g)
{
    // CheckSetupOrAbort_("Sync_");

    std::string field_name_temp;
    std::string field_name_temp2;

    // ---------------- Stage 1: FaceOnly (1D) ----------------
    if (g.do_physical)
    {
        bound_.ApplyPhysical(g.fields);
        bound_.ApplyPhysicalCornerDefault(g.fields); // 先补角区，保证 Edge halo 的输入一致
    }

    if (g.do_halo)
    {
        for (auto &fn : g.fields)
        {
            field_name_temp = fn;
            halo_->data_trans_1DCorner(field_name_temp);
        }
    }

    if (g.do_coupling)
    {
        for (auto &pr : g.coupling_pairs)
        {
            field_name_temp = pr.first;
            field_name_temp2 = pr.second;
            std::vector<int32_t> tmp_cids = g.fields_cids.at(pr);
            halo_->coupling_trans_1DCorner(field_name_temp, field_name_temp2, tmp_cids);
            bound_.ApplyCouplingPair_1DCorner(field_name_temp, field_name_temp2, tmp_cids);
        }
    }

    // ---------------- Stage 2: Edge (2D) ----------------
    if (g.halo_level >= HaloLevel::Edge)
    {
        if (g.do_physical)
        {
            // Security：Edge run corner default again
            bound_.ApplyPhysicalCornerDefault(g.fields);
        }

        if (g.do_halo)
        {
            for (auto &fn : g.fields)
            {
                field_name_temp = fn;
                halo_->data_trans_2DCorner(field_name_temp);
            }
        }

        if (g.do_coupling)
        {
            for (auto &pr : g.coupling_pairs)
            {
                field_name_temp = pr.first;
                field_name_temp2 = pr.second;
                std::vector<int32_t> tmp_cids = g.fields_cids.at(pr);
                halo_->coupling_trans_2DCorner(field_name_temp, field_name_temp2, tmp_cids);
                bound_.ApplyCouplingPair_2DCorner(field_name_temp, field_name_temp2, tmp_cids);
            }
        }
    }

    // ---------------- Stage 3: Vertex (3D) ----------------
    if (g.halo_level >= HaloLevel::Vertex)
    {
        // if (g.do_physical)
        // {
        //     // 最后再补一次角区，保证输出/算子读到的是最终一致状态
        //     bound_.ApplyPhysicalCornerDefault(g.fields);
        // }

        if (g.do_halo)
        {
            for (auto &fn : g.fields)
            {
                field_name_temp = fn;
                halo_->data_trans_3DCorner(field_name_temp);
            }
        }

        if (g.do_coupling)
        {
            for (auto &pr : g.coupling_pairs)
            {
                field_name_temp = pr.first;
                field_name_temp2 = pr.second;
                std::vector<int32_t> tmp_cids = g.fields_cids.at(pr);
                halo_->coupling_trans_3DCorner(field_name_temp, field_name_temp2, tmp_cids);
                bound_.ApplyCouplingPair_3DCorner(field_name_temp, field_name_temp2, tmp_cids);
            }
        }
    }
}