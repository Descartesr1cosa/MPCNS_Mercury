#include "1_Boundary.h"
#include "0_basic/Error.h"

void LunarBoundary::Sync_(const BoundGroup &g)
{
    // CheckSetupOrAbort_("Sync_");

    std::string field_name_temp;
    std::string field_name_temp2;

    // ---------------- Stage 1: 1DCorner ----------------
    if (g.do_physical)
    {
        bound_.ApplyPhysical(g.fields);
        bound_.ApplyPhysicalCornerDefault(g.fields); // 先补角区，保证 Edge halo 的输入一致
    }

    if (g.do_halo)
    {
        if (g.name == "Eface")
            halo_->data_trans_face_axis_triplet(g.fields, HaloLevel::Corner1D);
        else
            halo_->sync_group(g.name, HaloLevel::Corner1D);
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

    // ---------------- Stage 2: 2DCorner ----------------
    if (g.halo_level >= HaloLevel::Corner2D)
    {
        if (g.do_physical)
        {
            // Security: 2DCorner runs corner default again
            bound_.ApplyPhysicalCornerDefault(g.fields);
        }

        if (g.do_halo)
        {
            if (g.name == "Eface")
                halo_->data_trans_face_axis_triplet(g.fields, HaloLevel::Corner2D);
            else
                halo_->sync_group(g.name, HaloLevel::Corner2D);
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

    // ---------------- Stage 3: 3DCorner ----------------
    if (g.halo_level >= HaloLevel::Corner3D)
    {
        // if (g.do_physical)
        // {
        //     // 最后再补一次角区，保证输出/算子读到的是最终一致状态
        //     bound_.ApplyPhysicalCornerDefault(g.fields);
        // }

        if (g.do_halo)
        {
            if (g.name == "Eface")
                halo_->data_trans_face_axis_triplet(g.fields, HaloLevel::Corner3D);
            else
                halo_->sync_group(g.name, HaloLevel::Corner3D);
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

    // The staged halo API deliberately handles only one halo depth at a time;
    // unlike Halo::sync_group(group), it does not perform the final owner ->
    // alias reconciliation.  Do that after physical/coupling updates, so a
    // topological entity shared by several blocks has one canonical value.
    // This is essential at cubic-sphere singular edges (three incident blocks).
    // Lunar E/J physics is Fluid-only, but their Solid-side storage aliases
    // are populated by the active-Fluid reconciliation before Sync().  Thus
    // the generic topology owner remains a valid transport source regardless
    // of which material owns the entity.
    if (g.do_halo &&
        (g.name == "Bface" || g.name == "Badd" ||
         g.name == "dB" || g.name == "Eedge" ||
         g.name == "Ehall" || g.name == "Eres" ||
         g.name == "Jedge" || g.name == "dJ"))
    {
        halo_->sync_owner_alias_group(g.name);
    }
}
