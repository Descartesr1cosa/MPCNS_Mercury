#include "1_Boundary.h"
#include "0_basic/Error.h"

void MercuryBoundary::AddGroup(const BoundGroup &g)
{
    if (!halo_)
        ERROR::Abort("AddGroup: call Setup first");

    groups_[g.name] = g;

    for (auto &coupling_pair : g.coupling_pairs)
    {
        // 0) the definition of coupling src--dst should have been done + build_coupling_buffers
        if (!fld_->has_coupling_pair(coupling_pair.first, coupling_pair.second))
            return;

        const CouplingPairDesc &desc = fld_->coupling_pair(coupling_pair.first, coupling_pair.second); // description
        auto &ch = desc.channels;
        std::vector<int32_t> temp_cids;
        temp_cids.resize(0);

        for (auto field_name : g.fields)
        {
            int temp_num = -1;
            for (int cid = 0; ch.size(); cid++)
            {
                if (ch[cid].tag == field_name)
                {
                    temp_num = cid;
                    temp_cids.push_back(temp_num);
                    break;
                }
            }
            if (temp_num == -1)
            {
                std::cout << "Fatal Error, Can not find the cid in channel for field_name:  " << field_name << std::endl;
                exit(-1);
            }
        }

        groups_[g.name].fields_cids[coupling_pair] = temp_cids;
    }

    // if (g.do_halo)
    // {
    //     for (auto &fn : g.fields)
    //         halo_->register_halo_field(fn, g.halo_level);
    // }
}

void MercuryBoundary::RegisterPhysical_(const std::string &field, const std::string &region, PhysicalHandler h)
{
    bound_.RegisterPhysical(field, region, std::move(h));
}

void MercuryBoundary::RegisterCoupling_(const std::string &src, const std::string &dst,
                                        StaggerLocation loc,
                                        const std::string &channel_tag,
                                        const std::string &dst_field,
                                        CouplingHandler h)
{
    bound_.RegisterCoupling(src, dst, loc, channel_tag, dst_field, std::move(h));
}
