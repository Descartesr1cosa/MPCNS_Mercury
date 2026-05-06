#include "1_Boundary.h"
#include "0_basic/Error.h"
#include "4_halo/1_MPCNS_Halo_EdgeOwner.h"

#include <ostream>

namespace
{
const char *HaloLevelName(HaloLevel level)
{
    switch (level)
    {
    case HaloLevel::FaceOnly:
        return "FaceOnly";
    case HaloLevel::Edge:
        return "Edge";
    case HaloLevel::Vertex:
        return "Vertex";
    }
    return "Unknown";
}
} // namespace

void MercuryBoundary::AddGroup(const BoundGroup &g)
{
    if (!halo_)
        ERROR::Abort("AddGroup: call Setup first");

    BoundGroup configured = g;
    ConfigureOwnerSync_(configured);
    groups_[configured.name] = configured;

    for (auto &coupling_pair : configured.coupling_pairs)
    {
        // 0) the definition of coupling src--dst should have been done + build_coupling_buffers
        if (!fld_->has_coupling_pair(coupling_pair.first, coupling_pair.second))
            return;

        const CouplingPairDesc &desc = fld_->coupling_pair(coupling_pair.first, coupling_pair.second); // description
        auto &ch = desc.channels;
        std::vector<int32_t> temp_cids;
        temp_cids.resize(0);

        for (auto field_name : configured.fields)
        {
            int temp_num = -1;
            for (int cid = 0; cid < static_cast<int>(ch.size()); cid++)
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

        groups_[configured.name].fields_cids[coupling_pair] = temp_cids;
    }

    // if (g.do_halo)
    // {
    //     for (auto &fn : g.fields)
    //         halo_->register_halo_field(fn, g.halo_level);
    // }
}

void MercuryBoundary::AddStandardGroup_(const std::string &name,
                                        const std::vector<std::string> &fields,
                                        bool do_coupling,
                                        bool do_physical,
                                        bool do_halo,
                                        HaloLevel halo_level)
{
    BoundGroup group;
    group.name = name;
    group.fields = fields;
    group.do_coupling = do_coupling;
    group.do_physical = do_physical;
    group.do_halo = do_halo;
    group.halo_level = halo_level;
    if (do_coupling)
        group.coupling_pairs = {{"Solid", "Fluid"}, {"Fluid", "Solid"}};
    AddGroup(group);
}

void MercuryBoundary::DescribeGroups(std::ostream &os) const
{
    os << "[MercuryBoundary] Sync groups\n";
    for (const auto &entry : groups_)
    {
        const auto &g = entry.second;
        os << "  - " << g.name << ": fields={";
        for (std::size_t n = 0; n < g.fields.size(); ++n)
        {
            if (n)
                os << ", ";
            os << g.fields[n];
        }
        os << "} stages={";
        bool need_sep = false;
        if (g.do_physical)
        {
            os << "physical";
            need_sep = true;
        }
        if (g.do_owner_edge_sync)
        {
            os << (need_sep ? ", " : "") << "owner-edge:"
               << (g.owner_edge_is_1form ? "1form" : "vec");
            need_sep = true;
        }
        if (g.do_halo)
        {
            os << (need_sep ? ", " : "") << "halo:" << HaloLevelName(g.halo_level);
            need_sep = true;
        }
        if (g.do_coupling)
            os << (need_sep ? ", " : "") << "coupling";
        os << "}";

        if (!g.coupling_pairs.empty())
        {
            os << " pairs={";
            for (std::size_t n = 0; n < g.coupling_pairs.size(); ++n)
            {
                if (n)
                    os << ", ";
                os << g.coupling_pairs[n].first << "->" << g.coupling_pairs[n].second;
            }
            os << "}";
        }
        os << "\n";
    }
}

bool MercuryBoundary::IsEdgeTripletGroup_(const std::vector<std::string> &fields) const
{
    if (!fld_ || fields.size() != 3)
        return false;

    for (const auto &fn : fields)
    {
        if (!fld_->has_field(fn))
            return false;
    }

    const auto loc0 = fld_->descriptor(fld_->field_id(fields[0])).location;
    const auto loc1 = fld_->descriptor(fld_->field_id(fields[1])).location;
    const auto loc2 = fld_->descriptor(fld_->field_id(fields[2])).location;

    return loc0 == StaggerLocation::EdgeXi &&
           loc1 == StaggerLocation::EdgeEt &&
           loc2 == StaggerLocation::EdgeZe;
}

void MercuryBoundary::ConfigureOwnerSync_(BoundGroup &g) const
{
    if (!edge_owner_pat_)
        return;

    if (IsEdgeTripletGroup_(g.fields))
    {
        g.do_owner_edge_sync = true;
        g.owner_edge_is_1form = true;
    }
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
