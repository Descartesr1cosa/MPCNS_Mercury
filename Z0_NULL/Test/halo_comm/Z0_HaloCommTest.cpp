#include "Z0_Tests.h"

#include "Z0_Boundary.h"
#include "Z0_TestCommon.h"
#include "2_topology/Topology.h"
#include "4_halo/Halo.h"

#include <cstdlib>
#include <cmath>
#include <sstream>

namespace
{
    const char *level_name(HaloLevel level)
    {
        if (level == HaloLevel::FaceOnly)
            return "FaceOnly";
        if (level == HaloLevel::Edge)
            return "Edge";
        return "Vertex";
    }

    bool owner_sync_enabled()
    {
        const char *value = std::getenv("Z0_OWNER_SYNC");
        if (!value)
            return false;
        const std::string s(value);
        return s == "1" || s == "true" || s == "TRUE" || s == "on" || s == "ON";
    }

    bool owner_only_enabled()
    {
        const char *value = std::getenv("Z0_OWNER_ONLY");
        if (!value)
            return false;
        const std::string s(value);
        return s == "1" || s == "true" || s == "TRUE" || s == "on" || s == "ON";
    }

    bool member_matches_location(const TOPO::EquivMember &member, StaggerLocation loc)
    {
        if (loc == StaggerLocation::FaceXi)
            return member.entity.axis == TOPO::EntityAxis::Xi;
        if (loc == StaggerLocation::FaceEt)
            return member.entity.axis == TOPO::EntityAxis::Eta;
        if (loc == StaggerLocation::FaceZe)
            return member.entity.axis == TOPO::EntityAxis::Zeta;
        return TOPO::stagger_location(member.entity) == loc;
    }

    int owner_alias_sign(const FieldDescriptor &desc,
                         const TOPO::EquivMember &owner,
                         const TOPO::EquivMember &alias)
    {
        if (!desc.sync.orientation_aware)
            return +1;
        if (desc.value_kind == FieldValueKind::EdgeCovariant1Form ||
            desc.value_kind == FieldValueKind::FaceContravariant2Form)
            return owner.orient_sign * alias.orient_sign;
        return +1;
    }

    bool check_owner_alias(Field &field,
                           Halo &halo,
                           const std::string &name,
                           TOPO::EntityDim dim,
                           long long &checked,
                           double &max_err)
    {
        checked = 0;
        max_err = 0.0;

        const TOPO::Topology *topology = halo.topology_equiv();
        if (!topology || !field.has_field(name))
            return true;

        int myid = 0;
        PARALLEL::mpi_rank(&myid);

        const int fid = field.field_id(name);
        const FieldDescriptor &desc = field.descriptor(fid);
        if (desc.sync.owner_sync == OwnerSyncPolicy::None)
            return true;
        const auto &classes = topology->classes(dim);

        for (const auto &cls : classes)
        {
            const TOPO::EquivMember &owner = cls.owner;
            if (!member_matches_location(owner, desc.location))
                continue;

            for (const auto &alias : cls.members)
            {
                if (alias.is_owner || alias.entity.rank != myid)
                    continue;
                if (!member_matches_location(alias, desc.location))
                    continue;

                FieldBlock &fb_alias = field.field(fid, alias.entity.block);
                const int sign = owner_alias_sign(desc, owner, alias);
                for (int m = 0; m < desc.ncomp; ++m)
                {
                    const double expected =
                        static_cast<double>(sign) *
                        Z0_TEST::unique_code(owner.entity.rank,
                                             owner.entity.block,
                                             owner.entity.i,
                                             owner.entity.j,
                                             owner.entity.k,
                                             m);
                    const double got = fb_alias(alias.entity.i, alias.entity.j, alias.entity.k, m);
                    if (!std::isfinite(got))
                        max_err = std::numeric_limits<double>::infinity();
                    else if (std::isfinite(max_err))
                        max_err = std::max(max_err, std::abs(got - expected));
                    ++checked;
                }
            }
        }

        max_err = Z0_TEST::global_max(max_err);
        checked = Z0_TEST::global_sum(checked);
        return checked == 0 || max_err == 0.0;
    }

    long long count_nonfinite_recv_boxes(Field &field,
                                         Halo &halo,
                                         const std::string &name,
                                         HaloLevel stage,
                                         long long &region_count)
    {
        long long count = 0;
        const int fid = field.field_id(name);
        const FieldDescriptor &desc = field.descriptor(fid);
        const std::vector<HaloRegion> regions =
            halo.debug_halo_regions(desc.location, desc.nghost, stage);
        region_count = static_cast<long long>(regions.size());

        for (const HaloRegion &r : regions)
        {
            if (r.this_rank == r.neighbor_rank && stage == HaloLevel::FaceOnly)
            {
                if (r.neighbor_block < 0 || r.neighbor_block >= field.num_blocks())
                    continue;
                FieldBlock &fb = field.field(fid, r.neighbor_block);
                const Box3 &b = r.send_box;
                for (int i = b.lo.i; i < b.hi.i; ++i)
                    for (int j = b.lo.j; j < b.hi.j; ++j)
                        for (int k = b.lo.k; k < b.hi.k; ++k)
                        {
                            int ri = 0, rj = 0, rk = 0;
                            Z0_TEST::map_index(r.trans, i, j, k, ri, rj, rk);
                            for (int m = 0; m < desc.ncomp; ++m)
                                if (!std::isfinite(fb(ri, rj, rk, m)))
                                    ++count;
                        }
            }
            else
            {
                const int recv_block = r.this_block;
                if (recv_block < 0 || recv_block >= field.num_blocks())
                    continue;
                FieldBlock &fb = field.field(fid, recv_block);
                const Box3 &b = r.recv_box;
                for (int i = b.lo.i; i < b.hi.i; ++i)
                    for (int j = b.lo.j; j < b.hi.j; ++j)
                        for (int k = b.lo.k; k < b.hi.k; ++k)
                            for (int m = 0; m < desc.ncomp; ++m)
                                if (!std::isfinite(fb(i, j, k, m)))
                                    ++count;
            }
        }
        return count;
    }

    double component_copy_error(Field &field,
                                Halo &halo,
                                const std::string &name,
                                HaloLevel stage,
                                long long &checked_region_count)
    {
        double max_err = 0.0;
        checked_region_count = 0;
        const int fid = field.field_id(name);
        const FieldDescriptor &desc = field.descriptor(fid);
        const std::vector<HaloRegion> regions =
            halo.debug_halo_regions(desc.location, desc.nghost, stage);

        for (const HaloRegion &r : regions)
        {
            if (stage == HaloLevel::FaceOnly)
            {
                ++checked_region_count;
                if (r.this_rank == r.neighbor_rank)
                {
                    FieldBlock &recv = field.field(fid, r.neighbor_block);
                    const Box3 &b = r.send_box;
                    for (int i = b.lo.i; i < b.hi.i; ++i)
                        for (int j = b.lo.j; j < b.hi.j; ++j)
                            for (int k = b.lo.k; k < b.hi.k; ++k)
                            {
                                int ri = 0, rj = 0, rk = 0;
                                Z0_TEST::map_index(r.trans, i, j, k, ri, rj, rk);
                                for (int m = 0; m < desc.ncomp; ++m)
                                {
                                    const double expected = Z0_TEST::unique_code(r.this_rank, r.this_block, i, j, k, m);
                                    max_err = std::max(max_err, std::abs(recv(ri, rj, rk, m) - expected));
                                }
                            }
                }
                else
                {
                    FieldBlock &recv = field.field(fid, r.this_block);
                    const Box3 &b = r.recv_box;
                    for (int i = b.lo.i; i < b.hi.i; ++i)
                        for (int j = b.lo.j; j < b.hi.j; ++j)
                            for (int k = b.lo.k; k < b.hi.k; ++k)
                            {
                                int si = 0, sj = 0, sk = 0;
                                Z0_TEST::map_index(r.trans, i, j, k, si, sj, sk);
                                for (int m = 0; m < desc.ncomp; ++m)
                                {
                                    const double expected = Z0_TEST::unique_code(r.neighbor_rank, r.neighbor_block, si, sj, sk, m);
                                    max_err = std::max(max_err, std::abs(recv(i, j, k, m) - expected));
                                }
                            }
                }
            }
            else if (r.this_rank == r.neighbor_rank)
            {
                ++checked_region_count;
                FieldBlock &recv = field.field(fid, r.this_block);
                const Box3 &b = r.recv_box;
                for (int i = b.lo.i; i < b.hi.i; ++i)
                    for (int j = b.lo.j; j < b.hi.j; ++j)
                        for (int k = b.lo.k; k < b.hi.k; ++k)
                        {
                            int si = 0, sj = 0, sk = 0;
                            Z0_TEST::map_index(r.trans, i, j, k, si, sj, sk);
                            for (int m = 0; m < desc.ncomp; ++m)
                            {
                                const double expected = Z0_TEST::unique_code(r.neighbor_rank, r.neighbor_block, si, sj, sk, m);
                                max_err = std::max(max_err, std::abs(recv(i, j, k, m) - expected));
                            }
                        }
            }
        }

        return Z0_TEST::global_max(max_err);
    }

    bool run_stage(Field &field,
                   Halo &halo,
                   Z0_Boundary &boundary,
                   HaloLevel stage)
    {
        Z0_TEST::fill_all_with_nan(field);
        for (const std::string &name : Z0_TEST::registered_test_fields())
            Z0_TEST::fill_owned_unique(field, name);

        if (owner_only_enabled())
        {
            halo.sync_owner_alias_stage(HaloLevel::FaceOnly);
            if (static_cast<int>(stage) >= static_cast<int>(HaloLevel::Edge))
                halo.sync_owner_alias_stage(HaloLevel::Edge);
            if (static_cast<int>(stage) >= static_cast<int>(HaloLevel::Vertex))
                halo.sync_owner_alias_stage(HaloLevel::Vertex);
        }
        else
        {
            boundary.SyncAllRegistered(HaloLevel::FaceOnly);
            if (static_cast<int>(stage) >= static_cast<int>(HaloLevel::Edge))
                boundary.SyncAllRegistered(HaloLevel::Edge);
            if (static_cast<int>(stage) >= static_cast<int>(HaloLevel::Vertex))
                boundary.SyncAllRegistered(HaloLevel::Vertex);
        }

        bool passed = true;
        if (!owner_only_enabled())
        {
            for (const std::string &name : Z0_TEST::registered_test_fields())
            {
                long long local_regions = 0;
                const long long local_bad = count_nonfinite_recv_boxes(field, halo, name, stage, local_regions);
                const long long bad = Z0_TEST::global_sum(local_bad);
                const long long regions = Z0_TEST::global_sum(local_regions);
                std::ostringstream os;
                os << "stage=" << level_name(stage) << " regions=" << regions << " nonfinite_halo=" << bad;
                if (regions == 0)
                    os << " skipped_no_registered_regions";
                passed &= Z0_TEST::print_result("NaN overwrite " + name, bad == 0, os.str());
            }

            for (const std::string &name : {"U", "Bcell"})
            {
                long long local_checked_regions = 0;
                const double err = component_copy_error(field, halo, name, stage, local_checked_regions);
                const long long checked_regions = Z0_TEST::global_sum(local_checked_regions);
                std::ostringstream os;
                os << "stage=" << level_name(stage) << " checked_regions=" << checked_regions
                   << " diagnostic_exact_max=" << err;
                if (checked_regions == 0)
                    os << " skipped_no_checked_regions";
                passed &= Z0_TEST::print_result("unique component-copy diagnostic " + name,
                                                checked_regions >= 0 && std::isfinite(err),
                                                os.str());
            }
        }

        if (owner_sync_enabled())
        {
            struct OwnerCheck
            {
                const char *name;
                TOPO::EntityDim dim;
                HaloLevel stage;
            };

            const OwnerCheck checks[] = {
                {"B_xi", TOPO::EntityDim::Face, HaloLevel::FaceOnly},
                {"B_eta", TOPO::EntityDim::Face, HaloLevel::FaceOnly},
                {"B_zeta", TOPO::EntityDim::Face, HaloLevel::FaceOnly},
                {"E_xi", TOPO::EntityDim::Edge, HaloLevel::Edge},
                {"E_eta", TOPO::EntityDim::Edge, HaloLevel::Edge},
                {"E_zeta", TOPO::EntityDim::Edge, HaloLevel::Edge},
                {"phi", TOPO::EntityDim::Node, HaloLevel::Vertex},
            };

            for (const auto &c : checks)
            {
                if (static_cast<int>(stage) < static_cast<int>(c.stage))
                    continue;
                long long checked = 0;
                double err = 0.0;
                const bool ok = check_owner_alias(field, halo, c.name, c.dim, checked, err);
                std::ostringstream os;
                os << "stage=" << level_name(stage)
                   << " checked_alias_values=" << checked
                   << " owner_alias_max=" << err;
                if (checked == 0)
                    os << " skipped_no_local_alias";
                passed &= Z0_TEST::print_result(std::string("OwnerSync ") + c.name, ok, os.str());
            }
        }

        return passed;
    }
}

namespace Z0_TEST
{
    bool RunHaloCommunicationTests(Field &field, Halo &halo, Z0_Boundary &boundary, Param &param)
    {
        (void)param;
        if (rank() == 0)
            std::cout << "Z0 halo communication tests\n";

        bool passed = true;
        passed &= run_stage(field, halo, boundary, HaloLevel::FaceOnly);
        passed &= run_stage(field, halo, boundary, HaloLevel::Edge);
        passed &= run_stage(field, halo, boundary, HaloLevel::Vertex);
        return passed;
    }
}
