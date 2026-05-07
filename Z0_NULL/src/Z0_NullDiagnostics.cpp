#include "Z0_NullDiagnostics.h"

#include "0_basic/LayoutTraits.h"
#include "3_field/2_MPCNS_Field.h"
#include "4_halo/1_MPCNS_Halo.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

namespace
{
    bool member_matches_location(const TOPO::EquivMember &m, StaggerLocation loc)
    {
        return m.location == loc;
    }

    const char *null_halo_level_name(HaloLevel level)
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

    double member_value(const Field &fields,
                        int fid,
                        const TOPO::EquivMember &m,
                        int comp)
    {
        Field &mutable_fields = const_cast<Field &>(fields);
        return mutable_fields.field(fid, m.block)(m.i, m.j, m.k, comp);
    }

    void summarize_field(const Field &fields,
                         int fid,
                         int ib,
                         double sentinel,
                         std::ostream &os)
    {
        Field &mutable_fields = const_cast<Field &>(fields);
        FieldBlock &fb = mutable_fields.field(fid, ib);
        if (!fb.is_allocated())
        {
            os << "    block " << ib << " allocated=false\n";
            return;
        }

        const FieldDescriptor &desc = fields.descriptor(fid);
        const Int3 lo = fb.get_lo();
        const Int3 hi = fb.get_hi();

        double vmin = std::numeric_limits<double>::infinity();
        double vmax = -std::numeric_limits<double>::infinity();
        int64_t nan_count = 0;
        int64_t inf_count = 0;
        int64_t sentinel_count = 0;
        int64_t finite_count = 0;

        for (int i = lo.i; i < hi.i; ++i)
            for (int j = lo.j; j < hi.j; ++j)
                for (int k = lo.k; k < hi.k; ++k)
                    for (int m = 0; m < desc.ncomp; ++m)
                    {
                        const double v = fb(i, j, k, m);
                        if (std::isnan(v))
                        {
                            ++nan_count;
                            continue;
                        }
                        if (std::isinf(v))
                        {
                            ++inf_count;
                            continue;
                        }
                        if (v == sentinel)
                            ++sentinel_count;
                        vmin = std::min(vmin, v);
                        vmax = std::max(vmax, v);
                        ++finite_count;
                    }

        if (finite_count == 0)
        {
            vmin = 0.0;
            vmax = 0.0;
        }

        os << "    block " << ib
           << " allocated=true"
           << " min=" << vmin
           << " max=" << vmax
           << " nan_count=" << nan_count
           << " inf_count=" << inf_count
           << " sentinel_count=" << sentinel_count
           << "\n";
    }
}

namespace Z0_NULL
{
    void print_banner()
    {
        std::cout << "========== Z0_NULL Physics-Null Framework Template ==========\n"
                  << "This program demonstrates:\n"
                  << "  1. Register representative FieldDescriptor.\n"
                  << "  2. Initialize fields through user-replaceable functions.\n"
                  << "  3. Build Grid / Topology / TopologyEquiv / Halo.\n"
                  << "  4. Validate halo synchronization and owner alias sync.\n"
                  << "  5. Optionally write Tecplot output through IOModule.\n"
                  << "No numerical PDE solver is executed.\n"
                  << "============================================================\n"
                  << std::flush;
    }

    void print_diagnostics(const Field &fields,
                           const TOPO::TopologyEquiv &topology_equiv,
                           int dimension,
                           int nghost)
    {
        std::cout << "---------- Z0_NULL Diagnostics ----------\n"
                  << "dimension              = " << dimension << "\n"
                  << "ngg                    = " << nghost << "\n"
                  << "local blocks           = " << fields.num_blocks() << "\n"
                  << "registered fields      = " << fields.num_fields() << "\n"
                  << "halo requests          = " << fields.halo_requests().size() << "\n"
                  << "edge owner classes     = " << topology_equiv.edge_classes_general.size() << "\n"
                  << "face owner classes     = " << topology_equiv.face_classes.size() << "\n"
                  << "local edge owners      = " << topology_equiv.n_local_edge_owner << "\n"
                  << "global edge owners     = " << topology_equiv.n_global_edge_owner << "\n"
                  << "local face owners      = " << topology_equiv.n_local_face_owner << "\n"
                  << "global face owners     = " << topology_equiv.n_global_face_owner << "\n"
                  << std::flush;
    }

    void dump_halo_registry_if_requested(const Halo &halo, bool dump_registry)
    {
        if (dump_registry)
            halo.dump_sync_registry(std::cout);
    }

    namespace DIAG
    {
        void dump_field_catalog(const Field &fields,
                                int my_rank,
                                std::ostream &os)
        {
            os << "[rank " << my_rank << "] Field catalog:\n";
            int fid = 0;
            for (const auto &desc : fields.descriptors())
            {
                os << "  fid=" << fid++
                   << " name=" << desc.name
                   << " loc=" << LAYOUT::location_name(desc.location)
                   << " kind=" << field_value_kind_name(desc.value_kind)
                   << " ncomp=" << desc.ncomp
                   << " nghost=" << desc.nghost
                   << " group=" << desc.sync.group
                   << " do_halo=" << (desc.sync.do_halo ? "true" : "false")
                   << " halo_level=" << null_halo_level_name(desc.sync.halo_level)
                   << " orientation_aware=" << (desc.sync.orientation_aware ? "true" : "false")
                   << " owner_sync=" << owner_sync_policy_name(desc.sync.owner_sync)
                   << "\n";
            }
        }

        void dump_field_block_summary(const Field &fields,
                                      int my_rank,
                                      std::ostream &os)
        {
            constexpr double sentinel = -9.87654321e200;
            os << "[rank " << my_rank << "] Field block summary:\n";
            for (int fid = 0; fid < fields.num_fields(); ++fid)
            {
                os << "  field=" << fields.descriptor(fid).name << "\n";
                for (int ib = 0; ib < fields.num_blocks(); ++ib)
                    summarize_field(fields, fid, ib, sentinel, os);
            }
        }

        void dump_topology_equiv_summary(const TOPO::TopologyEquiv &equiv,
                                         int my_rank,
                                         std::ostream &os)
        {
            os << "[rank " << my_rank << "] TopologyEquiv summary:\n"
               << "  edge_classes=" << equiv.classes(TOPO::EquivDofKind::Edge).size() << "\n"
               << "  face_classes=" << equiv.classes(TOPO::EquivDofKind::Face).size() << "\n"
               << "  local_edge_owners=" << equiv.n_local_edge_owner << "\n"
               << "  global_edge_owners=" << equiv.n_global_edge_owner << "\n"
               << "  local_face_owners=" << equiv.n_local_face_owner << "\n"
               << "  global_face_owners=" << equiv.n_global_face_owner << "\n";
        }

        bool check_field_finite(const Field &fields,
                                const std::string &field_name,
                                int my_rank,
                                std::ostream &os)
        {
            if (!fields.has_field(field_name))
            {
                os << "[rank " << my_rank << "] finite check SKIP missing field=" << field_name << "\n";
                return true;
            }

            Field &mutable_fields = const_cast<Field &>(fields);
            const int fid = fields.field_id(field_name);
            const FieldDescriptor &desc = fields.descriptor(fid);
            int64_t bad_count = 0;

            for (int ib = 0; ib < fields.num_blocks(); ++ib)
            {
                FieldBlock &fb = mutable_fields.field(fid, ib);
                if (!fb.is_allocated())
                    continue;

                const Int3 lo = fb.get_lo();
                const Int3 hi = fb.get_hi();
                for (int i = lo.i; i < hi.i; ++i)
                    for (int j = lo.j; j < hi.j; ++j)
                        for (int k = lo.k; k < hi.k; ++k)
                            for (int m = 0; m < desc.ncomp; ++m)
                                if (!std::isfinite(fb(i, j, k, m)))
                                    ++bad_count;
            }

            os << "[rank " << my_rank << "] finite check field=" << field_name
               << " bad_count=" << bad_count
               << (bad_count == 0 ? " PASS\n" : " FAIL\n");
            return bad_count == 0;
        }

        double check_edge_owner_alias_error(const Field &fields,
                                            const TOPO::TopologyEquiv &equiv,
                                            const std::string &field_name,
                                            int my_rank,
                                            std::ostream &os)
        {
            if (!fields.has_field(field_name))
            {
                os << "[rank " << my_rank << "] edge owner alias SKIP missing field=" << field_name << "\n";
                return 0.0;
            }

            const int fid = fields.field_id(field_name);
            const FieldDescriptor &desc = fields.descriptor(fid);
            double max_err = 0.0;
            int64_t checked = 0;

            for (const auto &cls : equiv.classes(TOPO::EquivDofKind::Edge))
            {
                const auto &owner = cls.owner;
                if (owner.rank != my_rank || !member_matches_location(owner, desc.location))
                    continue;

                for (const auto &alias : cls.members)
                {
                    if (alias.is_owner || alias.rank != my_rank || !member_matches_location(alias, desc.location))
                        continue;

                    const int sign = alias.orient_sign * owner.orient_sign;
                    for (int m = 0; m < desc.ncomp; ++m)
                    {
                        const double ov = member_value(fields, fid, owner, m);
                        const double av = member_value(fields, fid, alias, m);
                        max_err = std::max(max_err, std::abs(av - static_cast<double>(sign) * ov));
                        ++checked;
                    }
                }
            }

            os << "[rank " << my_rank << "] edge owner alias field=" << field_name
               << " checked=" << checked
               << " max_error=" << max_err << "\n";
            return max_err;
        }

        double check_face_owner_alias_error(const Field &fields,
                                            const TOPO::TopologyEquiv &equiv,
                                            const std::string &field_name,
                                            int my_rank,
                                            std::ostream &os)
        {
            if (!fields.has_field(field_name))
            {
                os << "[rank " << my_rank << "] face owner alias SKIP missing field=" << field_name << "\n";
                return 0.0;
            }

            const int fid = fields.field_id(field_name);
            const FieldDescriptor &desc = fields.descriptor(fid);
            double max_err = 0.0;
            int64_t checked = 0;

            for (const auto &cls : equiv.classes(TOPO::EquivDofKind::Face))
            {
                const auto &owner = cls.owner;
                if (owner.rank != my_rank || !member_matches_location(owner, desc.location))
                    continue;

                for (const auto &alias : cls.members)
                {
                    if (alias.is_owner || alias.rank != my_rank || !member_matches_location(alias, desc.location))
                        continue;

                    const int sign = alias.orient_sign * owner.orient_sign;
                    for (int m = 0; m < desc.ncomp; ++m)
                    {
                        const double ov = member_value(fields, fid, owner, m);
                        const double av = member_value(fields, fid, alias, m);
                        max_err = std::max(max_err, std::abs(av - static_cast<double>(sign) * ov));
                        ++checked;
                    }
                }
            }

            os << "[rank " << my_rank << "] face owner alias field=" << field_name
               << " checked=" << checked
               << " max_error=" << max_err << "\n";
            return max_err;
        }

        bool run_basic_halo_validation(Field &fields,
                                       Halo &halo,
                                       const TOPO::TopologyEquiv &equiv,
                                       int my_rank,
                                       std::ostream &os)
        {
            bool pass = true;
            auto sync_field_if_present = [&](const std::string &name)
            {
                if (!fields.has_field(name))
                {
                    os << "[rank " << my_rank << "] sync_field SKIP missing field=" << name << "\n";
                    return;
                }
                halo.sync_field(name);
                pass = check_field_finite(fields, name, my_rank, os) && pass;
            };

            sync_field_if_present("phi");
            sync_field_if_present("U");
            sync_field_if_present("V_cart");

            if (fields.has_field("E_xi") && fields.has_field("E_eta") && fields.has_field("E_zeta"))
            {
                halo.sync_group("Eedge");
                pass = check_edge_owner_alias_error(fields, equiv, "E_xi", my_rank, os) == 0.0 && pass;
                pass = check_edge_owner_alias_error(fields, equiv, "E_eta", my_rank, os) == 0.0 && pass;
                pass = check_edge_owner_alias_error(fields, equiv, "E_zeta", my_rank, os) == 0.0 && pass;
            }
            else
            {
                os << "[rank " << my_rank << "] sync_group SKIP missing group=Eedge\n";
            }

            if (fields.has_field("B_xi") && fields.has_field("B_eta") && fields.has_field("B_zeta"))
            {
                halo.sync_group("Bface");
                pass = check_face_owner_alias_error(fields, equiv, "B_xi", my_rank, os) == 0.0 && pass;
                pass = check_face_owner_alias_error(fields, equiv, "B_eta", my_rank, os) == 0.0 && pass;
                pass = check_face_owner_alias_error(fields, equiv, "B_zeta", my_rank, os) == 0.0 && pass;
            }
            else
            {
                os << "[rank " << my_rank << "] sync_group SKIP missing group=Bface\n";
            }

            os << "[rank " << my_rank << "] basic halo validation "
               << (pass ? "PASS\n" : "FAIL\n");
            return pass;
        }
    }
}
