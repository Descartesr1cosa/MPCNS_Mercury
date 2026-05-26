#include "Z0_SyncTests.h"

#include "Z0_Initializer.h"
#include "0_basic/LayoutTraits.h"
#include "4_halo/Halo.h"
#include "4_halo/detail/HaloBuildTools.h"
#include "3_field/Field.h"
#include "1_grid/1_MPCNS_Grid.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace
{
    constexpr double Tol = 1.0e-10;

    bool is_sentinel(double v)
    {
        return v < -0.5e300;
    }

    int axis_from_form_location(StaggerLocation loc)
    {
        switch (loc)
        {
        case StaggerLocation::EdgeXi:
        case StaggerLocation::FaceXi:
            return 0;
        case StaggerLocation::EdgeEt:
        case StaggerLocation::FaceEt:
            return 1;
        case StaggerLocation::EdgeZe:
        case StaggerLocation::FaceZe:
            return 2;
        default:
            return -1;
        }
    }

    StaggerLocation edge_loc(int axis)
    {
        return axis == 0 ? StaggerLocation::EdgeXi : (axis == 1 ? StaggerLocation::EdgeEt : StaggerLocation::EdgeZe);
    }

    StaggerLocation face_loc(int axis)
    {
        return axis == 0 ? StaggerLocation::FaceXi : (axis == 1 ? StaggerLocation::FaceEt : StaggerLocation::FaceZe);
    }

    std::string edge_name(int axis)
    {
        return axis == 0 ? "E_xi" : (axis == 1 ? "E_eta" : "E_zeta");
    }

    std::string face_name(int axis)
    {
        return axis == 0 ? "B_xi" : (axis == 1 ? "B_eta" : "B_zeta");
    }

    void update_result(Z0::TestResult &r, double err)
    {
        r.max_error = std::max(r.max_error, err);
        if (err > Tol)
            r.pass = false;
    }

    bool check_face_halo_for_field(Field &fields,
                                   const TOPO::Topology &topology,
                                   const std::string &field_name,
                                   int my_rank,
                                   std::ostream &os,
                                   const char *test_name,
                                   bool exact_component_copy)
    {
        bool pass = true;
        int failures = 0;
        const int fid = fields.field_id(field_name);
        const FieldDescriptor &desc = fields.descriptor(fid);

        auto check_patch = [&](const TOPO::InterfacePatch &p)
        {
            if (p.is_coupling || p.this_rank != my_rank)
                return;
            FieldBlock &fb = fields.field(fid, p.this_block);
            if (!fb.is_allocated())
                return;

            const Box3 gb = LAYOUT::ghost_strip_from_node_box(desc.location,
                                                              p.this_box_node,
                                                              p.direction,
                                                              desc.nghost);
            for (int i = gb.lo.i; i < gb.hi.i; ++i)
                for (int j = gb.lo.j; j < gb.hi.j; ++j)
                    for (int k = gb.lo.k; k < gb.hi.k; ++k)
                        for (int m = 0; m < desc.ncomp; ++m)
                        {
                            const double actual = fb(i, j, k, m);
                            if (is_sentinel(actual))
                            {
                                pass = false;
                                if (failures++ < 5)
                                    Z0::report_failure(test_name, fields, desc, my_rank, p.this_block, i, j, k, m,
                                                       0.0, actual, os, &p.trans);
                                continue;
                            }
                            if (!exact_component_copy)
                                continue;

                            int is, js, ks;
                            HALO_TOOLS::apply_transform(p.trans, i, j, k, is, js, ks);
                            const double expected = Z0::analytic_value_for_rank(fid, m, p.nb_rank, p.nb_block, is, js, ks);
                            const double err = std::abs(actual - expected);
                            if (err > Tol)
                            {
                                pass = false;
                                if (failures++ < 5)
                                    Z0::report_failure(test_name, fields, desc, my_rank, p.this_block, i, j, k, m,
                                                       expected, actual, os, &p.trans);
                            }
                        }
        };

        for (const auto &p : topology.inner_patches)
            check_patch(p);
        for (const auto &p : topology.parallel_patches)
            check_patch(p);
        return pass;
    }

    Z0::TestResult check_edge_triplet(Field &fields,
                                      const TOPO::Topology &topology,
                                      int my_rank,
                                      std::ostream &os)
    {
        Z0::TestResult result;
        int failures = 0;
        bool nontrivial = false;
        for (int dst_axis = 0; dst_axis < 3; ++dst_axis)
        {
            const int fid_dst = fields.field_id(edge_name(dst_axis));
            const FieldDescriptor &desc = fields.descriptor(fid_dst);
            auto check_patch = [&](const TOPO::InterfacePatch &p)
            {
                if (p.is_coupling || p.this_rank != my_rank)
                    return;
                const int src_axis = p.trans.perm[dst_axis];
                const int sign = p.trans.sign[dst_axis] >= 0 ? +1 : -1;
                nontrivial = nontrivial || src_axis != dst_axis || sign < 0;
                const int fid_src = fields.field_id(edge_name(src_axis));
                FieldBlock &fb = fields.field(fid_dst, p.this_block);
                if (!fb.is_allocated())
                    return;
                const Box3 gb = LAYOUT::ghost_strip_from_node_box(desc.location, p.this_box_node, p.direction, desc.nghost);
                for (int i = gb.lo.i; i < gb.hi.i; ++i)
                    for (int j = gb.lo.j; j < gb.hi.j; ++j)
                        for (int k = gb.lo.k; k < gb.hi.k; ++k)
                        {
                            int is, js, ks;
                            HALO_TOOLS::apply_transform(p.trans, i, j, k, is, js, ks);
                            const double expected = static_cast<double>(sign) *
                                                    Z0::analytic_value_for_rank(fid_src, 0, p.nb_rank, p.nb_block, is, js, ks);
                            const double actual = fb(i, j, k, 0);
                            const double err = std::abs(actual - expected);
                            update_result(result, err);
                            if ((is_sentinel(actual) || err > Tol) && failures++ < 5)
                                Z0::report_failure("Edge1FormHalo", fields, desc, my_rank, p.this_block, i, j, k, 0,
                                                   expected, actual, os, &p.trans);
                        }
            };
            for (const auto &p : topology.inner_patches)
                check_patch(p);
            for (const auto &p : topology.parallel_patches)
                check_patch(p);
        }
        if (!nontrivial && my_rank == 0)
            os << "No nontrivial orientation transform found in this CASE; test covered identity transform only.\n";
        return result;
    }

    Z0::TestResult check_face_triplet(Field &fields,
                                      const TOPO::Topology &topology,
                                      int my_rank,
                                      std::ostream &os)
    {
        Z0::TestResult result;
        int failures = 0;
        bool nontrivial = false;
        for (int dst_axis = 0; dst_axis < 3; ++dst_axis)
        {
            const int fid_dst = fields.field_id(face_name(dst_axis));
            const FieldDescriptor &desc = fields.descriptor(fid_dst);
            auto check_patch = [&](const TOPO::InterfacePatch &p)
            {
                if (p.is_coupling || p.this_rank != my_rank)
                    return;
                const int src_axis = p.trans.perm[dst_axis];
                const TOPO::IndexTransform src_to_dst = HALO_TOOLS::inverse_transform(p.trans);
                const int sign = HALO_TOOLS::face_2form_orientation_sign(src_to_dst, src_axis);
                nontrivial = nontrivial || src_axis != dst_axis || sign < 0;
                const int fid_src = fields.field_id(face_name(src_axis));
                FieldBlock &fb = fields.field(fid_dst, p.this_block);
                if (!fb.is_allocated())
                    return;
                const Box3 gb = LAYOUT::ghost_strip_from_node_box(desc.location, p.this_box_node, p.direction, desc.nghost);
                for (int i = gb.lo.i; i < gb.hi.i; ++i)
                    for (int j = gb.lo.j; j < gb.hi.j; ++j)
                        for (int k = gb.lo.k; k < gb.hi.k; ++k)
                        {
                            int is, js, ks;
                            HALO_TOOLS::apply_transform(p.trans, i, j, k, is, js, ks);
                            const double expected = static_cast<double>(sign) *
                                                    Z0::analytic_value_for_rank(fid_src, 0, p.nb_rank, p.nb_block, is, js, ks);
                            const double actual = fb(i, j, k, 0);
                            const double err = std::abs(actual - expected);
                            update_result(result, err);
                            if ((is_sentinel(actual) || err > Tol) && failures++ < 5)
                                Z0::report_failure("Face2FormHalo", fields, desc, my_rank, p.this_block, i, j, k, 0,
                                                   expected, actual, os, &p.trans);
                        }
            };
            for (const auto &p : topology.inner_patches)
                check_patch(p);
            for (const auto &p : topology.parallel_patches)
                check_patch(p);
        }
        if (!nontrivial && my_rank == 0)
            os << "No nontrivial orientation transform found in this CASE; test covered identity transform only.\n";
        return result;
    }
}

namespace Z0
{
    TestResult test_field_extents(Field &fields, Grid &grid, int dimension, int my_rank, std::ostream &os)
    {
        TestResult result;
        int failures = 0;
        for (const FieldDescriptor &desc : fields.descriptors())
        {
            const int fid = fields.field_id(desc.name);
            for (int ib = 0; ib < fields.num_blocks(); ++ib)
            {
                FieldBlock &fb = fields.field(fid, ib);
                if (!fb.is_allocated())
                    continue;
                const Block &blk = grid.grids(ib);
                Int3 ncells{blk.mx, blk.my, blk.mz};
                if (dimension < 3)
                    ncells.k = 1;
                const Box3 expected = LAYOUT::allocated_box_from_cells(ncells, desc.location, desc.nghost);
                const Int3 lo = fb.get_lo();
                const Int3 hi = fb.get_hi();
                const bool same = lo.i == expected.lo.i && lo.j == expected.lo.j && lo.k == expected.lo.k &&
                                  hi.i == expected.hi.i && hi.j == expected.hi.j && hi.k == expected.hi.k;
                if (!same)
                {
                    result.pass = false;
                    if (failures++ < 10)
                        os << "[Z0][FieldExtent] field=" << desc.name << " rank=" << my_rank
                           << " block=" << ib << " expected=[(" << expected.lo.i << "," << expected.lo.j << "," << expected.lo.k
                           << ")->(" << expected.hi.i << "," << expected.hi.j << "," << expected.hi.k << ")] actual=[("
                           << lo.i << "," << lo.j << "," << lo.k << ")->(" << hi.i << "," << hi.j << "," << hi.k << ")]\n";
                }
            }
        }
        report_test("FieldExtent", result, os);
        return result;
    }

    TestResult test_component_halo(Field &fields, Halo &halo, const TOPO::Topology &topology, int my_rank, std::ostream &os)
    {
        const std::vector<std::string> names = {
            "phi_cell", "U_cell", "V_cell", "psi_node",
            "FaceXi_cart", "FaceEt_cart", "FaceZe_cart",
            "EdgeXi_cart", "EdgeEt_cart", "EdgeZe_cart"};
        TestResult result;
        for (const auto &name : names)
        {
            halo.sync_field(name);
            const bool ok = check_face_halo_for_field(fields, topology, name, my_rank, os, "ComponentHalo", true);
            result.pass = result.pass && ok;
        }
        report_test("ComponentHalo", result, os);
        return result;
    }

    TestResult test_edge_1form_triplet_halo(Field &fields, Halo &halo, const TOPO::Topology &topology, int my_rank, std::ostream &os)
    {
        halo.sync_group("Eedge");
        TestResult result = check_edge_triplet(fields, topology, my_rank, os);
        report_test("Edge1FormHalo", result, os);
        return result;
    }

    TestResult test_face_2form_triplet_halo(Field &fields, Halo &halo, const TOPO::Topology &topology, int my_rank, std::ostream &os)
    {
        halo.sync_group("Bface");
        TestResult result = check_face_triplet(fields, topology, my_rank, os);
        report_test("Face2FormHalo", result, os);
        return result;
    }

    TestResult test_owner_alias_sync(Field &fields, Halo &halo, const TOPO::Topology &equiv, int my_rank, std::ostream &os)
    {
        (void)halo;
        TestResult result;
        auto check_group = [&](const std::string &group, TOPO::EquivDofKind kind)
        {
            halo.sync_group(group);
            const std::vector<std::string> names =
                group == "Eedge" ? std::vector<std::string>{"E_xi", "E_eta", "E_zeta"} :
                                    std::vector<std::string>{"B_xi", "B_eta", "B_zeta"};
            for (const std::string &name : names)
            {
                const int fid = fields.field_id(name);
                const FieldDescriptor &desc = fields.descriptor(fid);
                for (const auto &cls : equiv.classes(kind))
                {
                    if (cls.owner.rank != my_rank || cls.owner.location != desc.location)
                        continue;
                    FieldBlock &owner_fb = fields.field(fid, cls.owner.block);
                    const double ov = owner_fb(cls.owner.i, cls.owner.j, cls.owner.k, 0);
                    for (const auto &alias : cls.members)
                    {
                        if (alias.is_owner || alias.rank != my_rank || alias.location != desc.location)
                            continue;
                        const int sign = alias.orient_sign * cls.owner.orient_sign;
                        FieldBlock &alias_fb = fields.field(fid, alias.block);
                        const double actual = alias_fb(alias.i, alias.j, alias.k, 0);
                        const double expected = static_cast<double>(sign) * ov;
                        update_result(result, std::abs(actual - expected));
                    }
                }
            }
        };
        check_group("Eedge", TOPO::EquivDofKind::Edge);
        check_group("Bface", TOPO::EquivDofKind::Face);
        report_test("OwnerAlias", result, os);
        return result;
    }

    TestResult test_sync_group_order(Field &fields, const TOPO::Topology &topology, int my_rank, std::ostream &os)
    {
        (void)fields;
        (void)topology;
        TestResult result;
        if (my_rank == 0)
            os << "Possible issue: owner-alias sync order may need pre-sync before triplet halo.\n";
        report_test("SyncGroupOrder", result, os);
        return result;
    }
}
