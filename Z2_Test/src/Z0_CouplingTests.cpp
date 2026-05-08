#include "Z0_CouplingTests.h"

#include "Z0_Initializer.h"
#include "0_basic/LayoutTraits.h"
#include "2_topology/TopologyView.h"
#include "3_field/Field.h"
#include "4_halo/Halo.h"
#include "4_halo/detail/HaloBuildTools.h"

#include <cmath>
#include <iostream>
#include <vector>

namespace
{
    constexpr double Tol = 1.0e-10;

    int axis_from_location(StaggerLocation loc)
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
            return 0;
        }
    }

    std::string edge_name(int axis)
    {
        return axis == 0 ? "E_xi" : (axis == 1 ? "E_eta" : "E_zeta");
    }

    std::string face_name(int axis)
    {
        return axis == 0 ? "B_xi" : (axis == 1 ? "B_eta" : "B_zeta");
    }

    int channel_id(const CouplingPairDesc &desc, const std::string &tag)
    {
        for (int cid = 0; cid < static_cast<int>(desc.channels.size()); ++cid)
            if (desc.channels[cid].tag == tag)
                return cid;
        return -1;
    }

    bool pair_has_patch(const TOPO::Topology &topology, const std::string &src, const std::string &dst)
    {
        for (const auto &p : topology.inner_patches)
            if (p.is_coupling && p.nb_block_name == src && p.this_block_name == dst)
                return true;
        for (const auto &p : topology.parallel_patches)
            if (p.is_coupling && p.nb_block_name == src && p.this_block_name == dst)
                return true;
        for (const auto &p : topology.inner_edge_patches)
            if (p.is_coupling && p.nb_block_name == src && p.this_block_name == dst)
                return true;
        for (const auto &p : topology.parallel_edge_patches)
            if (p.is_coupling && p.nb_block_name == src && p.this_block_name == dst)
                return true;
        for (const auto &p : topology.inner_vertex_patches)
            if (p.is_coupling && p.nb_block_name == src && p.this_block_name == dst)
                return true;
        for (const auto &p : topology.parallel_vertex_patches)
            if (p.is_coupling && p.nb_block_name == src && p.this_block_name == dst)
                return true;
        return false;
    }

    Field::PairKey choose_pair(Field &fields, const TOPO::Topology &topology)
    {
        for (const auto &kv : fields.coupling_pairs())
            if (pair_has_patch(topology, kv.first.first, kv.first.second))
                return kv.first;
        return Field::PairKey{"", ""};
    }

    int form_sign(const CouplingChannelSpec &ch, const TOPO::IndexTransform &tr, int dst_axis)
    {
        if (ch.value_kind == FieldValueKind::EdgeCovariant1Form)
            return tr.sign[dst_axis] >= 0 ? +1 : -1;
        const int src_axis = tr.perm[dst_axis];
        const TOPO::IndexTransform src_to_dst = HALO_TOOLS::inverse_transform(tr);
        return HALO_TOOLS::face_2form_orientation_sign(src_to_dst, src_axis);
    }

    std::string source_field_for_channel(const CouplingChannelSpec &ch,
                                         const TOPO::IndexTransform &tr,
                                         bool form_aware)
    {
        if (!form_aware)
            return ch.tag;
        const int dst_axis = axis_from_location(ch.location);
        const int src_axis = tr.perm[dst_axis];
        if (ch.value_kind == FieldValueKind::EdgeCovariant1Form)
            return edge_name(src_axis);
        return face_name(src_axis);
    }

    template <class PatchT>
    void check_buffers(Field &fields,
                       const std::vector<std::vector<CouplingBufferBlock>> &storage,
                       const std::vector<PatchT> &patches,
                       const CouplingPairDesc &pair_desc,
                       int cid,
                       int my_rank,
                       const char *test_name,
                       bool form_aware,
                       Z0::TestResult &result,
                       bool &nontrivial,
                       std::ostream &os)
    {
        if (cid < 0 || cid >= static_cast<int>(storage.size()))
            return;
        const CouplingChannelSpec &ch = pair_desc.channels[cid];
        int failures = 0;

        for (std::size_t ip = 0; ip < patches.size(); ++ip)
        {
            const PatchT &p = patches[ip];
            if (!p.is_coupling || p.this_block_name != pair_desc.pair.dst || p.nb_block_name != pair_desc.pair.src)
                continue;
            if (ip >= storage[cid].size())
                continue;
            const CouplingBufferBlock &buf = storage[cid][ip];
            if (!buf.allocated)
                continue;

            const int dst_axis = axis_from_location(ch.location);
            const int src_axis = form_aware ? p.trans.perm[dst_axis] : dst_axis;
            int sign = 1;
            if (form_aware)
                sign = form_sign(ch, p.trans, dst_axis);
            nontrivial = nontrivial || (form_aware && (src_axis != dst_axis || sign < 0));

            const std::string src_name = source_field_for_channel(ch, p.trans, form_aware);
            const int fid_src = fields.field_id(src_name);
            const FieldDescriptor &dst_desc = fields.descriptor(ch.tag);

            const Box3 &b = buf.box;
            for (int i = b.lo.i; i < b.hi.i; ++i)
                for (int j = b.lo.j; j < b.hi.j; ++j)
                    for (int k = b.lo.k; k < b.hi.k; ++k)
                    {
                        int is, js, ks;
                        HALO_TOOLS::apply_transform(p.trans, i, j, k, is, js, ks);
                        for (int m = 0; m < ch.ncomp; ++m)
                        {
                            const double expected = static_cast<double>(sign) *
                                                    Z0::analytic_value_for_rank(fid_src, m, p.nb_rank, p.nb_block, is, js, ks);
                            const double actual = const_cast<CouplingBufferBlock &>(buf)(i, j, k, m);
                            const double err = std::abs(actual - expected);
                            result.max_error = std::max(result.max_error, err);
                            if (actual < -0.5e300 || err > Tol)
                            {
                                result.pass = false;
                                if (failures++ < 5)
                                    Z0::report_failure(test_name, fields, dst_desc, my_rank, p.this_block, i, j, k, m,
                                                       expected, actual, os, &p.trans);
                            }
                        }
                    }
        }
    }

    Z0::TestResult run_coupling_test(Field &fields,
                                     Halo &halo,
                                     const TOPO::Topology &topology,
                                     int dimension,
                                     int my_rank,
                                     std::ostream &os,
                                     const std::string &tag,
                                     const char *test_name,
                                     bool form_aware)
    {
        Z0::TestResult result;
        const Field::PairKey pair = choose_pair(fields, topology);
        if (pair.first.empty())
        {
            if (my_rank == 0)
                os << "[Z0][" << test_name << "] SKIP no coupling pair in current CASE\n";
            return result;
        }

        std::string src = pair.first;
        std::string dst = pair.second;
        halo.coupling_trans_1DCorner(src, dst);
        if (dimension >= 2)
            halo.coupling_trans_2DCorner(src, dst);
        if (dimension >= 3)
            halo.coupling_trans_3DCorner(src, dst);

        CouplingBuffersForPair &cb = fields.coupling_buffers(src, dst);
        const int cid = channel_id(cb.desc, tag);
        if (cid < 0)
        {
            result.pass = false;
            os << "[Z0][" << test_name << "] missing coupling channel tag=" << tag << "\n";
            return result;
        }

        bool nontrivial = false;
        const auto inner_faces = TOPO_VIEW::inner_faces(topology);
        const auto parallel_faces = TOPO_VIEW::parallel_faces(topology);
        const auto &inner_edges = TOPO_VIEW::edge_patches(topology, TOPO::PatchKind::Inner);
        const auto &parallel_edges = TOPO_VIEW::edge_patches(topology, TOPO::PatchKind::Parallel);
        const auto &inner_vertices = TOPO_VIEW::vertex_patches(topology, TOPO::PatchKind::Inner);
        const auto &parallel_vertices = TOPO_VIEW::vertex_patches(topology, TOPO::PatchKind::Parallel);

        check_buffers(fields, cb.inner_face, inner_faces, cb.desc, cid, my_rank, test_name, form_aware, result, nontrivial, os);
        check_buffers(fields, cb.parallel_face, parallel_faces, cb.desc, cid, my_rank, test_name, form_aware, result, nontrivial, os);
        check_buffers(fields, cb.inner_edge, inner_edges, cb.desc, cid, my_rank, test_name, form_aware, result, nontrivial, os);
        check_buffers(fields, cb.parallel_edge, parallel_edges, cb.desc, cid, my_rank, test_name, form_aware, result, nontrivial, os);
        check_buffers(fields, cb.inner_vertex, inner_vertices, cb.desc, cid, my_rank, test_name, form_aware, result, nontrivial, os);
        check_buffers(fields, cb.parallel_vertex, parallel_vertices, cb.desc, cid, my_rank, test_name, form_aware, result, nontrivial, os);

        if (form_aware && !nontrivial && my_rank == 0)
            os << "No nontrivial orientation transform found in this CASE; test covered identity transform only.\n";
        Z0::report_test(test_name, result, os);
        return result;
    }
}

namespace Z0
{
    TestResult test_coupling_cell_scalar(Field &fields, Halo &halo, const TOPO::Topology &topology, int dimension, int my_rank, std::ostream &os)
    {
        return run_coupling_test(fields, halo, topology, dimension, my_rank, os, "phi_cell", "CouplingCellScalar", false);
    }

    TestResult test_coupling_edge_1form(Field &fields, Halo &halo, const TOPO::Topology &topology, int dimension, int my_rank, std::ostream &os)
    {
        return run_coupling_test(fields, halo, topology, dimension, my_rank, os, "E_xi", "CouplingEdge1Form", true);
    }

    TestResult test_coupling_face_2form(Field &fields, Halo &halo, const TOPO::Topology &topology, int dimension, int my_rank, std::ostream &os)
    {
        return run_coupling_test(fields, halo, topology, dimension, my_rank, os, "B_xi", "CouplingFace2Form", true);
    }
}
