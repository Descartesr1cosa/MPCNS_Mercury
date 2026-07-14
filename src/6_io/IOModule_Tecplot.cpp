#include "6_io/IOModule.h"
#include "7_metric/Metric.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <set>
#include <stdexcept>
#include <unordered_set>

void IOModule::TecWriteI32_(FILE *fp, int32_t v) { std::fwrite(&v, 4, 1, fp); }

void IOModule::TecWriteF32_(FILE *fp, float v) { std::fwrite(&v, 4, 1, fp); }

void IOModule::TecWriteF64_(FILE *fp, double v) { std::fwrite(&v, 8, 1, fp); }

void IOModule::TecWriteStr_(FILE *fp, const std::string &s)
{
    // Tecplot binary string: each char written as int32, terminated by '\0'
    for (char c : s)
    {
        int32_t v = static_cast<int32_t>(c);
        std::fwrite(&v, sizeof(int32_t), 1, fp);
    }
    int32_t z = 0;
    std::fwrite(&z, sizeof(int32_t), 1, fp);
}

IOModule::TecplotMode IOModule::NormalizedTecplotMode_() const
{
    if (tec_mode_ == TecplotMode::CellToNode)
        return TecplotMode::AllNode;
    if (tec_mode_ == TecplotMode::Mixed)
        return TecplotMode::CellAndNode;
    return tec_mode_;
}

bool IOModule::TecplotBlockSelected_(const std::string &block_name) const
{
    return tec_block_.empty() ||
           std::find(tec_block_.begin(), tec_block_.end(), block_name) != tec_block_.end();
}

bool IOModule::BuildTecFormTriplet_(const std::string &name,
                                    int &fid_xi,
                                    int &fid_eta,
                                    int &fid_zeta,
                                    bool &is_face,
                                    std::string &output_name) const
{
    fid_xi = -1;
    fid_eta = -1;
    fid_zeta = -1;
    is_face = false;
    output_name = name;

    std::string group_name = name;
    if (fld_->has_field(name))
    {
        const FieldDescriptor &direct = fld_->descriptor(name);
        const bool direct_form =
            direct.value_kind == FieldValueKind::FaceContravariant2Form ||
            direct.value_kind == FieldValueKind::EdgeCovariant1Form;
        if (direct_form && !direct.sync.group.empty())
        {
            group_name = direct.sync.group;
            output_name = direct.sync.group;
        }
    }

    auto accept = [&](const FieldDescriptor &d) -> bool
    {
        if (d.sync.group == group_name)
            return true;
        if (d.name == group_name)
            return true;
        return false;
    };

    for (const FieldDescriptor &d : fld_->descriptors())
    {
        if (!accept(d))
            continue;

        const bool face = (d.value_kind == FieldValueKind::FaceContravariant2Form);
        const bool edge = (d.value_kind == FieldValueKind::EdgeCovariant1Form);
        if (!face && !edge)
            continue;

        if (d.ncomp != 1)
            continue;

        if (d.sync.group.empty() && d.name != group_name)
            continue;

        if (fid_xi >= 0 && is_face != face)
            return false;

        is_face = face;
        if (!d.sync.group.empty())
            output_name = d.sync.group;

        const int fid = fld_->field_id(d.name);
        if (d.location == StaggerLocation::FaceXi || d.location == StaggerLocation::EdgeXi)
            fid_xi = fid;
        else if (d.location == StaggerLocation::FaceEt || d.location == StaggerLocation::EdgeEt)
            fid_eta = fid;
        else if (d.location == StaggerLocation::FaceZe || d.location == StaggerLocation::EdgeZe)
            fid_zeta = fid;
    }

    return fid_xi >= 0 && fid_eta >= 0 && fid_zeta >= 0;
}

bool IOModule::AddTecFieldOrGroup_(const std::string &name,
                                   std::vector<TecVar> &vars,
                                   std::unordered_set<std::string> &seen_scalars,
                                   std::unordered_set<std::string> &seen_forms) const
{
    if (name.empty())
        return true;

    int form_xi = -1, form_eta = -1, form_zeta = -1;
    bool form_is_face = false;
    std::string form_output_name;
    if (BuildTecFormTriplet_(name, form_xi, form_eta, form_zeta, form_is_face, form_output_name))
    {
        const std::string key = (form_is_face ? "face:" : "edge:") + form_output_name;
        if (!seen_forms.insert(key).second)
            return true;

        const bool all_node = NormalizedTecplotMode_() == TecplotMode::AllNode ||
                              tec_form_reconstruction_ == TecplotFormReconstruction::ToNode;
        const int32_t loc = all_node ? 0 : 1;

        auto custom = tec_comp_names_.find(form_output_name);
        if (custom == tec_comp_names_.end())
            custom = tec_comp_names_.find(name);
        const std::vector<std::string> default_names = {
            form_output_name + "_x",
            form_output_name + "_y",
            form_output_name + "_z"};

        for (int c = 0; c < 3; ++c)
        {
            TecVar tv;
            tv.kind = TecVar::Kind::ReconstructedForm;
            tv.comp = c;
            tv.loc = loc;
            tv.form_fid_xi = form_xi;
            tv.form_fid_eta = form_eta;
            tv.form_fid_zeta = form_zeta;
            tv.form_is_face = form_is_face;
            if (custom != tec_comp_names_.end() && c < static_cast<int>(custom->second.size()))
                tv.name = custom->second[c];
            else
                tv.name = default_names[c];
            vars.push_back(std::move(tv));
        }
        return true;
    }

    if (!fld_->has_field(name))
    {
        bool expanded_group = false;
        for (const FieldDescriptor &gd : fld_->descriptors())
        {
            if (gd.sync.group != name)
                continue;
            if (gd.location != StaggerLocation::Cell && gd.location != StaggerLocation::Node)
                continue;

            expanded_group = true;
            AddTecFieldOrGroup_(gd.name, vars, seen_scalars, seen_forms);
        }

        if (expanded_group)
            return true;

        std::fprintf(stderr, "[IOModule][Tecplot] skip missing field or group: %s\n", name.c_str());
        return false;
    }

    int fid = fld_->field_id(name);
    const auto &d = fld_->descriptor(fid);

    if (d.location != StaggerLocation::Cell && d.location != StaggerLocation::Node)
    {
        std::fprintf(stderr, "[IOModule][Tecplot] skip unsupported standalone field=%s loc=%d kind=%s\n",
                     name.c_str(), static_cast<int>(d.location), field_value_kind_name(d.value_kind));
        return false;
    }

    if (!seen_scalars.insert(d.name).second)
        return true;

    auto it = tec_comp_names_.find(name);
    const bool has_custom_names = (it != tec_comp_names_.end());
    const TecplotMode mode = NormalizedTecplotMode_();

    for (int c = 0; c < d.ncomp; ++c)
    {
        TecVar tv;
        tv.kind = TecVar::Kind::Field;
        tv.fid = fid;
        tv.comp = c;

        if (mode == TecplotMode::CellAndNode)
            tv.loc = (d.location == StaggerLocation::Cell) ? 1 : 0;
        else
            tv.loc = 0;

        if (has_custom_names && c < static_cast<int>(it->second.size()))
            tv.name = it->second[c];
        else
            tv.name = name + "_" + std::to_string(c);

        vars.push_back(std::move(tv));
    }

    return true;
}

std::vector<IOModule::TecVar> IOModule::BuildTecVars_() const
{
    std::vector<TecVar> vars;

    // 先放坐标（名字固定）
    for (const std::string &coord : {"x", "y", "z"})
    {
        TecVar tv;
        tv.kind = TecVar::Kind::Coordinate;
        tv.comp = static_cast<int>(vars.size());
        tv.loc = 0;
        tv.name = coord;
        vars.push_back(std::move(tv));
    }

    std::unordered_set<std::string> seen_scalars;
    std::unordered_set<std::string> seen_forms;
    for (const auto &fname : tec_fields_)
        AddTecFieldOrGroup_(fname, vars, seen_scalars, seen_forms);

    return vars;
}

void IOModule::EvalCoord_Node_(int ib, int i, int j, int k, double &x, double &y, double &z) const
{
    Block &blk = grd_->grids(ib);
    x = blk.x(i, j, k);
    y = blk.y(i, j, k);
    z = blk.z(i, j, k);
}

void IOModule::EvalCoord_Cell_(int ib, int i, int j, int k, double &x, double &y, double &z) const
{
    // cell center: dual_x(i+1,j+1,k+1) matches your grid convention
    Block &blk = grd_->grids(ib);
    x = blk.dual_x(i + 1, j + 1, k + 1);
    y = blk.dual_y(i + 1, j + 1, k + 1);
    z = blk.dual_z(i + 1, j + 1, k + 1);
}

double IOModule::EvalValue_CellAsNode_(const TecVar &tv, int ib, int i, int j, int k) const
{
    if (tv.kind == TecVar::Kind::Coordinate)
    {
        double x, y, z;
        EvalCoord_Cell_(ib, i, j, k, x, y, z);
        if (tv.name == "x")
            return x;
        if (tv.name == "y")
            return y;
        return z;
    }

    if (tv.kind == TecVar::Kind::ReconstructedForm)
        return EvalValue_ReconstructedFormAtCell_(tv, ib, i, j, k);

    const auto &d = fld_->descriptor(tv.fid);
    FieldBlock &fb = fld_->field(tv.fid, ib);
    if (!fb.is_allocated())
        return 0.0;

    // CellAsNode 的网格点就是 cell center
    if (d.location == StaggerLocation::Cell)
    {
        return fb(i, j, k, tv.comp);
    }

    // Node -> Cell 平均（取 2^dim 个节点，边界自动裁剪）
    // cell(i,j,k) uses nodes (i,i+1)*(j,j+1)*(k,k+1)
    const Block &blk = grd_->grids(ib);
    const int dim = blk.dimension;

    const int NiN = blk.mx + 1;
    const int NjN = blk.my + 1;
    const int NkN = (dim == 2) ? 1 : (blk.mz + 1);

    double sum = 0.0;
    int cnt = 0;

    auto add_node = [&](int ii, int jj, int kk)
    {
        if (ii < 0 || ii >= NiN)
            return;
        if (jj < 0 || jj >= NjN)
            return;
        if (kk < 0 || kk >= NkN)
            return;
        sum += fb(ii, jj, kk, tv.comp);
        cnt += 1;
    };

    add_node(i, j, k);
    add_node(i + 1, j, k);
    add_node(i, j + 1, k);
    add_node(i + 1, j + 1, k);

    if (dim == 3)
    {
        add_node(i, j, k + 1);
        add_node(i + 1, j, k + 1);
        add_node(i, j + 1, k + 1);
        add_node(i + 1, j + 1, k + 1);
    }

    return (cnt > 0) ? (sum / cnt) : 0.0;
}

double IOModule::EvalValue_CellToNode_(const TecVar &tv, int ib, int i, int j, int k) const
{
    if (tv.kind == TecVar::Kind::Coordinate)
    {
        double x, y, z;
        EvalCoord_Node_(ib, i, j, k, x, y, z);
        if (tv.name == "x")
            return x;
        if (tv.name == "y")
            return y;
        return z;
    }

    if (tv.kind == TecVar::Kind::ReconstructedForm)
        return EvalValue_ReconstructedFormAtNode_(tv, ib, i, j, k);

    const auto &d = fld_->descriptor(tv.fid);
    FieldBlock &fb = fld_->field(tv.fid, ib);
    if (!fb.is_allocated())
        return 0.0;

    // Node field: direct
    if (d.location == StaggerLocation::Node)
        return fb(i, j, k, tv.comp);

    // Singular nodes must be reconstructed from the quotient-mesh incident
    // cells.  Their block-local ghost cells are degenerate aliases and are not
    // trustworthy interpolation samples.
    double singular_value = 0.0;
    if (EvalTecplotSingularNodeAverage_(tv, ib, i, j, k, singular_value))
        return singular_value;

    // Cell -> Node volume-weighted average over the complete 2^dim cell
    // stencil around the node.  |Jac| is the physical cell-volume weight.
    // node.  At a block interface the cells on the other side are available
    // through the synchronized ghost layers.  Using them here is essential:
    // clipping to the block-interior range turns an interface-node value into
    // two unrelated one-sided averages, one in each Tecplot zone, and creates
    // an artificial jump even when the cell data are continuous.
    const Block &blk = grd_->grids(ib);
    const int dim = blk.dimension;

    const Int3 cell_lo = fb.get_lo();
    const Int3 cell_hi = fb.get_hi();
    FieldBlock &jac = fld_->field("Jac", ib);

    double weighted_sum = 0.0;
    double weight_sum = 0.0;

    auto add_cell = [&](int ii, int jj, int kk)
    {
        if (ii < cell_lo.i || ii >= cell_hi.i)
            return;
        if (jj < cell_lo.j || jj >= cell_hi.j)
            return;
        if (kk < cell_lo.k || kk >= cell_hi.k)
            return;
        if (!jac.is_allocated())
            return;
        const double weight = std::abs(jac(ii, jj, kk, 0));
        if (!std::isfinite(weight) || weight <= 0.0)
            return;
        weighted_sum += weight * fb(ii, jj, kk, tv.comp);
        weight_sum += weight;
    };

    add_cell(i - 1, j - 1, k - 1);
    add_cell(i, j - 1, k - 1);
    add_cell(i - 1, j, k - 1);
    add_cell(i, j, k - 1);

    if (dim == 3)
    {
        add_cell(i - 1, j - 1, k);
        add_cell(i, j - 1, k);
        add_cell(i - 1, j, k);
        add_cell(i, j, k);
    }
    else
    {
        // 2D: k 固定 0
        weighted_sum = 0.0;
        weight_sum = 0.0;
        add_cell(i - 1, j - 1, 0);
        add_cell(i, j - 1, 0);
        add_cell(i - 1, j, 0);
        add_cell(i, j, 0);
    }

    return (weight_sum > 0.0) ? (weighted_sum / weight_sum) : 0.0;
}

std::uint64_t IOModule::TecplotFieldComponentKey_(int fid, int comp)
{
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(fid)) << 32) |
           static_cast<std::uint32_t>(comp);
}

void IOModule::BuildTecplotSingularNodeStencils_()
{
    tec_singular_node_stencils_.clear();
    tec_singular_node_alias_to_stencil_.clear();
    tec_singular_node_averages_.clear();

    if (!tec_topology_ || !tec_singular_edges_ || tec_singular_edges_->empty())
    {
        tec_singular_stencil_ready_ = true;
        return;
    }

    const auto &topology = *tec_topology_;

    // Node classes are globally gathered.  Build a lookup for both local and
    // remote aliases so endpoints of every singular-edge alias resolve to the
    // same quotient-node id.
    std::unordered_map<TOPO::EntityKey, int, TOPO::EntityKey::Hash> node_gid;
    for (const auto &cls : topology.nodes.classes)
    {
        if (cls.global_id < 0)
            continue;
        node_gid[cls.owner.entity] = cls.global_id;
        for (const auto &member : cls.members)
            node_gid[member.entity] = cls.global_id;
    }
    for (const auto &[node, rep] : topology.nodes.local_to_rep)
    {
        const auto gid_it = topology.nodes.rep_to_qid.find(rep);
        if (gid_it != topology.nodes.rep_to_qid.end())
            node_gid[node] = gid_it->second;
    }

    std::set<int> singular_gids;
    std::unordered_map<int,
                       std::set<TOPO::EntityKey>,
                       std::hash<int>> local_cells_by_gid;

    for (const auto &edge : tec_singular_edges_->entries())
    {
        std::set<int> endpoint_gids;
        for (const auto &alias : edge.aliases)
        {
            const auto endpoints = TOPO::endpoints(alias.edge);
            for (const TOPO::EntityKey *endpoint : {&endpoints.first, &endpoints.second})
            {
                const auto gid_it = node_gid.find(*endpoint);
                if (gid_it != node_gid.end())
                    endpoint_gids.insert(gid_it->second);
            }
        }

        if (endpoint_gids.size() != 2)
        {
            throw std::runtime_error(
                "[IOModule][Tecplot] singular edge does not resolve to exactly two quotient nodes: gid=" +
                std::to_string(edge.global_id));
        }

        for (const int gid : endpoint_gids)
        {
            singular_gids.insert(gid);
            auto &cells = local_cells_by_gid[gid];
            for (const auto &incident : edge.local_incident_cells)
                cells.insert(incident.entity);
        }
    }

    tec_singular_node_stencils_.reserve(singular_gids.size());
    std::unordered_map<int, std::size_t> gid_to_stencil;
    for (const int gid : singular_gids)
    {
        TecSingularNodeStencil stencil;
        stencil.node_gid = gid;
        const auto cells_it = local_cells_by_gid.find(gid);
        if (cells_it != local_cells_by_gid.end())
            stencil.local_cells.assign(cells_it->second.begin(), cells_it->second.end());
        gid_to_stencil[gid] = tec_singular_node_stencils_.size();
        tec_singular_node_stencils_.push_back(std::move(stencil));
    }

    // Every local Tecplot node alias points at the one physical-node stencil.
    for (const auto &[node, rep] : topology.nodes.local_to_rep)
    {
        const auto gid_it = topology.nodes.rep_to_qid.find(rep);
        if (gid_it == topology.nodes.rep_to_qid.end())
            continue;
        const auto stencil_it = gid_to_stencil.find(gid_it->second);
        if (stencil_it != gid_to_stencil.end())
            tec_singular_node_alias_to_stencil_[node] = stencil_it->second;
    }

    tec_singular_stencil_ready_ = true;
}

void IOModule::PrepareTecplotSingularNodeAverages_(const std::vector<TecVar> &vars)
{
    if (!tec_singular_stencil_ready_)
        BuildTecplotSingularNodeStencils_();

    tec_singular_node_averages_.clear();
    const std::size_t nnode = tec_singular_node_stencils_.size();
    if (nnode == 0)
        return;

    std::vector<std::uint64_t> keys;
    std::unordered_set<std::uint64_t> seen;
    for (const TecVar &tv : vars)
    {
        if (tv.kind != TecVar::Kind::Field || tv.fid < 0)
            continue;
        const auto &desc = fld_->descriptor(tv.fid);
        if (desc.location != StaggerLocation::Cell)
            continue;
        const std::uint64_t key = TecplotFieldComponentKey_(tv.fid, tv.comp);
        if (seen.insert(key).second)
            keys.push_back(key);
    }

    const std::size_t nkey = keys.size();
    if (nkey == 0)
        return;

    std::vector<double> local_sum(nkey * nnode, 0.0);
    std::vector<double> global_sum(nkey * nnode, 0.0);
    std::vector<double> local_weight(nkey * nnode, 0.0);
    std::vector<double> global_weight(nkey * nnode, 0.0);

    for (std::size_t q = 0; q < nkey; ++q)
    {
        const int fid = static_cast<int>(keys[q] >> 32);
        const int comp = static_cast<int>(keys[q] & 0xffffffffu);
        for (std::size_t n = 0; n < nnode; ++n)
        {
            const std::size_t out = q * nnode + n;
            for (const TOPO::EntityKey &cell : tec_singular_node_stencils_[n].local_cells)
            {
                if (cell.rank != myid_ || cell.block < 0 || cell.block >= fld_->num_blocks())
                    continue;
                FieldBlock &fb = fld_->field(fid, cell.block);
                if (!fb.is_allocated())
                    continue;
                FieldBlock &jac = fld_->field("Jac", cell.block);
                if (!jac.is_allocated())
                    continue;
                const Int3 lo = fb.inner_lo();
                const Int3 hi = fb.inner_hi();
                if (cell.i < lo.i || cell.i >= hi.i ||
                    cell.j < lo.j || cell.j >= hi.j ||
                    cell.k < lo.k || cell.k >= hi.k)
                    continue;
                const double weight = std::abs(jac(cell.i, cell.j, cell.k, 0));
                if (!std::isfinite(weight) || weight <= 0.0)
                    continue;
                local_sum[out] += weight * fb(cell.i, cell.j, cell.k, comp);
                local_weight[out] += weight;
            }
        }
    }

    MPI_Allreduce(local_sum.data(), global_sum.data(),
                  static_cast<int>(local_sum.size()), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(local_weight.data(), global_weight.data(),
                  static_cast<int>(local_weight.size()), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for (std::size_t q = 0; q < nkey; ++q)
    {
        auto &average = tec_singular_node_averages_[keys[q]];
        average.resize(nnode, 0.0);
        for (std::size_t n = 0; n < nnode; ++n)
        {
            const std::size_t in = q * nnode + n;
            if (global_weight[in] > 0.0)
                average[n] = global_sum[in] / global_weight[in];
        }
    }
}

bool IOModule::EvalTecplotSingularNodeAverage_(const TecVar &tv,
                                               int ib, int i, int j, int k,
                                               double &value) const
{
    if (!tec_topology_ || tv.kind != TecVar::Kind::Field)
        return false;

    const TOPO::EntityKey node = TOPO::make_node(myid_, ib, i, j, k);
    const auto stencil_it = tec_singular_node_alias_to_stencil_.find(node);
    if (stencil_it == tec_singular_node_alias_to_stencil_.end())
        return false;

    const auto value_it = tec_singular_node_averages_.find(
        TecplotFieldComponentKey_(tv.fid, tv.comp));
    if (value_it == tec_singular_node_averages_.end() ||
        stencil_it->second >= value_it->second.size())
        return false;

    value = value_it->second[stencil_it->second];
    return true;
}

double IOModule::EvalValue_Mixed_(const TecVar &tv, int ib, int i, int j, int k) const
{
    // Mixed: coords are node coords, variables follow their own location
    if (tv.kind == TecVar::Kind::Coordinate)
    {
        double x, y, z;
        EvalCoord_Node_(ib, i, j, k, x, y, z);
        if (tv.name == "x")
            return x;
        if (tv.name == "y")
            return y;
        return z;
    }

    if (tv.kind == TecVar::Kind::ReconstructedForm)
    {
        if (tv.loc == 0)
            return EvalValue_ReconstructedFormAtNode_(tv, ib, i, j, k);
        return EvalValue_ReconstructedFormAtCell_(tv, ib, i, j, k);
    }

    FieldBlock &fb = fld_->field(tv.fid, ib);
    if (!fb.is_allocated())
        return 0.0;

    return fb(i, j, k, tv.comp);
}

double IOModule::EvalValue_ReconstructedFormAtCell_(const TecVar &tv, int ib, int i, int j, int k) const
{
    double xyz[3] = {0.0, 0.0, 0.0};
    bool ok = false;
    if (tv.form_is_face)
        ok = METRIC::reconstruct_face_2form_to_cell(*fld_,
                                                    tv.form_fid_xi, tv.form_fid_eta, tv.form_fid_zeta,
                                                    ib, i, j, k, xyz);
    else
        ok = METRIC::reconstruct_edge_1form_to_cell(*fld_,
                                                    tv.form_fid_xi, tv.form_fid_eta, tv.form_fid_zeta,
                                                    ib, i, j, k, xyz);
    return ok ? xyz[tv.comp] : 0.0;
}

double IOModule::EvalValue_ReconstructedFormAtNode_(const TecVar &tv, int ib, int i, int j, int k) const
{
    double xyz[3] = {0.0, 0.0, 0.0};
    bool ok = false;
    if (tv.form_is_face)
        ok = METRIC::reconstruct_face_2form_to_node(*fld_,
                                                    tv.form_fid_xi, tv.form_fid_eta, tv.form_fid_zeta,
                                                    ib, i, j, k, xyz);
    else
        ok = METRIC::reconstruct_edge_1form_to_node(*fld_,
                                                    tv.form_fid_xi, tv.form_fid_eta, tv.form_fid_zeta,
                                                    ib, i, j, k, xyz);
    return ok ? xyz[tv.comp] : 0.0;
}

void IOModule::plt_write_header_(FILE *fp, const std::string &title,
                                 const std::vector<std::string> &var_names)
{
    const char magic[8] = {'#', '!', 'T', 'D', 'V', '1', '1', '2'};
    std::fwrite(magic, 8, 1, fp);

    TecWriteI32_(fp, 1); // little endian
    TecWriteI32_(fp, 0); // FULL

    // 你原先用 plt_write_str_，这里直接用 TecWriteStr_ 即可
    TecWriteStr_(fp, title);

    TecWriteI32_(fp, (int32_t)var_names.size());
    for (size_t i = 0; i < var_names.size(); ++i)
        TecWriteStr_(fp, var_names[i]);
}

void IOModule::plt_write_zone_header_(FILE *fp, const std::string &zone_name,
                                      int IMax, int JMax, int KMax,
                                      double sol_time,
                                      const std::vector<int32_t> &var_locs)
{
    TecWriteF32_(fp, 299.0f);
    TecWriteStr_(fp, zone_name);

    TecWriteI32_(fp, -1); // parent
    TecWriteI32_(fp, -2); // strand
    TecWriteF64_(fp, sol_time);
    TecWriteI32_(fp, -1); // zone color
    TecWriteI32_(fp, 0);  // ORDERED

    TecWriteI32_(fp, 1); // specify var location = true
    for (size_t i = 0; i < var_locs.size(); ++i)
        TecWriteI32_(fp, var_locs[i]); // 0 nodal / 1 cell

    TecWriteI32_(fp, 0); // raw face neighbors
    TecWriteI32_(fp, 0); // misc neighbors

    TecWriteI32_(fp, IMax);
    TecWriteI32_(fp, JMax);
    TecWriteI32_(fp, KMax);

    TecWriteI32_(fp, 0); // AuxData count
}

void IOModule::plt_write_data_section_header_(FILE *fp, int nvar)
{
    TecWriteF32_(fp, 299.0f);

    // all variables type = float
    for (int v = 0; v < nvar; ++v)
        TecWriteI32_(fp, 1);

    TecWriteI32_(fp, 0);  // passive vars
    TecWriteI32_(fp, 0);  // shared vars
    TecWriteI32_(fp, -1); // shared connectivity
}

void IOModule::WriteTecplotBinFile(int step, double time)
{
    if (!par_ || !grd_ || !fld_)
        Fail_("[IOModule][Tecplot] Setup() must be called before WriteTecplotFile()");

    const std::string path = tecplot_path_;
    const TecplotMode output_mode = NormalizedTecplotMode_();

    FILE *fp = std::fopen(path.c_str(), "wb");
    if (!fp)
        Fail_("[IOModule][Tecplot] cannot open: " + path);

    // -------- build variable list --------
    const std::vector<TecVar> vars = BuildTecVars_();
    if (output_mode == TecplotMode::AllNode)
        PrepareTecplotSingularNodeAverages_(vars);

    std::vector<std::string> var_names;
    std::vector<int32_t> var_locs;

    var_names.reserve(vars.size());
    var_locs.reserve(vars.size());

    for (size_t v = 0; v < vars.size(); ++v)
    {
        var_names.push_back(vars[v].name);

        // coords always nodal
        if (v < 3)
            var_locs.push_back(0);
        else
            var_locs.push_back(vars[v].loc);
    }

    // -------- Header section --------
    plt_write_header_(fp, "flow_field", var_names);

    // -------- Zone headers (one per selected block) --------
    const int nblock = grd_->nblock;
    for (int ib = 0; ib < nblock; ++ib)
    {
        Block &blk = grd_->grids(ib);
        if (!TecplotBlockSelected_(blk.block_name))
            continue;

        const int dim = blk.dimension;

        int IMax = 0, JMax = 0, KMax = 0;

        if (output_mode == TecplotMode::CellAsNode)
        {
            IMax = blk.mx;
            JMax = blk.my;
            KMax = (dim == 2) ? 1 : blk.mz;
        }
        else
        {
            IMax = blk.mx + 1;
            JMax = blk.my + 1;
            KMax = (dim == 2) ? 1 : (blk.mz + 1);
        }

        std::string zone_name = "Zone" + std::to_string(ib) + "_" + blk.block_name;
        plt_write_zone_header_(fp, zone_name, IMax, JMax, KMax, time, var_locs);
    }

    // end of header marker
    TecWriteF32_(fp, 357.0f);

    // -------- Data section per zone --------
    for (int ib = 0; ib < nblock; ++ib)
    {
        Block &blk = grd_->grids(ib);
        if (!TecplotBlockSelected_(blk.block_name))
            continue;

        const int dim = blk.dimension;

        const int NiC = blk.mx;
        const int NjC = blk.my;
        const int NkC = (dim == 2) ? 1 : blk.mz;

        const int NiN = blk.mx + 1;
        const int NjN = blk.my + 1;
        const int NkN = (dim == 2) ? 1 : (blk.mz + 1);
        // Tecplot binary pads ordered cell-centered data along the fastest
        // indices: write IMax*JMax*(KMax-1) values for IJK zones.
        const int NkCellWrite = (dim == 2) ? 1 : (NkN - 1);

        int IMax = 0, JMax = 0, KMax = 0;
        if (output_mode == TecplotMode::CellAsNode)
        {
            IMax = NiC;
            JMax = NjC;
            KMax = NkC;
        }
        else
        {
            IMax = NiN;
            JMax = NjN;
            KMax = NkN;
        }

        // 1 section header
        plt_write_data_section_header_(fp, (int)vars.size());

        // 2 min/max for each var
        for (size_t v = 0; v < vars.size(); ++v)
        {
            double vmin = 0.0, vmax = 0.0;
            bool first = true;

            const TecVar &tv = vars[v];

            // decide loop extents by mode + loc
            if (output_mode == TecplotMode::CellAsNode)
            {
                for (int k = 0; k < NkC; ++k)
                    for (int j = 0; j < NjC; ++j)
                        for (int i = 0; i < NiC; ++i)
                        {
                            double val = EvalValue_CellAsNode_(tv, ib, i, j, k);
                            if (first)
                            {
                                vmin = vmax = val;
                                first = false;
                            }
                            else
                            {
                                if (val < vmin)
                                    vmin = val;
                                if (val > vmax)
                                    vmax = val;
                            }
                        }
            }
            else if (output_mode == TecplotMode::AllNode)
            {
                for (int k = 0; k < NkN; ++k)
                    for (int j = 0; j < NjN; ++j)
                        for (int i = 0; i < NiN; ++i)
                        {
                            double val = EvalValue_CellToNode_(tv, ib, i, j, k);
                            if (first)
                            {
                                vmin = vmax = val;
                                first = false;
                            }
                            else
                            {
                                if (val < vmin)
                                    vmin = val;
                                if (val > vmax)
                                    vmax = val;
                            }
                        }
            }
            else // Mixed
            {
                const bool is_cell_var = (v >= 3 && tv.loc == 1);

                if (is_cell_var)
                {
                    for (int k = 0; k < NkCellWrite; ++k)
                        for (int j = 0; j < NjN; ++j)
                            for (int i = 0; i < NiN; ++i)
                            {
                                const bool is_real_cell = (i < NiC && j < NjC && k < NkC);
                                double val = is_real_cell ? EvalValue_Mixed_(tv, ib, i, j, k) : 0.0;
                                if (first)
                                {
                                    vmin = vmax = val;
                                    first = false;
                                }
                                else
                                {
                                    if (val < vmin)
                                        vmin = val;
                                    if (val > vmax)
                                        vmax = val;
                                }
                            }
                }
                else
                {
                    for (int k = 0; k < NkN; ++k)
                        for (int j = 0; j < NjN; ++j)
                            for (int i = 0; i < NiN; ++i)
                            {
                                double val = EvalValue_Mixed_(tv, ib, i, j, k);
                                if (first)
                                {
                                    vmin = vmax = val;
                                    first = false;
                                }
                                else
                                {
                                    if (val < vmin)
                                        vmin = val;
                                    if (val > vmax)
                                        vmax = val;
                                }
                            }
                }
            }

            TecWriteF64_(fp, vmin);
            TecWriteF64_(fp, vmax);
        }

        // 6.3 write data variable-major, float
        for (size_t v = 0; v < vars.size(); ++v)
        {
            const TecVar &tv = vars[v];

            if (output_mode == TecplotMode::CellAsNode)
            {
                for (int k = 0; k < NkC; ++k)
                    for (int j = 0; j < NjC; ++j)
                        for (int i = 0; i < NiC; ++i)
                        {
                            float fv = (float)EvalValue_CellAsNode_(tv, ib, i, j, k);
                            TecWriteF32_(fp, fv);
                        }
            }
            else if (output_mode == TecplotMode::AllNode)
            {
                for (int k = 0; k < NkN; ++k)
                    for (int j = 0; j < NjN; ++j)
                        for (int i = 0; i < NiN; ++i)
                        {
                            float fv = (float)EvalValue_CellToNode_(tv, ib, i, j, k);
                            TecWriteF32_(fp, fv);
                        }
            }
            else // Mixed
            {
                const bool is_cell_var = (v >= 3 && tv.loc == 1);

                if (is_cell_var)
                {
                    for (int k = 0; k < NkCellWrite; ++k)
                        for (int j = 0; j < NjN; ++j)
                            for (int i = 0; i < NiN; ++i)
                            {
                                const bool is_real_cell = (i < NiC && j < NjC && k < NkC);
                                float fv = is_real_cell ? (float)EvalValue_Mixed_(tv, ib, i, j, k) : 0.0f;
                                TecWriteF32_(fp, fv);
                            }
                }
                else
                {
                    for (int k = 0; k < NkN; ++k)
                        for (int j = 0; j < NjN; ++j)
                            for (int i = 0; i < NiN; ++i)
                            {
                                float fv = (float)EvalValue_Mixed_(tv, ib, i, j, k);
                                TecWriteF32_(fp, fv);
                            }
                }
            }
        }
    }

    std::fclose(fp);

    // if (par_->GetInt("myid") == 0)
    //     std::printf("[IOModule][Tecplot] wrote %s\n", path.c_str());
}
