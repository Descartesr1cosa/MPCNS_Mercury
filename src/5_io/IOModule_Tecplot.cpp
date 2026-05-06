#include "5_io/IOModule.h"

#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
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

std::vector<IOModule::TecVar> IOModule::BuildTecVars_() const
{
    std::vector<TecVar> vars;

    // 先放坐标（名字固定）
    vars.push_back(TecVar{-1, 0, 0, "x"});
    vars.push_back(TecVar{-1, 0, 0, "y"});
    vars.push_back(TecVar{-1, 0, 0, "z"});

    // 用户白名单字段
    for (const auto &fname : tec_fields_)
    {
        if (!fld_->has_field(fname))
        {
            std::fprintf(stderr, "[IOModule][Tecplot] skip missing field: %s\n", fname.c_str());
            continue;
        }

        int fid = fld_->field_id(fname);
        const auto &d = fld_->descriptor(fid);

        // 目前只支持 Cell/Node
        if (d.location != StaggerLocation::Cell && d.location != StaggerLocation::Node)
        {
            std::fprintf(stderr, "[IOModule][Tecplot] skip unsupported location field=%s loc=%d\n",
                         fname.c_str(), static_cast<int>(d.location));
            continue;
        }

        // 组件名映射（可选）
        auto it = tec_comp_names_.find(fname);
        const bool has_custom_names = (it != tec_comp_names_.end());

        for (int c = 0; c < d.ncomp; ++c)
        {
            TecVar tv;
            tv.fid = fid;
            tv.comp = c;

            // loc 取决于模式
            if (tec_mode_ == TecplotMode::Mixed)
            {
                tv.loc = (d.location == StaggerLocation::Cell) ? 1 : 0;
            }
            else
            {
                tv.loc = 0; // CellAsNode / CellToNode 都是全 nodal
            }

            if (has_custom_names && c < static_cast<int>(it->second.size()))
                tv.name = it->second[c];
            else
                tv.name = fname + "_" + std::to_string(c);

            vars.push_back(std::move(tv));
        }
    }

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
    if (tv.fid < 0)
    {
        double x, y, z;
        EvalCoord_Cell_(ib, i, j, k, x, y, z);
        if (tv.name == "x")
            return x;
        if (tv.name == "y")
            return y;
        return z;
    }

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
        // if (ii < 0 || ii >= NiN)
        //     return;
        // if (jj < 0 || jj >= NjN)
        //     return;
        // if (kk < 0 || kk >= NkN)
        //     return;
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
    if (tv.fid < 0)
    {
        double x, y, z;
        EvalCoord_Node_(ib, i, j, k, x, y, z);
        if (tv.name == "x")
            return x;
        if (tv.name == "y")
            return y;
        return z;
    }

    const auto &d = fld_->descriptor(tv.fid);
    FieldBlock &fb = fld_->field(tv.fid, ib);
    if (!fb.is_allocated())
        return 0.0;

    // Node field: direct
    if (d.location == StaggerLocation::Node)
        return fb(i, j, k, tv.comp);

    // Cell -> Node 平均（相邻 cell: (i-1,i)*(j-1,j)*(k-1,k)，裁剪到有效 cell 范围）
    const Block &blk = grd_->grids(ib);
    const int dim = blk.dimension;

    const int NiC = blk.mx;
    const int NjC = blk.my;
    const int NkC = (dim == 2) ? 1 : blk.mz;

    double sum = 0.0;
    int cnt = 0;

    auto add_cell = [&](int ii, int jj, int kk)
    {
        // if (ii < 0 || ii >= NiC)
        //     return;
        // if (jj < 0 || jj >= NjC)
        //     return;
        // if (kk < 0 || kk >= NkC)
        //     return;
        sum += fb(ii, jj, kk, tv.comp);
        cnt += 1;
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
        sum = 0.0;
        cnt = 0;
        add_cell(i - 1, j - 1, 0);
        add_cell(i, j - 1, 0);
        add_cell(i - 1, j, 0);
        add_cell(i, j, 0);
    }

    return (cnt > 0) ? (sum / cnt) : 0.0;
}

double IOModule::EvalValue_Mixed_(const TecVar &tv, int ib, int i, int j, int k) const
{
    // Mixed: coords are node coords, variables follow their own location
    if (tv.fid < 0)
    {
        double x, y, z;
        EvalCoord_Node_(ib, i, j, k, x, y, z);
        if (tv.name == "x")
            return x;
        if (tv.name == "y")
            return y;
        return z;
    }

    FieldBlock &fb = fld_->field(tv.fid, ib);
    if (!fb.is_allocated())
        return 0.0;

    return fb(i, j, k, tv.comp);
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

    FILE *fp = std::fopen(path.c_str(), "wb");
    if (!fp)
        Fail_("[IOModule][Tecplot] cannot open: " + path);

    // -------- build variable list --------
    const std::vector<TecVar> vars = BuildTecVars_();

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
        const int dim = blk.dimension;

        int IMax = 0, JMax = 0, KMax = 0;

        if (tec_mode_ == TecplotMode::CellAsNode)
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
        const int dim = blk.dimension;

        const int NiC = blk.mx;
        const int NjC = blk.my;
        const int NkC = (dim == 2) ? 1 : blk.mz;

        const int NiN = blk.mx + 1;
        const int NjN = blk.my + 1;
        const int NkN = (dim == 2) ? 1 : (blk.mz + 1);

        int IMax = 0, JMax = 0, KMax = 0;
        if (tec_mode_ == TecplotMode::CellAsNode)
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
            if (tec_mode_ == TecplotMode::CellAsNode)
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
            else if (tec_mode_ == TecplotMode::CellToNode)
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
                    for (int k = 0; k < NkC; ++k)
                        for (int j = 0; j < NjC; ++j)
                            for (int i = 0; i < NiC; ++i)
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

            if (tec_mode_ == TecplotMode::CellAsNode)
            {
                for (int k = 0; k < NkC; ++k)
                    for (int j = 0; j < NjC; ++j)
                        for (int i = 0; i < NiC; ++i)
                        {
                            float fv = (float)EvalValue_CellAsNode_(tv, ib, i, j, k);
                            TecWriteF32_(fp, fv);
                        }
            }
            else if (tec_mode_ == TecplotMode::CellToNode)
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
                    for (int k = 0; k < NkC; ++k)
                        for (int j = 0; j < NjC; ++j)
                            for (int i = 0; i < NiC; ++i)
                            {
                                float fv = (float)EvalValue_Mixed_(tv, ib, i, j, k);
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
