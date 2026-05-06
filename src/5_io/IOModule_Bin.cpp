#include "5_io/IOModule.h"

#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <unordered_set>

void IOModule::WriteI32_(std::ofstream &out, int32_t v) { out.write(reinterpret_cast<const char *>(&v), sizeof(v)); }
void IOModule::WriteI64_(std::ofstream &out, int64_t v) { out.write(reinterpret_cast<const char *>(&v), sizeof(v)); }
void IOModule::WriteF64_(std::ofstream &out, double v) { out.write(reinterpret_cast<const char *>(&v), sizeof(v)); }

int32_t IOModule::ReadI32_(std::ifstream &in)
{
    int32_t v;
    in.read(reinterpret_cast<char *>(&v), sizeof(v));
    if (!in)
        Fail_("[IOModule] ReadI32 failed");
    return v;
}
int64_t IOModule::ReadI64_(std::ifstream &in)
{
    int64_t v;
    in.read(reinterpret_cast<char *>(&v), sizeof(v));
    if (!in)
        Fail_("[IOModule] ReadI64 failed");
    return v;
}
double IOModule::ReadF64_(std::ifstream &in)
{
    double v;
    in.read(reinterpret_cast<char *>(&v), sizeof(v));
    if (!in)
        Fail_("[IOModule] ReadF64 failed");
    return v;
}

void IOModule::WriteStr_(std::ofstream &out, const std::string &s)
{
    // length(int32) + bytes(no '\0')
    int32_t n = static_cast<int32_t>(s.size());
    WriteI32_(out, n);
    if (n > 0)
        out.write(s.data(), n);
}

std::string IOModule::ReadStr_(std::ifstream &in)
{
    int32_t n = ReadI32_(in);
    if (n < 0)
        Fail_("[IOModule] ReadStr got negative length");
    std::string s;
    s.resize(static_cast<size_t>(n));
    if (n > 0)
        in.read(&s[0], n);
    if (!in)
        Fail_("[IOModule] ReadStr payload failed");
    return s;
}

int32_t IOModule::LocToI32_(StaggerLocation loc)
{
    switch (loc)
    {
    case StaggerLocation::Cell:
        return 0;
    case StaggerLocation::FaceXi:
        return 1;
    case StaggerLocation::FaceEt:
        return 2;
    case StaggerLocation::FaceZe:
        return 3;
    case StaggerLocation::EdgeXi:
        return 4;
    case StaggerLocation::EdgeEt:
        return 5;
    case StaggerLocation::EdgeZe:
        return 6;
    case StaggerLocation::Node:
        return 7;
    default:
        return -999;
    }
}

StaggerLocation IOModule::I32ToLoc_(int32_t v)
{
    switch (v)
    {
    case 0:
        return StaggerLocation::Cell;
    case 1:
        return StaggerLocation::FaceXi;
    case 2:
        return StaggerLocation::FaceEt;
    case 3:
        return StaggerLocation::FaceZe;
    case 4:
        return StaggerLocation::EdgeXi;
    case 5:
        return StaggerLocation::EdgeEt;
    case 6:
        return StaggerLocation::EdgeZe;
    case 7:
        return StaggerLocation::Node;
    default:
        return StaggerLocation::Cell;
    }
}

//=========================================================================

void IOModule::SetRestartFields(const std::vector<std::string> &names)
{
    restart_field_names_.clear();
    restart_field_names_.reserve(names.size());

    std::unordered_set<std::string> seen;
    for (const auto &nm : names)
    {
        if (nm.empty())
            continue;
        if (seen.insert(nm).second)
            restart_field_names_.push_back(nm);
    }
}

void IOModule::AddRestartField(const std::string &name)
{
    if (name.empty())
        return;

    for (auto &s : restart_field_names_)
        if (s == name)
            return;

    restart_field_names_.push_back(name);
}

void IOModule::ClearRestartFields()
{
    restart_field_names_.clear();
}

std::vector<int> IOModule::BuildRestartFids_() const
{
    if (restart_field_names_.empty())
    {
        Fail_("[IOModule] Restart field whitelist is empty. "
              "Call SetRestartFields(...) before WriteRestart().");
    }

    std::vector<int> fids;
    fids.reserve(restart_field_names_.size());

    for (const auto &nm : restart_field_names_)
    {
        if (!fld_->has_field(nm))
        {
            Fail_("[IOModule] Restart whitelist contains unknown field: " + nm);
        }
        fids.push_back(fld_->field_id(nm));
    }
    return fids;
}

void IOModule::WriteRestartBinFile(int step, double time)
{
    std::string path = restart_path_;
    std::ofstream out(path, std::ios::binary);
    if (!out)
        Fail_("[IOModule] WriteRestart: cannot open " + path);

    // Header
    const char magic[8] = {'M', 'P', 'C', 'N', 'S', 'R', 'S', 'T'};
    out.write(magic, 8);
    WriteI32_(out, 1); // version
    WriteI32_(out, static_cast<int32_t>(step));
    WriteF64_(out, time);

    const int32_t nblock = static_cast<int32_t>(fld_->num_blocks());

    // 关键：只输出白名单字段
    const std::vector<int> fids = BuildRestartFids_();
    const int32_t nfield = static_cast<int32_t>(fids.size());

    WriteI32_(out, nblock);
    WriteI32_(out, nfield);

    for (int fid : fids)
    {
        const auto &d = fld_->descriptor(fid);

        WriteStr_(out, d.name);
        WriteI32_(out, LocToI32_(d.location));
        WriteI32_(out, static_cast<int32_t>(d.ncomp));
        WriteI32_(out, static_cast<int32_t>(d.nghost));

        for (int ib = 0; ib < nblock; ++ib)
        {
            FieldBlock &fb = fld_->field(fid, ib);

            const Int3 lo = fb.get_lo();
            const Int3 hi = fb.get_hi();
            WriteI32_(out, lo.i);
            WriteI32_(out, lo.j);
            WriteI32_(out, lo.k);
            WriteI32_(out, hi.i);
            WriteI32_(out, hi.j);
            WriteI32_(out, hi.k);

            const bool active = fb.is_allocated();
            WriteI32_(out, active ? 1 : 0); // active flag
            if (!active)
                continue;

            const int ncomp_i = d.ncomp;

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                        for (int m = 0; m < ncomp_i; ++m)
                        {
                            double v = fb(i, j, k, m);
                            WriteF64_(out, v);
                        }
        }
    }

    out.close();
}

void IOModule::ReadRestartBinFile()
{
    std::string path = restart_path_;
    std::ifstream in(path, std::ios::binary);
    if (!in)
        Fail_("[IOModule] ReadRestart: cannot open " + path);

    // ---- Header ----
    char magic[8];
    in.read(magic, 8);
    if (!in)
        Fail_("[IOModule] ReadRestart: bad header (magic read failed)");

    const char expect[8] = {'M', 'P', 'C', 'N', 'S', 'R', 'S', 'T'};
    if (std::memcmp(magic, expect, 8) != 0)
        Fail_("[IOModule] ReadRestart: magic mismatch");

    int32_t version = ReadI32_(in);
    if (version != 1)
        Fail_("[IOModule] ReadRestart: unsupported version=" + std::to_string(version));

    int32_t step = ReadI32_(in);
    double time = ReadF64_(in);

    int32_t nblock_file = ReadI32_(in);
    int32_t nfield_file = ReadI32_(in);

    const int32_t nblock = static_cast<int32_t>(fld_->num_blocks());
    if (nblock_file != nblock)
        Fail_("[IOModule] ReadRestart: nblock mismatch file=" + std::to_string(nblock_file) +
              " current=" + std::to_string(nblock));

    // 可选：写回参数（便于续算控制）
    par_->AddParam("Nstep", (int)step);
    par_->AddParam("Physic_Time", time);

    // ---- Fields ----
    for (int32_t ifld = 0; ifld < nfield_file; ++ifld)
    {
        std::string fname = ReadStr_(in);
        int32_t loc_i32 = ReadI32_(in);
        int32_t ncomp = ReadI32_(in);
        int32_t nghost = ReadI32_(in);

        const bool exists = fld_->has_field(fname);
        int fid = exists ? fld_->field_id(fname) : -1;

        // 若存在，做 descriptor 强校验
        if (exists)
        {
            const auto &d = fld_->descriptor(fid);

            if (LocToI32_(d.location) != loc_i32)
                Fail_("[IOModule] ReadRestart: location mismatch field=" + fname);

            if (d.ncomp != ncomp)
                Fail_("[IOModule] ReadRestart: ncomp mismatch field=" + fname +
                      " file=" + std::to_string(ncomp) +
                      " current=" + std::to_string(d.ncomp));

            if (d.nghost != nghost)
                Fail_("[IOModule] ReadRestart: nghost mismatch field=" + fname +
                      " file=" + std::to_string(nghost) +
                      " current=" + std::to_string(d.nghost));
        }

        for (int32_t ib = 0; ib < nblock_file; ++ib)
        {
            int32_t lo_i = ReadI32_(in);
            int32_t lo_j = ReadI32_(in);
            int32_t lo_k = ReadI32_(in);
            int32_t hi_i = ReadI32_(in);
            int32_t hi_j = ReadI32_(in);
            int32_t hi_k = ReadI32_(in);

            const int Ni = hi_i - lo_i;
            const int Nj = hi_j - lo_j;
            const int Nk = hi_k - lo_k;

            if (Ni < 0 || Nj < 0 || Nk < 0)
                Fail_("[IOModule] ReadRestart: negative extent field=" + fname +
                      " block=" + std::to_string((int)ib));

            // 新格式：每个 block 都有 active flag（但 version 仍为 1）
            int32_t file_active = ReadI32_(in); // 0/1

            const std::size_t nval =
                static_cast<std::size_t>(Ni) *
                static_cast<std::size_t>(Nj) *
                static_cast<std::size_t>(Nk) *
                static_cast<std::size_t>(ncomp);

            // --- field 不存在：跳过（若 file_active==1 才有数据段） ---
            if (!exists)
            {
                if (file_active == 1)
                {
                    in.seekg(static_cast<std::streamoff>(nval * sizeof(double)), std::ios::cur);
                    if (!in)
                        Fail_("[IOModule] ReadRestart: seek failed skipping field=" + fname);
                }
                continue;
            }

            FieldBlock &fb = fld_->field(fid, ib);

            // lo/hi 必须一致（您要求）
            const Int3 lo = fb.get_lo();
            const Int3 hi = fb.get_hi();
            if (lo.i != lo_i || lo.j != lo_j || lo.k != lo_k ||
                hi.i != hi_i || hi.j != hi_j || hi.k != hi_k)
            {
                Fail_("[IOModule] ReadRestart: lo/hi mismatch field=" + fname +
                      " block=" + std::to_string((int)ib));
            }

            const bool mem_active = fb.is_allocated();

            // 文件 inactive：不应该有数据段
            if (file_active == 0)
            {
                // 严格模式：文件没有数据但内存 block 是 active -> 无法恢复，直接报错
                if (mem_active)
                {
                    Fail_("[IOModule] ReadRestart: file block inactive but memory block active. field=" +
                          fname + " block=" + std::to_string((int)ib));
                }
                continue;
            }

            // 文件 active==1：有数据段
            if (!mem_active)
            {
                // 内存未分配：跳过数据段
                in.seekg(static_cast<std::streamoff>(nval * sizeof(double)), std::ios::cur);
                if (!in)
                    Fail_("[IOModule] ReadRestart: seek failed discarding data field=" + fname);
                continue;
            }

            // mem_active && file_active==1：逐点读入
            for (int i = lo_i; i < hi_i; ++i)
                for (int j = lo_j; j < hi_j; ++j)
                    for (int k = lo_k; k < hi_k; ++k)
                        for (int m = 0; m < ncomp; ++m)
                        {
                            double v = ReadF64_(in);
                            fb(i, j, k, m) = v;
                        }
        }
    }

    in.close();
}
