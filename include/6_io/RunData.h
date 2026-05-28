#pragma once

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <string>
#include <vector>
class RunData
{
public:
    // ---------- core dynamic state ----------
    int32_t step = 0;
    double time = 0.0;
    double dt = 0.0;

    // residual reference: size = nvar (defined by solver)
    std::vector<double> residual_ref;

    // typical global residual indicator for stop/print
    double residual_max = 0.0;

    // ---------- variable monitors (multi-physics extendable) ----------
    // key example: "Fluid/H/rho_min", "Chem/Na/photo_max", "MHD/Bmax"
    std::map<std::string, double> monitor_d;
    std::map<std::string, int32_t> monitor_i;

    // output path
    std::string rundata_path_;

public:
    RunData()
    {
        rundata_path_ = "./DATA/RunData.bin";
    };
    ~RunData() = default;

    void ResetMonitors()
    {
        monitor_d.clear();
        monitor_i.clear();
    }

    void ResizeResidualRef(int nvar, double init_value = 1.0)
    {
        if (nvar < 0)
            nvar = 0;
        residual_ref.assign((size_t)nvar, init_value);
    }

    // -------- set / get monitors --------
    void SetMonitor(const std::string &key, double value)
    {
        // 可选：禁止同名跨类型复用（如需开启，取消注释）
        // if (monitor_i.find(key) != monitor_i.end())
        //     Fail_("RunData::SetMonitor(double) key already exists in monitor_i: " + key);

        monitor_d[key] = value;
    }

    void SetMonitor(const std::string &key, int value)
    {
        // 可选：禁止同名跨类型复用（如需开启，取消注释）
        // if (monitor_d.find(key) != monitor_d.end())
        //     Fail_("RunData::SetMonitor(int) key already exists in monitor_d: " + key);

        monitor_i[key] = static_cast<int32_t>(value);
    }

    bool HasMonitorD(const std::string &key) const
    {
        return monitor_d.find(key) != monitor_d.end();
    }

    bool HasMonitorI(const std::string &key) const
    {
        return monitor_i.find(key) != monitor_i.end();
    }

    double GetMonitorD(const std::string &key, double default_value = 0.0) const
    {
        std::map<std::string, double>::const_iterator it = monitor_d.find(key);
        if (it == monitor_d.end())
            return default_value;
        return it->second;
    }

    int32_t GetMonitorI(const std::string &key, int32_t default_value = 0) const
    {
        std::map<std::string, int32_t>::const_iterator it = monitor_i.find(key);
        if (it == monitor_i.end())
            return default_value;
        return it->second;
    }

    // ---------------------- binary I/O (single file) ----------------------
    // Format (version=1):
    // magic(8) + version(i32)
    // step(i32) + time(f64) + dt(f64)
    // nref(i32) + residual_ref(f64*nref)
    // residual_max(f64)
    // nd(i32) + [ key(str) + val(f64) ] * nd
    // ni(i32) + [ key(str) + val(i32) ] * ni
    void WriteBinary() const
    {
        std::string path = rundata_path_;
        std::ofstream out(path, std::ios::binary);
        if (!out)
            Fail_("RunData::WriteBinary cannot open " + path);

        const char magic[8] = {'M', 'P', 'C', 'N', 'S', 'R', 'U', 'N'};
        out.write(magic, 8);
        if (!out)
            Fail_("RunData::WriteBinary write magic failed");

        WriteI32_(out, 1); // version

        WriteI32_(out, step);
        WriteF64_(out, time);
        WriteF64_(out, dt);

        WriteI32_(out, (int32_t)residual_ref.size());
        for (size_t i = 0; i < residual_ref.size(); ++i)
            WriteF64_(out, residual_ref[i]);

        WriteF64_(out, residual_max);

        // monitor_d
        WriteI32_(out, (int32_t)monitor_d.size());
        for (std::map<std::string, double>::const_iterator it = monitor_d.begin();
             it != monitor_d.end(); ++it)
        {
            WriteStr_(out, it->first);
            WriteF64_(out, it->second);
        }

        // monitor_i
        WriteI32_(out, (int32_t)monitor_i.size());
        for (std::map<std::string, int32_t>::const_iterator it = monitor_i.begin();
             it != monitor_i.end(); ++it)
        {
            WriteStr_(out, it->first);
            WriteI32_(out, it->second);
        }

        out.close();
    }

    // expected_nref >= 0: enforce residual_ref size match
    void ReadBinary(int expected_nref = -1)
    {
        std::string path = rundata_path_;
        std::ifstream in(path, std::ios::binary);
        if (!in)
            Fail_("RunData::ReadBinary cannot open " + path);

        char magic[8];
        in.read(magic, 8);
        if (!in)
            Fail_("RunData::ReadBinary read magic failed: " + path);

        const char expect[8] = {'M', 'P', 'C', 'N', 'S', 'R', 'U', 'N'};
        for (int i = 0; i < 8; ++i)
        {
            if (magic[i] != expect[i])
                Fail_("RunData::ReadBinary bad magic: " + path);
        }

        int32_t ver = ReadI32_(in);
        if (ver != 1)
            Fail_("RunData::ReadBinary unsupported version: " + path);

        step = ReadI32_(in);
        time = ReadF64_(in);
        dt = ReadF64_(in);

        int32_t nref = ReadI32_(in);
        if (nref < 0)
            nref = 0;

        if (expected_nref >= 0 && nref != expected_nref)
            Fail_("RunData::ReadBinary residual_ref size mismatch");

        residual_ref.assign((size_t)nref, 0.0);
        for (int32_t i = 0; i < nref; ++i)
            residual_ref[(size_t)i] = ReadF64_(in);

        residual_max = ReadF64_(in);

        // monitor_d
        monitor_d.clear();
        int32_t nd = ReadI32_(in);
        if (nd < 0)
            nd = 0;

        for (int32_t i = 0; i < nd; ++i)
        {
            std::string key = ReadStr_(in);
            double val = ReadF64_(in);
            monitor_d[key] = val;
        }

        // monitor_i
        monitor_i.clear();
        int32_t ni = ReadI32_(in);
        if (ni < 0)
            ni = 0;

        for (int32_t i = 0; i < ni; ++i)
        {
            std::string key = ReadStr_(in);
            int32_t val = ReadI32_(in);
            monitor_i[key] = val;
        }

        in.close();
    }

private:
    static void Fail_(const std::string &msg)
    {
        std::fprintf(stderr, "%s\n", msg.c_str());
        std::abort();
    }

    static void WriteI32_(std::ofstream &out, int32_t v)
    {
        out.write(reinterpret_cast<const char *>(&v), sizeof(int32_t));
        if (!out)
            Fail_("RunData write int32 failed");
    }

    static void WriteF64_(std::ofstream &out, double v)
    {
        out.write(reinterpret_cast<const char *>(&v), sizeof(double));
        if (!out)
            Fail_("RunData write f64 failed");
    }

    static void WriteStr_(std::ofstream &out, const std::string &s)
    {
        int32_t n = (int32_t)s.size();
        WriteI32_(out, n);
        if (n > 0)
        {
            out.write(s.data(), n);
            if (!out)
                Fail_("RunData write string bytes failed");
        }
    }

    static int32_t ReadI32_(std::ifstream &in)
    {
        int32_t v = 0;
        in.read(reinterpret_cast<char *>(&v), sizeof(int32_t));
        if (!in)
            Fail_("RunData read int32 failed");
        return v;
    }

    static double ReadF64_(std::ifstream &in)
    {
        double v = 0.0;
        in.read(reinterpret_cast<char *>(&v), sizeof(double));
        if (!in)
            Fail_("RunData read f64 failed");
        return v;
    }

    static std::string ReadStr_(std::ifstream &in)
    {
        int32_t n = ReadI32_(in);
        if (n < 0)
            n = 0;

        std::string s;
        s.resize((size_t)n);

        if (n > 0)
        {
            in.read(&s[0], n);
            if (!in)
                Fail_("RunData read string bytes failed");
        }
        return s;
    }
};