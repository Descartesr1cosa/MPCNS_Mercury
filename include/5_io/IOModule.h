#pragma once

#include "0_basic/MPI_WRAPPER.h"
#include "1_grid/1_MPCNS_Grid.h"
#include "3_field/Field.h"
#include "4_halo/Halo.h"
#include "5_io/RunData.h"
#include "5_io/RuntimeMonitor.h"

#include <string>
#include <unordered_map>
#include <vector>

class IOModule
{
public:
    void Setup(Param *par, Grid *grd, Field *fld, int nvar);

    //=========================================================================
    // 1) Restart Binary File IO
    // 设置白名单（只输出这些字段）
    void SetRestartFields(const std::vector<std::string> &names);
    void AddRestartField(const std::string &name);
    void ClearRestartFields();
    // 读入写出field的bin文件
    void WriteRestartBinFile(int step, double time);
    void ReadRestartBinFile();
    // binary restart_path_;
    const std::string &Get_RestartPath() const { return restart_path_; }
    void Change_RestartPath(std::string bin_path) { restart_path_ = bin_path; }
    //=========================================================================

    //=========================================================================
    // 2) Tecplot File Output
    void WriteTecplotBinFile(int step, double time);

    // Tecplot controls
    enum class TecplotMode
    {
        CellAsNode, // cell data written as nodal at dual_x/y/z
        CellToNode, // interpolate cell->node and write nodal at x/y/z
        Mixed       // x/y/z nodal, node fields nodal, cell fields cell-centered
    };
    void SetTecplotMode(TecplotMode m) { tec_mode_ = m; }
    void SetTecplotFields(const std::vector<std::string> &names) { tec_fields_ = names; }
    void SetTecplotFieldComponentNames(const std::string &field_name,
                                       const std::vector<std::string> &comp_names)
    {
        tec_comp_names_[field_name] = comp_names;
    }
    // 可选：只输出某些 block_name（空=全部）
    void SetTecplotBlock(const std::vector<std::string> &bn) { tec_block_ = bn; }
    //=========================================================================

    //=========================================================================
    // 3) ParaView / VTK File Output
    void WriteParaViewFile();

    void SetParaViewFields(const std::vector<std::string> &names);
    void ClearParaViewFields();

    void SetParaViewPath(const std::string &path);
    void SetParaViewIncludeGhost(bool flag);
    //=========================================================================

    //=========================================================================
    // 4) Run info / diagnostics
    RunData &Run() { return run_; }
    RuntimeMonitor &Runtime() { return runtime_; }
    // 所有进程读同一个文件
    void ReadRunDataFile();
    // 仅 rank0 调用写出
    void WriteRunDataFile();
    //=========================================================================

    //=========================================================================
    // 5) Output Checkpoint protection and archive
    void BackupDataDirectory(const std::string &backup_dir = "./DATA_backup");
    void MaybeArchiveDataDirectory(int step, double time);
    //=========================================================================

private:
    //=========================================================================
    // ----- binary datas -----
    std::vector<std::string> restart_field_names_;
    std::string restart_path_;
    // ----- binary helpers -----
    static void WriteI32_(std::ofstream &out, int32_t v);
    static void WriteI64_(std::ofstream &out, int64_t v);
    static void WriteF64_(std::ofstream &out, double v);
    static void WriteStr_(std::ofstream &out, const std::string &s);
    static int32_t ReadI32_(std::ifstream &in);
    static int64_t ReadI64_(std::ifstream &in);
    static double ReadF64_(std::ifstream &in);
    static std::string ReadStr_(std::ifstream &in);
    static void Fail_(const std::string &msg);
    static int32_t LocToI32_(StaggerLocation loc);
    static StaggerLocation I32ToLoc_(int32_t v);
    std::vector<int> BuildRestartFids_() const;
    //=========================================================================

    //=========================================================================
    // ----- tecplot datas -----
    std::vector<std::string> tec_fields_;
    std::vector<std::string> tec_block_;
    std::unordered_map<std::string, std::vector<std::string>> tec_comp_names_;
    TecplotMode tec_mode_{TecplotMode::Mixed};
    std::string tecplot_path_;
    // ----- tecplot helpers -----
    static void TecWriteStr_(FILE *fp, const std::string &s);
    static void TecWriteI32_(FILE *fp, int32_t v);
    static void TecWriteF32_(FILE *fp, float v);
    static void TecWriteF64_(FILE *fp, double v);
    static void plt_write_header_(FILE *fp, const std::string &title,
                                  const std::vector<std::string> &var_names);
    static void plt_write_zone_header_(FILE *fp, const std::string &zone_name,
                                       int IMax, int JMax, int KMax,
                                       double sol_time,
                                       const std::vector<int32_t> &var_locs);
    static void plt_write_data_section_header_(FILE *fp, int nvar);
    // ----- tecplot vars -----
    struct TecVar
    {
        int fid = -1;
        int comp = -1;
        int32_t loc = 0; // 0 nodal, 1 cell-centered (Tecplot ZoneHeader var-location)
        std::string name;
    };
    std::vector<TecVar> BuildTecVars_() const; // Build variable list from Field descriptors
    // ---- evaluators ----
    // coordinate evaluators
    void EvalCoord_Node_(int ib, int i, int j, int k, double &x, double &y, double &z) const;
    void EvalCoord_Cell_(int ib, int i, int j, int k, double &x, double &y, double &z) const;
    double EvalValue_CellAsNode_(const TecVar &tv, int ib, int i, int j, int k) const;
    double EvalValue_CellToNode_(const TecVar &tv, int ib, int i, int j, int k) const;
    double EvalValue_Mixed_(const TecVar &tv, int ib, int i, int j, int k) const;
    //=========================================================================

    //=========================================================================
    // ----- paraview / vtk datas -----
    std::vector<std::string> paraview_fields_;
    std::string paraview_path_ = "./DATA/VTK";
    bool paraview_include_ghost_ = false;
    //=========================================================================

    //=========================================================================
    // run time data for residual control and restart
    RunData run_;
    RuntimeMonitor runtime_;
    int myid_;
    //=========================================================================

    //=========================================================================
    // checkpoint backup / archive
    double archive_output_interval_s_ = 0.0;
    double last_archive_wall_s_ = 0.0;

    std::string data_dir_ = "./DATA";
    std::string archive_root_dir_ = "./DATA_archive";

    void CopyCurrentCheckpointToDirectory_(const std::string &dst_dir) const;
    std::string MakeArchiveDirectoryName_(int step, double time) const;
    //=========================================================================

    //=========================================================================
    Param *par_{};
    Grid *grd_{};
    Field *fld_{};
};
