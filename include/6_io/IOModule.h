#pragma once

#include "0_basic/MPI_WRAPPER.h"
#include "1_grid/1_MPCNS_Grid.h"
#include "2_topology/Topology.h"
#include "3_field/Field.h"
#include "4_halo/Halo.h"
#include "6_io/RunData.h"
#include "6_io/RuntimeMonitor.h"
#include "7_metric/SingularEdgeRegistry.h"
#include "6_io/PostDataWriter.h"

#include <string>
#include <unordered_map>
#include <unordered_set>
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
    // When false (the production default), restart files containing any
    // solver-native Jedge component are rejected. When true, both legacy
    // files without Jedge and debug files with the complete triplet are read.
    void SetRestartDecJedgeAllowed(bool allow) { restart_dec_jedge_allowed_ = allow; }
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
    void WriteTecplotBinary(int step, double time) { WriteTecplotBinFile(step, time); }

    // Tecplot controls
    enum class TecplotMode
    {
        CellAndNode, // x/y/z nodal, node fields nodal, cell fields cell-centered
        AllNode,     // x/y/z nodal, cell fields interpolated to nodes

        CellAsNode, // legacy: cell data written as nodal at dual_x/y/z
        CellToNode, // legacy alias of AllNode
        Mixed       // legacy alias of CellAndNode
    };
    enum class TecplotFormReconstruction
    {
        ToCell,
        ToNode
    };
    void SetTecplotMode(TecplotMode m) { tec_mode_ = m; }
    void SetTecplotOutputMode(TecplotMode m) { SetTecplotMode(m); }
    void SetTecplotFields(const std::vector<std::string> &names) { tec_fields_ = names; }
    void SetTecplotPhysicalOutputs(const std::vector<std::string> &names) { SetTecplotFields(names); }
    void SetTecplotFieldComponentNames(const std::string &field_name,
                                       const std::vector<std::string> &comp_names)
    {
        tec_comp_names_[field_name] = comp_names;
    }
    // 可选：只输出某些 block_name（空=全部）
    void SetTecplotBlock(const std::vector<std::string> &bn) { tec_block_ = bn; }
    void SetTecplotBlocks(const std::vector<std::string> &bn) { SetTecplotBlock(bn); }
    void SetTecplotFormReconstruction(TecplotFormReconstruction mode) { tec_form_reconstruction_ = mode; }
    void SetTecplotSingularNodeContext(const TOPO::Topology *topology,
                                       const METRIC::SingularEdgeRegistry *singular_edges)
    {
        tec_topology_ = topology;
        tec_singular_edges_ = singular_edges;
        tec_singular_stencil_ready_ = false;
    }
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
    const RunData &Run() const { return run_; }
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

    // 6) On-demand offline post-data export. No output occurs until one of
    // these methods is called explicitly.
    void SetPostDataContext(const TOPO::Topology *topology,
                            const METRIC::SingularEdgeRegistry *singular_edges = nullptr)
    {
        post_topology_ = topology;
        post_singular_edges_ = singular_edges;
    }
    void WritePostStaticData(const std::string &output_directory,
                             POST::WriteOptions options = {}) const;
    //=========================================================================

private:
    //=========================================================================
    // ----- binary datas -----
    std::vector<std::string> restart_field_names_;
    std::string restart_path_;
    bool restart_dec_jedge_allowed_{false};
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
    TecplotFormReconstruction tec_form_reconstruction_{TecplotFormReconstruction::ToCell};
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
        enum class Kind
        {
            Coordinate,
            Field,
            ReconstructedForm
        };

        Kind kind = Kind::Field;
        int fid = -1;
        int comp = -1;
        int32_t loc = 0; // 0 nodal, 1 cell-centered (Tecplot ZoneHeader var-location)
        int form_fid_xi = -1;
        int form_fid_eta = -1;
        int form_fid_zeta = -1;
        bool form_is_face = false;
        std::string name;
    };
    struct TecSingularNodeStencil
    {
        int node_gid = -1;
        std::vector<TOPO::EntityKey> local_cells;
    };
    std::vector<TecVar> BuildTecVars_() const; // Build variable list from Field descriptors
    TecplotMode NormalizedTecplotMode_() const;
    bool TecplotBlockSelected_(const std::string &block_name) const;
    bool BuildTecFormTriplet_(const std::string &name,
                              int &fid_xi,
                              int &fid_eta,
                              int &fid_zeta,
                              bool &is_face,
                              std::string &output_name) const;
    bool AddTecFieldOrGroup_(const std::string &name,
                             std::vector<TecVar> &vars,
                             std::unordered_set<std::string> &seen_scalars,
                             std::unordered_set<std::string> &seen_forms) const;
    // ---- evaluators ----
    // coordinate evaluators
    void EvalCoord_Node_(int ib, int i, int j, int k, double &x, double &y, double &z) const;
    void EvalCoord_Cell_(int ib, int i, int j, int k, double &x, double &y, double &z) const;
    double EvalValue_ReconstructedFormAtCell_(const TecVar &tv, int ib, int i, int j, int k) const;
    double EvalValue_ReconstructedFormAtNode_(const TecVar &tv, int ib, int i, int j, int k) const;
    double EvalValue_CellAsNode_(const TecVar &tv, int ib, int i, int j, int k) const;
    double EvalValue_CellToNode_(const TecVar &tv, int ib, int i, int j, int k) const;
    double EvalValue_Mixed_(const TecVar &tv, int ib, int i, int j, int k) const;
    void BuildTecplotSingularNodeStencils_();
    void PrepareTecplotSingularNodeAverages_(const std::vector<TecVar> &vars);
    bool EvalTecplotSingularNodeAverage_(const TecVar &tv,
                                         int ib, int i, int j, int k,
                                         double &value) const;
    static std::uint64_t TecplotFieldComponentKey_(int fid, int comp);

    const TOPO::Topology *tec_topology_ = nullptr;
    const METRIC::SingularEdgeRegistry *tec_singular_edges_ = nullptr;
    bool tec_singular_stencil_ready_ = false;
    std::vector<TecSingularNodeStencil> tec_singular_node_stencils_;
    std::unordered_map<TOPO::EntityKey, std::size_t, TOPO::EntityKey::Hash>
        tec_singular_node_alias_to_stencil_;
    std::unordered_map<std::uint64_t, std::vector<double>>
        tec_singular_node_averages_;
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
    const TOPO::Topology *post_topology_ = nullptr;
    const METRIC::SingularEdgeRegistry *post_singular_edges_ = nullptr;
};
