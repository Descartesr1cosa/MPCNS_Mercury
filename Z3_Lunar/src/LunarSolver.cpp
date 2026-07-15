// Core
#include "1_grid/1_MPCNS_Grid.h"
#include "2_topology/TopologyBuilder.h"
#include "3_field/Field.h"
#include "4_halo/Halo.h"

// Z3_Lunar
#include "LunarSolver.h"
#include "0_LunarFieldCatalog.h"

LunarSolver::LunarSolver(Grid *grd, TOPO::Topology *topo, Field *fld, Halo *halo,
                             Param *par,
                             TOPO::Topology *topo_equiv,
                             HALO_OWNER::EdgeOwnerSyncPattern *edge_owner_pat,
                             METRIC::SingularEdgeRegistry *singular_edges)
    : grd_(grd),
      topo_(topo),
      fld_(fld),
      halo_(halo),
      par_(par),
      topo_equiv_(topo_equiv),
      edge_owner_pat_(edge_owner_pat),
      singular_edges_(singular_edges)
{
    // ---- Cache field ids ----
    fid_.Init(fld_);

    // ---- Build IO Module ----
    constexpr int NRES = 8; // H+ conservative variables and induced magnetic field
    io_.Setup(par_, grd_, fld_, NRES);

    {
        std::vector<std::string> bin_name = {"U_H", "B_xi", "B_eta", "B_zeta"};
        io_.ClearRestartFields();
        io_.SetRestartFields(bin_name);

        // Visualization output: write every variable at grid nodes.  Cell
        // fields are averaged to nodes by IOModule using the synchronized
        // ghost-cell stencil, so duplicated nodes on adjacent block zones see
        // the same two-sided data.  CellAndNode remains available in IOModule
        // for applications that need Tecplot's mixed variable locations.
        io_.SetTecplotOutputMode(IOModule::TecplotMode::AllNode);
        io_.SetTecplotSingularNodeContext(topo_, singular_edges_);
        io_.SetPostDataContext(topo_equiv_, singular_edges_);
        std::vector<std::string> tec_block_name = {}; // 全部物理块输出
        io_.SetTecplotBlocks(tec_block_name);

        std::vector<std::string> plt_name = {"PV_H", "B_cell", "Bind_cell", "J_cell"};
        io_.SetTecplotPhysicalOutputs(plt_name);

        std::string fld_name = "PV_H";
        std::vector<std::string> var_name = {"u_H", "v_H", "w_H", "p_H", "T_H"};
        io_.SetTecplotFieldComponentNames(fld_name, var_name);

        fld_name = "B_cell";
        var_name = {"Bx", "By", "Bz"};
        io_.SetTecplotFieldComponentNames(fld_name, var_name);

        fld_name = "Bind_cell";
        var_name = {"bx_ind", "by_ind", "bz_ind"};
        io_.SetTecplotFieldComponentNames(fld_name, var_name);

        fld_name = "J_cell";
        var_name = {"Jx", "Jy", "Jz"};
        io_.SetTecplotFieldComponentNames(fld_name, var_name);

        // fld_name = "U_plus";
        // var_name = {"Up_u", "Up_v", "Up_w"};
        // io_.SetTecplotFieldComponentNames(fld_name, var_name);
    }

    {
        std::vector<std::string> vtk_name = {
            // cell volume output
            "PV_H",
            "B_cell",
            "Bind_cell",
            "J_cell",
            "divB",

            // face 2-form debug output
            "B_xi",
            "B_eta",
            "B_zeta",
            "Badd_xi",
            "Badd_eta",
            "Badd_zeta",

            // edge 1-form debug output
            "E_xi",
            "E_eta",
            "E_zeta",
            "J_xi",
            "J_eta",
            "J_zeta"};

        io_.SetParaViewFields(vtk_name);
        io_.SetParaViewPath("./DATA/VTK");
        io_.SetParaViewIncludeGhost(false);
    }

    // ---- Calc Constants ----
    calc_physical_constant(par_);

    // ---- Boundary ----
    {
        std::vector<std::string> bnd_fields = fld_->boundary_field_names();
        lunar_bound_.Setup(grd_, fld_, topo_, halo_, par_, bnd_fields);
    }

    // ---- Initialization ----
    if (par_->GetBoo("continue_calc"))
    {
        io_.ReadRestartBinFile();
        io_.ReadRunDataFile();
    }
    initial_.Initialization(fld_, fid_);
    // Also sanitize a restarted state before its first CFL evaluation.
    ApplyDensityFloor_();

    // Apply the one-way lunar-surface state before the first derived-field
    // reconstruction; U_plus then supplies the limited one-sided velocity to
    // the boundary EMF assembly from the first step onward.
    lunar_bound_.Sync("Ucell");
    lunar_bound_.Sync("Badd");
    // Restarted halo/alias values may be stale.  Synchronize Bface before the
    // first derived-field/Jedge evaluation so the initial implicit RHS uses
    // one globally reconciled Bind field.
    lunar_bound_.Sync("Bface");

    // ---- components ----
    control_.Setup(par_);

    run_data_ = &io_.Run();
    runtime_data_ = &io_.Runtime();

    auto count_global_cells = [&]() -> int64_t
    {
        int64_t local = 0;
        for (int ib = 0; ib < fld_->num_blocks(); ++ib)
        {
            auto &U = fld_->field(fid_.fid_U_H, ib); // 任意 cell-centered 场都行
            if (!U.is_allocated())
                continue;
            const Int3 lo = U.inner_lo();
            const Int3 hi = U.inner_hi();
            local += int64_t(hi.i - lo.i) * int64_t(hi.j - lo.j) * int64_t(hi.k - lo.k);
        }
        double local_d = (double)local;
        double global = 0.0;
        PARALLEL::mpi_sum(&local_d, &global, 1);
        return (int64_t)std::llround(global);
    };

    runtime_data_->Begin(*run_data_, par_, count_global_cells());

    ambipolar_control.enabled = par_->GetBoo("is_Ambipolar_Efield");
    // Keep legacy lunar cases Hall-enabled when the new switch is absent.
    hall_enabled_ = !par_->HasBoo("is_Hall_Efield") ||
                    par_->GetBoo("is_Hall_Efield");
    hall_taper_r_min = par_->GetDou("r_min");
    hall_taper_r_max = par_->GetDou("r_max");

    arti_resist_control.eta_max = par_->GetDou("arti_eta_max");
    arti_resist_control.J_range_start = par_->GetDou("J_range_start");
    arti_resist_control.J_range_on = par_->GetDou("J_range_on");
    arti_resist_control.local_enabled = par_->GetBoo("is_local_artificial_resistivity");
    arti_resist_control.local_eta_max = par_->GetDou("local_arti_eta_max");
    arti_resist_control.local_center[0] = par_->GetDou("local_arti_x0");
    arti_resist_control.local_center[1] = par_->GetDou("local_arti_y0");
    arti_resist_control.local_center[2] = par_->GetDou("local_arti_z0");
    arti_resist_control.local_r_decay = par_->GetDou("local_arti_r_decay");
    arti_resist_control.local_r_cutoff = par_->GetDou("local_arti_r_cutoff");

    SetupHallFaceScratch_();
    SetupCellReconstructionWeights_();

    // ---- Optional offline post-data output ----
    // Missing keys deliberately mean disabled, preserving older CASE files.
    post_static_output_enabled_ = par_->GetBoo("post_static_output");
    post_output_path_ = "./DATA_bin";

    post_write_options_.constant_fields = {
        "Badd_xi", "Badd_eta", "Badd_zeta"};
    post_write_options_.normalization = {
        {"length_ref", L_ref},
        {"time_ref", L_ref / U_ref},
        {"density_ref", rho_ref},
        {"velocity_ref", U_ref},
        {"pressure_ref", rho_ref * U_ref * U_ref},
        {"magnetic_field_ref", B_ref},
        {"electric_field_ref", U_ref * B_ref},
        {"current_density_ref", B_ref / (mu0 * L_ref)},
        {"temperature_ref", T_ref}};
    post_write_options_.physical_constants = {
        {"gamma", gamma_}, {"mu0", mu0}, {"molar_mass_H", M_H}, {"particle_mass_H", m_H}};
    post_write_options_.species = {"H"};
    post_write_options_.existing_flow_fields = {
        "U_H", "B_xi", "B_eta", "B_zeta"};

    if (post_static_output_enabled_)
        io_.WritePostStaticData(post_output_path_, post_write_options_);
}

LunarSolver::~LunarSolver() = default;

void LunarSolver::WritePostStaticData(const std::string &output_directory,
                                        POST::WriteOptions options) const
{
    io_.WritePostStaticData(output_directory, std::move(options));
}
