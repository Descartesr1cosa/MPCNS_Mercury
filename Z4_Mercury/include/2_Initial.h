#pragma once

#include "3_field/Field.h"

#include "00_Mercury_Const.h"
#include "0_SolverFields.h"
#include "operators/Dipole.h"

class Mercury_Initial
{
public:
    double qinf[5], q_pv_inf[5], B_imf[3]; // SW: conservative quanities, primitive quantities, IMF
    double qinfs[5], q_pv_infs[5];         // Na+: seed initial state

    DipoleField dipB;

    // Remark：Initial here will not be responsible for restart reading；restart has been done by IOModule
    void Initialization(Field *fld, SolverFields &fid)
    {
        Param *par = fld->par;

        // 1) Background state（IMF / Inflow etc.）
        build_background_state(par);

        // 2) Initialization from the beginning: U will be calculated；restart: U will be read by IOModule
        const bool if_continue = par->GetBoo("continue_calc");
        if (!if_continue)
            Common_Initial(fld, fid);

        // 3) No matter from the beginning/restarting，rebuild not writed(saved) dependent fields
        Initial_Badd(fld, fid);
        Initial_Na_And_PhotoRate_(fld, fid);

        if (par->GetInt("myid") == 0)
        {
            std::cout << "********************Finish the Initial Process! !*********************\n\n\n";
            std::cout << "====================================="
                      << "===============================================================\n";
            std::cout << "\t \t Multi-Physics Coupling Numerical Simulation Solver Begins NOW!\n";
            std::cout << "====================================="
                      << "===============================================================\n\n\n";
        }
    }

    void build_background_state(Param *par)
    {
        // ---- Constants ----
        double gamma = par->GetDou_List("constant").data["gamma"];
        double NA = par->GetDou_List("constant").data["NA"];
        double R_uni = par->GetDou_List("constant").data["R_uni"];
        double q_e = par->GetDou_List("constant").data["q_e"];
        double k_Boltz = R_uni / NA;
        double mu_mag = par->GetDou_List("constant").data["mu_mag"];

        List<double> ini = par->GetDou_List("INITIAL");
        List<double> ref = par->GetDou_List("REF");

        // ---- IMF（Tesla）----
        double Bx_phy = ini.data["Bx"];
        double By_phy = ini.data["By"];
        double Bz_phy = ini.data["Bz"];

        double B_ref = ref.data["B_ref"];

        // Nondimensional IMF
        double Bx = Bx_phy / B_ref;
        double By = By_phy / B_ref;
        double Bz = Bz_phy / B_ref;

        // Fluid reference physical quantities
        double L_ref = ref.data["L_ref"];
        double U_ref = ref.data["U"];
        double n_ref = ref.data["n"];
        double T_ref = ref.data["T"];
        double Molecular_mass_ref = ref.data["Molecular_mass"];
        double rho_ref = Molecular_mass_ref / NA * n_ref;

        // ---- nondimensional primitive vars ----
        double rho0 = (ini.data["n"] / n_ref) * (ini.data["Molecular_mass"] / Molecular_mass_ref);

        double c_y = ini.data["c_y"];
        double c_z = ini.data["c_z"];
        double c_x = -std::sqrt(1.0 - c_y * c_y - c_z * c_z);
        double u0 = c_x * ini.data["U"] / U_ref;
        double v0 = c_y * ini.data["U"] / U_ref;
        double w0 = c_z * ini.data["U"] / U_ref;

        double p_ini = ini.data["n"] * k_Boltz * ini.data["T"];
        double p0 = p_ini / (rho_ref * U_ref * U_ref);

        double T0 = ini.data["T"] / T_ref;

        // Solar Wind state
        q_pv_inf[0] = u0;
        q_pv_inf[1] = v0;
        q_pv_inf[2] = w0;
        q_pv_inf[3] = p0;
        q_pv_inf[4] = T0;

        qinf[0] = rho0;
        qinf[1] = rho0 * u0;
        qinf[2] = rho0 * v0;
        qinf[3] = rho0 * w0;
        qinf[4] = 0.5 * rho0 * (u0 * u0 + v0 * v0 + w0 * w0) // Kinetic energy
                  + p0 / (gamma - 1.0);                      // Inertial energy

        // IMF
        B_imf[0] = Bx;
        B_imf[1] = By;
        B_imf[2] = Bz;

        // qinfs[5], q_pv_infs[5];// Na+: seed initial state

        q_pv_infs[0] = 0.0;                   // Na seeds are assumed to be static
        q_pv_infs[1] = 0.0;                   // Na seeds are assumed to be static
        q_pv_infs[2] = 0.0;                   // Na seeds are assumed to be static
        q_pv_infs[3] = p0 * rho_small / 23.0; // Very low background pressure
        q_pv_infs[4] = T0;                    // Temperature is the same as inflow

        qinfs[0] = rho_small * rho0;
        qinfs[1] = 0.0;
        qinfs[2] = 0.0;
        qinfs[3] = 0.0;
        qinfs[4] = q_pv_infs[3] / (gamma - 1.0) + 0.5 * qinfs[0] * (q_pv_infs[0] * q_pv_infs[0] + q_pv_infs[1] * q_pv_infs[1] + q_pv_infs[2] * q_pv_infs[2]);

        dipB.load_from_param(par);
    }

public:
    Mercury_Initial() {};
    ~Mercury_Initial() = default;

private:
    void Common_Initial(Field *fld, SolverFields &fid)
    {
        Param *par = fld->par;
        if (par->GetInt("myid") == 0)
            std::cout << "---->Starting the Initial Process...\n";

        for (int32_t iblock = 0; iblock < fld->num_blocks(); iblock++)
        {
            FieldBlock &UH = fld->field(fid.fid_U_H, iblock);
            if (!UH.is_allocated())
                continue;
            const Int3 &sub = UH.get_lo();
            const Int3 &sup = UH.get_hi();
            for (int i = sub.i; i < sup.i; i++)
                for (int j = sub.j; j < sup.j; j++)
                    for (int k = sub.k; k < sup.k; k++)
                    {
                        for (int32_t ll = 0; ll < UH.descriptor().ncomp; ll++)
                            UH(i, j, k, ll) = qinf[ll];
                    }
        }

        for (int32_t iblock = 0; iblock < fld->num_blocks(); iblock++)
        {
            FieldBlock &UNa = fld->field(fid.fid_U_Na, iblock);
            if (!UNa.is_allocated())
                continue;
            const Int3 &sub = UNa.get_lo();
            const Int3 &sup = UNa.get_hi();
            for (int i = sub.i; i < sup.i; i++)
                for (int j = sub.j; j < sup.j; j++)
                    for (int k = sub.k; k < sup.k; k++)
                    {
                        for (int32_t ll = 0; ll < UNa.descriptor().ncomp; ll++)
                            UNa(i, j, k, ll) = qinfs[ll];
                    }
        }

        for (int32_t iblock = 0; iblock < fld->num_blocks(); iblock++)
        {
            FieldBlock &Ub_xi = fld->field(fid.fid_B.xi, iblock);
            FieldBlock &Ub_et = fld->field(fid.fid_B.eta, iblock);
            FieldBlock &Ub_ze = fld->field(fid.fid_B.zeta, iblock);
            if (!Ub_xi.is_allocated() || !Ub_et.is_allocated() || !Ub_ze.is_allocated())
                continue;
            {
                const Int3 &sub = Ub_xi.get_lo();
                const Int3 &sup = Ub_xi.get_hi();
                for (int i = sub.i; i < sup.i; i++)
                    for (int j = sub.j; j < sup.j; j++)
                        for (int k = sub.k; k < sup.k; k++)
                        {
                            for (int32_t ll = 0; ll < Ub_xi.descriptor().ncomp; ll++)
                                Ub_xi(i, j, k, ll) = 0.0; // induced Magnetic Field
                        }
            }

            {
                const Int3 &sub = Ub_et.get_lo();
                const Int3 &sup = Ub_et.get_hi();
                for (int i = sub.i; i < sup.i; i++)
                    for (int j = sub.j; j < sup.j; j++)
                        for (int k = sub.k; k < sup.k; k++)
                        {
                            for (int32_t ll = 0; ll < Ub_et.descriptor().ncomp; ll++)
                                Ub_et(i, j, k, ll) = 0.0; // induced Magnetic Field
                        }
            }

            {
                const Int3 &sub = Ub_ze.get_lo();
                const Int3 &sup = Ub_ze.get_hi();
                for (int i = sub.i; i < sup.i; i++)
                    for (int j = sub.j; j < sup.j; j++)
                        for (int k = sub.k; k < sup.k; k++)
                        {
                            for (int32_t ll = 0; ll < Ub_ze.descriptor().ncomp; ll++)
                                Ub_ze(i, j, k, ll) = 0.0; // induced Magnetic Field
                        }
            }
        }
    }

    // Build B_add：IMF + local Dipoles
    void Initial_Badd(Field *fld, const SolverFields &fid)
    {
        // IMF
        // Dipoles
        dipB.Build_Badd_FaceFlux(fld->grd, fld, fld->par,
                                 fid.fid_Badd.xi, fid.fid_Badd.eta, fid.fid_Badd.zeta);
    }

    // Build neutral Na density and photoionization rate
    void Initial_Na_And_PhotoRate_(Field *fld, const SolverFields &fid)
    {
        Param *par = fld->par;

        // Photoionization frequencies.
        // sk1: illuminated region, sk2: geometric shadow region.
        const double sk1 = 5.0e-5; // [s^-1]
        const double sk2 = 1.0e-5; // [s^-1]

        // Mercury radius [m]. Coordinates below are assumed nondimensionalized by R_M.
        const double Rm = par->GetDou_List("REF").data["L_ref"];

        const double pi = std::acos(-1.0);
        const double deg = pi / 180.0;

        // -----------------------------
        // Exner et al. (2020), Table 2
        // -----------------------------
        // MIV: radially symmetric
        const double n_miv0 = 7.84e6; // [m^-3]
        const double H_miv = 431.0e3; // [m]

        // TD: latitude = 0 deg, longitude = 0 deg, angular width = 15 deg
        const double n_td0 = 8.86e9; // [m^-3]
        const double H_td = 100.0e3; // [m]
        const double H_theta_td = 15.0 * deg;

        // PSD: latitude = +/-50 deg, longitude = 0 deg, angular width = 20 deg
        const double n_psd0 = 4.06e10; // [m^-3]
        const double H_psd = 232.0e3;  // [m]
        const double H_theta_psd = 20.0 * deg;
        const double lat_psd_n = 50.0 * deg;
        const double lat_psd_s = -50.0 * deg;

        // SP dayside: latitude = +/-80 deg, longitude = 0 deg, angular width = 15 deg
        // SP nightside: latitude = +/-15 deg, longitude = 180 deg, angular width = 10 deg
        const double n_sp0 = 5.67e6; // [m^-3]
        const double H_sp = 748.0e3; // [m]
        const double H_theta_sp_day = 15.0 * deg;
        const double H_theta_sp_night = 10.0 * deg;
        const double lat_sp_day_n = 80.0 * deg;
        const double lat_sp_day_s = -80.0 * deg;
        const double lat_sp_night_n = 15.0 * deg;
        const double lat_sp_night_s = -15.0 * deg;

        // In this code convention x > 0 is dayside.
        // Longitude is measured from +x toward +y.
        const double lon_day = 0.0;
        const double lon_night = pi;

        auto gauss2 = [](double dtheta, double Htheta) -> double
        {
            if (Htheta <= 0.0)
                return 0.0;
            const double a = dtheta / Htheta;
            return std::exp(-(a * a));
        };

        auto clamp = [](double a) -> double
        {
            return std::max(-1.0, std::min(1.0, a));
        };

        // Convert source-center latitude/longitude to local zenith unit vector.
        auto peak_vec = [](double lat, double lon,
                           double &px, double &py, double &pz)
        {
            const double clat = std::cos(lat);
            px = clat * std::cos(lon);
            py = clat * std::sin(lon);
            pz = std::sin(lat);
        };

        // Angular distance between current local zenith direction r_hat
        // and source peak direction.
        auto delta_angle = [&](double rx, double ry, double rz,
                               double lat0, double lon0) -> double
        {
            double px, py, pz;
            peak_vec(lat0, lon0, px, py, pz);

            const double dot = clamp(rx * px + ry * py + rz * pz);
            return std::acos(dot);
        };

        const int nblock = fld->num_blocks();

        for (int ib = 0; ib < nblock; ++ib)
        {
            FieldBlock &Na = fld->field(fid.fid_Na, ib);       // neutral Na density [cm^-3]
            FieldBlock &Photo = fld->field(fid.fid_Photo, ib); // production rate [cm^-3 s^-1]

            if (!Na.is_allocated() || !Photo.is_allocated())
                continue;

            Block &blk = fld->grd->grids(ib);

            // It is okay to initialize all allocated cells.
            // For source integration / RHS addition, use only physical inner or owner cells.
            const Int3 lo = Na.get_lo();
            const Int3 hi = Na.get_hi();

            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const double x = blk.dual_x(i + 1, j + 1, k + 1);
                        const double y = blk.dual_y(i + 1, j + 1, k + 1);
                        const double z = blk.dual_z(i + 1, j + 1, k + 1);

                        const double r_nd = std::sqrt(x * x + y * y + z * z);

                        // Inside/on Mercury's surface: no neutral exosphere source.
                        if (!(r_nd > 1.0) || !std::isfinite(r_nd))
                        {
                            Na(i, j, k, 0) = 0.0;
                            Photo(i, j, k, 0) = 0.0;
                            continue;
                        }

                        // Unit local zenith direction.
                        const double inv_r = 1.0 / r_nd;
                        const double rx = x * inv_r;
                        const double ry = y * inv_r;
                        const double rz = z * inv_r;

                        // Physical altitude above surface [m].
                        const double r_phys = Rm * r_nd;
                        const double dr = r_phys - Rm;

                        // Exner radial factor (R_M / r)^2.
                        const double r_factor = 1.0 / (r_nd * r_nd);

                        // ------------------------
                        // 1) Neutral Na density [m^-3]
                        // ------------------------

                        // MIV: radially symmetric.
                        const double n_miv =
                            n_miv0 * std::exp(-dr / H_miv);

                        // TD: local Gaussian bell centered at lat=0, lon=0.
                        const double d_td =
                            delta_angle(rx, ry, rz, 0.0, lon_day);

                        const double n_td =
                            n_td0 * std::exp(-dr / H_td) *
                            gauss2(d_td, H_theta_td);

                        // PSD: two local Gaussian bells centered at lat=+/-50 deg, lon=0.
                        const double radial_psd =
                            n_psd0 * std::exp(-dr / H_psd);

                        const double d_psd_n =
                            delta_angle(rx, ry, rz, lat_psd_n, lon_day);

                        const double d_psd_s =
                            delta_angle(rx, ry, rz, lat_psd_s, lon_day);

                        const double n_psd =
                            radial_psd *
                            (gauss2(d_psd_n, H_theta_psd) +
                             gauss2(d_psd_s, H_theta_psd));

                        // SP dayside: two local Gaussian bells centered at lat=+/-80 deg, lon=0.
                        const double radial_sp =
                            n_sp0 * std::exp(-dr / H_sp);

                        const double d_sp_day_n =
                            delta_angle(rx, ry, rz, lat_sp_day_n, lon_day);

                        const double d_sp_day_s =
                            delta_angle(rx, ry, rz, lat_sp_day_s, lon_day);

                        const double n_sp_day =
                            radial_sp *
                            (gauss2(d_sp_day_n, H_theta_sp_day) +
                             gauss2(d_sp_day_s, H_theta_sp_day));

                        // SP nightside: two local Gaussian bells centered at lat=+/-15 deg, lon=180 deg.
                        const double d_sp_night_n =
                            delta_angle(rx, ry, rz, lat_sp_night_n, lon_night);

                        const double d_sp_night_s =
                            delta_angle(rx, ry, rz, lat_sp_night_s, lon_night);

                        const double n_sp_night =
                            radial_sp *
                            (gauss2(d_sp_night_n, H_theta_sp_night) +
                             gauss2(d_sp_night_s, H_theta_sp_night));

                        double n_neutral_m3 =
                            r_factor * (n_miv + n_td + n_psd + n_sp_day + n_sp_night);

                        if (!std::isfinite(n_neutral_m3) || n_neutral_m3 < 0.0)
                            n_neutral_m3 = 0.0;

                        // Store neutral Na density in [cm^-3].
                        const double n_neutral_cm3 = n_neutral_m3 * 1.0e-6;
                        Na(i, j, k, 0) = n_neutral_cm3;

                        // ------------------------
                        // 2) Na+ production rate q = nu * n_Na [cm^-3 s^-1]
                        // ------------------------
                        // Geometric shadow: behind the planet and inside the cylindrical shadow.
                        // This assumes sunlight comes from +x toward the origin, with x > 0 dayside.
                        const double rho_cyl = std::sqrt(y * y + z * z);
                        const bool in_shadow = (x < 0.0 && rho_cyl < 1.0);

                        const double nu = in_shadow ? sk2 : sk1;

                        double qmm = nu * n_neutral_cm3;

                        if (!std::isfinite(qmm) || qmm < 0.0)
                            qmm = 0.0;

                        Photo(i, j, k, 0) = qmm;
                    }
        }

        if (par->GetInt("myid") == 0)
        {
            double sl8 = par->GetDou_List("REF").data["L_ref"];
            std::cout << std::scientific
                      << "[NaInit] L_ref/sl8 = " << 1.0 << "\n"
                      << "[NaInit] H_td/RM   = " << H_td / sl8 << "\n"
                      << "[NaInit] H_psd/RM  = " << H_psd / sl8 << "\n"
                      << "[NaInit] H_miv/RM  = " << H_miv / sl8 << "\n"
                      << "[NaInit] H_sp/RM   = " << H_sp / sl8 << "\n";
        }
    }

    // // Obsolete reference implementation, not compiled.
    // // The active implementation above uses sk2 only in Mercury's geometric shadow.
    // // Build neutral Na density and photoionization rate
    // void Initial_Na_And_PhotoRate_(Field *fld, const SolverFields &fid)
    // {
    //     Param *par = fld->par;

    //     // ---- from Fortran common1.f: sk1, sk2 ----
    //     const double sk1 = 5e-5; // par->GetDou("sk1"); // e.g. 5e-5
    //     const double sk2 = 1e-5; // par->GetDou("sk2"); // e.g. 1e-5

    //     // ---- from Fortran common1.f: sl8 = R_M (m) ----
    //     //  REF.L_ref as Mercury radius（m）
    //     const double sl8 = par->GetDou_List("REF").data["L_ref"];

    //     const double pi = std::acos(-1.0);

    //     // ---- Na model constants (from 03_initiation.for) ----
    //     const double n_miv0 = 7.84e6;
    //     const double H_miv = 431.0e3;

    //     const double n_td0 = 8.86e9;
    //     const double H_td = 100.0e3;
    //     const double H_theta_td = 15.0 * pi / 180.0;

    //     const double n_psd0 = 4.06e10;
    //     const double H_psd = 232.0e3;
    //     const double H_theta_psd = 20.0 * pi / 180.0;
    //     const double theta_psd_n = 50.0 * pi / 180.0;
    //     const double theta_psd_s = -50.0 * pi / 180.0;

    //     const double n_sp0 = 5.67e6;
    //     const double H_sp = 748.0e3;
    //     const double H_theta_sp = 15.0 * pi / 180.0;
    //     const double H_theta_sp_night = 10.0 * pi / 180.0;
    //     const double theta_sp_day = 80.0 * pi / 180.0;
    //     const double theta_sp_night = 165.0 * pi / 180.0;

    //     auto gauss2 = [](double d, double H) -> double
    //     {
    //         if (H <= 0.0)
    //             return 0.0;
    //         const double a = d / H;
    //         return std::exp(-(a * a));
    //     };

    //     auto clamp = [](double a)
    //     {
    //         return std::max(-1.0, std::min(1.0, a));
    //     };

    //     auto peak_vec = [](double lat, double lon,
    //                        double &px, double &py, double &pz)
    //     {
    //         px = std::cos(lat) * std::cos(lon);
    //         py = std::cos(lat) * std::sin(lon);
    //         pz = std::sin(lat);
    //     };

    //     const int nblock = fld->num_blocks();

    //     for (int ib = 0; ib < nblock; ++ib)
    //     {
    //         FieldBlock &Na = fld->field(fid.fid_Na, ib);       // (Cell,1) neutral Na density
    //         FieldBlock &Photo = fld->field(fid.fid_Photo, ib); // (Cell,1) photo-ionization rate

    //         // For Fluid
    //         if (!Na.is_allocated() || !Photo.is_allocated())
    //             continue;

    //         // cell-center coordinates (assumed nondimensional in R_M units)
    //         Block &blk = fld->grd->grids(ib);

    //         const Int3 lo = Na.get_lo();
    //         const Int3 hi = Na.get_hi();

    //         for (int i = lo.i; i < hi.i; ++i)
    //             for (int j = lo.j; j < hi.j; ++j)
    //                 for (int k = lo.k; k < hi.k; ++k)
    //                 {
    //                     const double x = blk.dual_x(i + 1, j + 1, k + 1);
    //                     const double y = blk.dual_y(i + 1, j + 1, k + 1);
    //                     const double z = blk.dual_z(i + 1, j + 1, k + 1);

    //                     // hxz0 = r/R_M  (dimensionless)
    //                     const double hxz0 = std::sqrt(x * x + y * y + z * z);

    //                     auto delta_angle = [&](double lat0, double lon0) -> double
    //                     {
    //                         const double r = std::sqrt(x * x + y * y + z * z);
    //                         const double rx = x / r;
    //                         const double ry = y / r;
    //                         const double rz = z / r;

    //                         double px, py, pz;
    //                         peak_vec(lat0, lon0, px, py, pz);

    //                         const double dot = clamp(rx * px + ry * py + rz * pz);
    //                         return std::acos(dot);
    //                     };

    //                     // inside/on surface: set to zero (Fortran: hxz0 <= 1)
    //                     if (!(hxz0 > 1.0))
    //                     {
    //                         Na(i, j, k, 0) = 0.0;
    //                         Photo(i, j, k, 0) = 0.0;
    //                         continue;
    //                     }

    //                     // physical radius r_norm (m)
    //                     const double r_norm = sl8 * hxz0;

    //                     // r_factor = (R_M/r)^2 = (sl8/r_norm)^2 = 1/hxz0^2
    //                     const double r_factor = 1.0 / (hxz0 * hxz0);

    //                     // theta = atan2(sqrt(y^2+z^2), x)  (0..pi)
    //                     const double r_perp = std::sqrt(y * y + z * z);
    //                     const double theta = std::atan2(r_perp, x);

    //                     // ------------------------
    //                     // 1) Neutral Na: roenum_neu
    //                     // ------------------------
    //                     const double dr = (r_norm - sl8); // (m), altitude above surface in meters

    //                     // MIV
    //                     const double n_miv = n_miv0 * std::exp(-dr / H_miv);

    //                     // TD
    //                     const double n_td = n_td0 * std::exp(-dr / H_td) * gauss2(std::fabs(theta), H_theta_td);

    //                     // PSD (north + south)
    //                     const double radial_psd = n_psd0 * std::exp(-dr / H_psd);
    //                     // const double n_psd = radial_psd *
    //                     //                      (gauss2(std::fabs(theta - theta_psd_n), H_theta_psd) + gauss2(std::fabs(theta - theta_psd_s), H_theta_psd));
    //                     const double lat = std::asin(z / hxz0);

    //                     const double n_psd = radial_psd *
    //                                          (gauss2(lat - theta_psd_n, H_theta_psd) + gauss2(lat - theta_psd_s, H_theta_psd));

    //                     // SP day (two lobes: +/- theta_sp_day)
    //                     const double radial_sp = n_sp0 * std::exp(-dr / H_sp);
    //                     const double n_sp_day = radial_sp *
    //                                             (gauss2(std::fabs(theta - theta_sp_day), H_theta_sp) + gauss2(std::fabs(theta + theta_sp_day), H_theta_sp));

    //                     // SP night (only x < 0)
    //                     const double n_sp_night = (x < 0.0)
    //                                                   ? (radial_sp * (gauss2(std::fabs(theta - theta_sp_night), H_theta_sp_night) +
    //                                                                   gauss2(std::fabs(theta + theta_sp_night), H_theta_sp_night)))
    //                                                   : 0.0;

    //                     // roenum_neu in m^-3
    //                     double roenum_m3 = r_factor * (n_miv + n_td + n_psd + n_sp_day + n_sp_night);
    //                     // double roenum_m3 = r_factor * (n_miv + n_td + n_sp_day + n_sp_night);
    //                     // double roenum_m3 = r_factor * (n_td + n_psd + n_sp_day + n_sp_night);
    //                     if (roenum_m3 < 0.0)
    //                         roenum_m3 = 0.0;

    //                     // convert to cm^-3 (Fortran: *1e-6)
    //                     const double roenum_cm3 = roenum_m3 * 1.0e-6;
    //                     Na(i, j, k, 0) = roenum_cm3;

    //                     // ------------------------
    //                     // 2) Photo-ionization rate: qmm
    //                     // ------------------------
    //                     // Fortran's cosxz = |x| / sqrt(x^2+y^2+z^2) = |cos(theta)|
    //                     const double cosxz = std::fabs(x) / hxz0; // in [0,1]

    //                     double qmm = 0.0;
    //                     if (x > 0.0)
    //                     {
    //                         // day side
    //                         qmm = (((sk1 - sk2) * cosxz) + sk2) * roenum_cm3;
    //                     }
    //                     else
    //                     {
    //                         // night side
    //                         qmm = sk2 * (1.0 - cosxz) * roenum_cm3;
    //                     }

    //                     if (qmm < 0.0)
    //                         qmm = 0.0;
    //                     Photo(i, j, k, 0) = qmm;
    //                 }
    //     }
    // }
};
