#pragma once

#include "3_field/FieldArray.h"

// Scratch storage used by the explicit Hall face reconstruction.
struct HallFaceScratchBlock
{
    Vector Ehc;
    Scalar beta;
};

// Number-density data shared by the Hall and two-fluid source terms.
struct NumInfo
{
    double nH_true{0.0};
    double nNa_true{0.0};
    double ne_true{0.0};
    double ne_eff{0.0};

    double chiH{1.0};
    double chiNa{0.0};

    double wH_mhd{0.0};
    double wNa_mhd{0.0};
    double mhd_taper{0.0};
};

struct ResistiveEdgeEMFControl
{
    bool is_Mercury_resistance{false};
    double radial_inner{0.84};
    double radial_outer{1.04};
    double radial_width{0.02};

    double implicit_ksp_rtol{1.0e-8};
    double implicit_ksp_atol{1.0e-12};
    int implicit_ksp_max_it{200};
};

struct ArtificialResistivityControl
{
    double eta_max{0.0};
    double J_range_start{0.0};
    double J_range_on{0.0};

    bool local_enabled{false};
    double local_eta_max{0.0};
    double local_center[3]{0.0, 0.0, 0.0};
    double local_r_decay{1.0};
    double local_r_cutoff{0.0};
};

struct AmbipolarEdgeEMFControl
{
    bool enabled{false};
};
