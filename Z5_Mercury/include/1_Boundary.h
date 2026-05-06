#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <functional>
#include <stdexcept>
#include <utility>
#include <iosfwd>
#include <map>

#include "6_boundary/Boundary.h"   // BoundaryCore, BOUND::PhysicalRegion, CouplingBufferBlock...
#include "4_halo/1_MPCNS_Halo.h"   // Halo, HaloLevel
#include "3_field/2_MPCNS_Field.h" // Field, FieldBlock, StaggerLocation
// Grid/Topology/Param forward decl
class Grid;
namespace TOPO
{
    class Topology;
}
class Param;

class MercuryBoundary
{
public:
    MercuryBoundary() = default;

    // ----------------------------- Lifecycle --------------------------------
    // bind pointers + initialize BoundaryCore patterns cache (per location)
    void Setup(Grid *grd, Field *fld, TOPO::Topology *topo, Halo *halo, Param *par, const std::vector<std::string> &boundary_fields);

    // ----------------------------- Run-time API ------------------------------
    void Sync(const std::string &group_name);
    void DescribeGroups(std::ostream &os) const;

private:
    // ----------------------------- Internal types ----------------------------
    using PhysicalHandler = std::function<void(FieldBlock &, Field *, const BOUND::PhysicalRegion &, int ngh)>;
    using CouplingHandler = std::function<void(FieldBlock &, Field *, CouplingBufferBlock &,
                                               const std::string &src, const std::string &dst,
                                               const std::string &tag)>;

    struct BCState
    {
        double q_pv_inf[5], q_pv_infs[5], qinf[5], qinfs[5], B_imf[3];
        double gamma{0.0};
        double dt{0.0};
        double M_H{0.0};
        double M_Na{0.0};
        double state_coeff_H{0.0};
        double state_coeff_Na{0.0};
        double CFL{0.0};
    };

    struct BoundGroup
    {
        std::string name;
        std::vector<std::string> fields;
        std::map<std::pair<std::string, std::string>, std::vector<int32_t>> fields_cids;

        bool do_coupling = false; // whether to apply coupling before physical BC
        bool do_physical = true;  // whether to apply physical BC
        bool do_halo = true;      // whether to halo-exchange

        HaloLevel halo_level = HaloLevel::Vertex;
        int ngh = 0; // reserved: boundary ngh (if BoundaryCore later supports per-call ngh)

        // directed pairs for coupling stage (e.g., {{"Solid","Fluid"},{"Fluid","Solid"}})
        std::vector<std::pair<std::string, std::string>> coupling_pairs;
    };

private:
    // pointers
    Grid *grd_{nullptr};
    Field *fld_{nullptr};
    TOPO::Topology *topo_{nullptr};
    Halo *halo_{nullptr};
    Param *par_{nullptr};
    //  boundary core
    BoundaryCore bound_;

    // bookkeeping
    std::vector<std::string> boundary_fields_;
    std::unordered_map<std::string, BoundGroup> groups_;

    // flag whether this is built successfully
    bool built_{false};

    BCState bc_state_;

private:
    // ---------------------------- Boundary build ----------------------------
    // 0) BCstate
    void InitBCStateFromParam_();

    // 1) install all physical/coupling handlers (table-driven inside this class)
    void InstallHandlers();
    void InstallPhysicalHandlers_();
    void InstallDefaultPhysicalHandlers_();
    void InstallFarfieldPhysicalHandlers_();
    void InstallCoupledPhysicalHandlers_();
    void InstallPolePhysicalHandlers_();
    void InstallCouplingHandlers_();

    // 2) install default groups (Ucell/Bface/Baddface/...); user can extend in cpp later
    void InstallDefaultGroups();

    // 3) finalize: build halo patterns once, and optionally check handler completeness
    void Build(bool strict_check = true);

    // ----------------------------- Group ops (build-time) --------------------
    // Define a sync group (Ucell, Bface, Baddface, Eedge, Jedge, ...).
    void AddGroup(const BoundGroup &g);
    void AddStandardGroup_(const std::string &name,
                           const std::vector<std::string> &fields,
                           bool do_coupling,
                           bool do_physical = true,
                           bool do_halo = true,
                           HaloLevel halo_level = HaloLevel::Vertex);
    // void CheckSetupOrAbort_(const char *where) const;

    // handler registration (build-time)
    void RegisterPhysical_(const std::string &field, const std::string &region, PhysicalHandler h);
    void RegisterCoupling_(const std::string &src, const std::string &dst,
                           StaggerLocation loc,
                           const std::string &channel_tag, const std::string &dst_field,
                           CouplingHandler h);

    // ----------------------------- Sync pipeline (run-time) ------------------
    void Sync_(const BoundGroup &g);

    // -------------------------------  handlers -------------------------------
    void BC_UH_Farfield_H_(FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh);
    void BC_UH_Farfield_Na_(FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh);
    void BC_Solid_Surface_(FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh);
    void BC_Solid_Surface_Eface_(FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh);
    void BC_Solid_Surface_Eedge_(FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh);
    void BC_Pole_Eedge_(FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh);
    void BC_Pole_Eedge_Zero(FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh);
    void BC_Pole_Cell_(FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh);
    void BC_Pole_Eedge_RegulateKAndCopyGhost_(FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh);

    void BC_Pole_Bface_Collapse_(FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh);
    void BC_Pole_Bcell_Collapse_(FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh);
    void BC_Pole_Jcell_Collapse_(FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh);
    void BC_Solid_Surface_Jcell(FieldBlock &U, Field *fld, const BOUND::PhysicalRegion &r, int ngh);
};
