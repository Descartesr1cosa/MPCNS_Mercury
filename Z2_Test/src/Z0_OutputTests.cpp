#include "Z0_OutputTests.h"

#include "0_basic/1_MPCNS_Parameter.h"
#include "0_basic/MPI_WRAPPER.h"
#include "1_grid/1_MPCNS_Grid.h"
#include "3_field/Field.h"
#include "5_io/IOModule.h"

#include <iostream>

namespace Z0
{
    TestResult test_location_output_smoke(Param &par,
                                          Grid &grid,
                                          Field &fields,
                                          int my_rank,
                                          std::ostream &os)
    {
        TestResult result;
        try
        {
            IOModule io;
            io.Setup(&par, &grid, &fields, fields.num_fields());
            io.SetTecplotMode(IOModule::TecplotMode::Mixed);
            io.SetTecplotFields({"phi_cell", "U_cell", "V_cell", "psi_node"});
            io.SetTecplotFieldComponentNames("U_cell", {"rho", "rho_u", "rho_v", "rho_w", "rho_E"});
            io.SetTecplotFieldComponentNames("V_cell", {"Vx", "Vy", "Vz"});
            io.WriteTecplotBinFile(0, 0.0);
            io.SetParaViewFields({"phi_cell", "U_cell", "V_cell", "psi_node",
                                  "E_xi", "E_eta", "E_zeta",
                                  "B_xi", "B_eta", "B_zeta",
                                  "EdgeXi_cart", "EdgeEt_cart", "EdgeZe_cart",
                                  "FaceXi_cart", "FaceEt_cart", "FaceZe_cart"});
            io.SetParaViewIncludeGhost(false);
            io.WriteParaViewFile();
            if (my_rank == 0)
                os << "[Z0][LocationOutput] Tecplot and ParaView smoke writers were called.\n";
        }
        catch (...)
        {
            result.pass = false;
        }
        report_test("LocationOutput", result, os);
        return result;
    }
}
