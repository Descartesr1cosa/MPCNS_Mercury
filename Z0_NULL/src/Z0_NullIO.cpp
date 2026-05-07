#include "Z0_NullIO.h"

#include "0_basic/1_MPCNS_Parameter.h"
#include "0_basic/MPI_WRAPPER.h"
#include "1_grid/1_MPCNS_Grid.h"
#include "3_field/2_MPCNS_Field.h"
#include "5_io/IOModule.h"

#include <iostream>

namespace Z0_NULL
{
    void write_tecplot_output(Param &par,
                              Grid &grid,
                              Field &fields,
                              int step,
                              double time)
    {
        IOModule io;
        io.Setup(&par, &grid, &fields, fields.num_fields());
        io.SetTecplotMode(IOModule::TecplotMode::Mixed);
        io.SetTecplotFields({"phi", "U", "V_cart", "tmp"});
        io.SetTecplotFieldComponentNames("U", {"rho", "rho_u", "rho_v", "rho_w", "rho_E"});
        io.SetTecplotFieldComponentNames("V_cart", {"Vx", "Vy", "Vz"});
        io.WriteTecplotBinFile(step, time);

        int myid = 0;
        PARALLEL::mpi_rank(&myid);
        if (myid == 0)
        {
            std::cout << "Z0_NULL Tecplot output written to ./DATA/flow_field####.plt\n"
                      << std::flush;
        }
    }
}
