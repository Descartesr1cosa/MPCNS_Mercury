//==============================================================================
//-------------->>>Multi-Physics Coupling Numerical Simulation<<<---------------
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>M P C N S<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//==============================================================================

//==============================================================================
#include "1_grid/1_MPCNS_Grid.h"
#include "0_basic/MPI_WRAPPER.h"
#include "2_topology/2_MPCNS_Topology.h"
#include "3_field/2_MPCNS_Field.h"
#include "3_field/3_MPCNS_Halo.h"
//==============================================================================

//==============================================================================

//==============================================================================

int main(int arg, char **argv)
{
    //=============================================================================================
    // MPI initialization
    int myid;
    PARALLEL::mpi_initial(arg, argv);
    PARALLEL::mpi_rank(&myid);
    //=============================================================================================

    //=============================================================================================
    //--------------------------------------------------------------------------
    // 读入控制参数
    Param *par = new Param;
    par->ReadParam(myid);
    //--------------------------------------------------------------------------
    // 读入网格并作预处理
    Grid *grd = new Grid;
    grd->Grid_Preprocess(par);
    //--------------------------------------------------------------------------
    // 建立topology
    TOPO::Topology topology = TOPO::build_topology(*grd, myid, par->GetInt("dimension"));
    //--------------------------------------------------------------------------
    // 建立Field
    Field *fld = new Field(grd, par);
    //-------------------------------------
    // 加入求解物理场
    fld->register_field(
        FieldDescriptor{"U_", StaggerLocation::Cell, 5, par->GetInt("ngg")});
    fld->register_field(
        FieldDescriptor{"B_xi", StaggerLocation::FaceXi, 1, par->GetInt("ngg")});
    fld->register_field(
        FieldDescriptor{"B_eta", StaggerLocation::FaceEt, 1, par->GetInt("ngg")});
    fld->register_field(
        FieldDescriptor{"B_zeta", StaggerLocation::FaceZe, 1, par->GetInt("ngg")});
    fld->register_field(
        FieldDescriptor{"PV_", StaggerLocation::Cell, 4, par->GetInt("ngg")});
    //--------------------------------------------------------------------------
    // 建立Halo通信
    Halo *hal = new Halo(fld, &topology);
    //=============================================================================================

    hal->exchange_inner("U_");
    std::cout << "Finish inner communication of rank \t" << myid << "\n"
              << std::flush;
    PARALLEL::mpi_barrier();
    hal->exchange_parallel("U_");
    std::cout << "Finish para communication of rank \t" << myid << "\n"
              << std::flush;
    //=============================================================================================
    //--------------------------------------------------------------------------
    // MPI终止
    PARALLEL::mpi_finalize();
    //--------------------------------------------------------------------------
    // 释放所分配的空间，建议按照创建顺序逆序释放
    delete hal;
    delete fld;
    delete par;
    delete grd;
    //=============================================================================================
    return 0;
}