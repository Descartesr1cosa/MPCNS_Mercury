#pragma once
#include "0_basic/MPI_WRAPPER.h"
#include "1_grid/Grid_Boundary.h"
#include "1_grid/1_MPCNS_Grid.h"

namespace GRID_TRANS
{
    // 结构网格用于传递数据的缓冲数据结构
    class Para_data_transfer
    {
    public:
        double **buf_send; // 所有的块一共 num_parallel_face个指针数组,存储发送数据
        double **buf_recv; // 所有的块一共 num_parallel_face个指针数组,存储接收数据
        int32_t *length;   // 每一个发送面的发送、接受数据长度

        MPI_Request *request_s;
        MPI_Request *request_r;
        MPI_Status *status_s;
        MPI_Status *status_r;

        int32_t num_parallel_face;  // 本进程中一共有多少个并行通信面
        int32_t depth;              // 传递通信面的深度，例如ngg=1，对于计算需要depth=1，第0层会自动传
        Array1D<Block> *grids_info; // 指向本进程网格的指针，以获取网格信息
        bool whether_allocate;      // 标记是否需要释放空间
    public:
        /**
         * @brief 分配用于所有界面虚网格传递的缓存空间，需要网格指针、所需要传递变量的个数
         * @remark 周期边界条件跳过，不开辟空间
         * @remark 不同物理块面耦合界面也会分配缓存空间，但是网格预处理不会传递
         */
        void Initial_Allocate_Grid(Grid *grd, int32_t var_num);

        /**
         * @brief 分配用于不同物理块【面耦合界面】虚网格几何信息传递的缓存空间，需要网格指针、所需要传递变量的个数
         * @remark Faker.is_multi_phys在判断中被使用了，需要注意构造Faker区域后调用
         */
        void Initial_Allocate_Couple_Grid(Grid *grd, int32_t var_num);

        /**
         * @brief 分配用于所有界面网格并行通信传递的缓存空间，需要网格指针
         */
        void Initial_Allocate_Corner_Grid(Grid *grd, int32_t var_num);

        /**
         * @brief 分配用于所有界面3D角区网格并行通信传递的缓存空间，需要网格指针
         */
        void Initial_Allocate_Corner3D_Grid(Grid *grd, int32_t var_num);

        /**
         * @brief   释放所分配的空间
         */
        void Deallocate();

        Para_data_transfer();
        ~Para_data_transfer();
    };

    //=============================================================================================
    /*
     * @brief Parallel_send_scalar与并行传数据无直接关联，但是提供了一个接口，
     *          能够利用输入的Parallel_Boundary信息将具有相同结构（mx my mz ngg）的scalar数组对应Parallel_Boundary
     *          中范围的数据取出存入buf_send中，注意取出方式为Transform的顺序。
     */
    void Parallel_send_scalar(Parallel_Boundary *bound, int depth, double3D &scalar, double *buf_send);

    /*
     * @brief Parallel_recv_scalar与并行传数据无直接关联，但是提供了一个接口，
     *          能够利用输入的Parallel_Boundary信息使具有相同结构（mx my mz ngg）的scalar数组对应Parallel_Boundary
     *          中范围读入存在buf_recv中的数据，注意取出方式为Transform的顺序。
     */
    void Parallel_recv_scalar(Parallel_Boundary *bound, int depth, double3D &scalar, double *buf_recv);

    // 这里的Inner_trans每调用一次就实现了交界面两边的传值过程，【请勿重复调】用否则会对界面第0层（或第max层）造成影响
    void Inner_trans_scalar(Inner_Boundary *bound, int depth, double3D &scalar, double3D &scalar_tar);

    //=============================================================================================

    /*
     * @brief Parallel_send_scalar与并行传数据无直接关联，但是提供了一个接口，
     *          能够利用输入的Parallel_Boundary信息将具有相同结构（mx my mz ngg l）的vec数组对应Parallel_Boundary
     *          中范围的数据取出存入buf_send中，注意取出方式为Transform的顺序。
     */
    void Parallel_send_tensor(Parallel_Boundary *bound, int depth, Phy_Tensor &ten, double *buf_send);

    // 将目标块的MPI传过来的缓冲数组传递到Faker区域
    void Parallel_flush_recv_scalar(Parallel_Boundary *bound, double3D &Faker_Bnd, int depth, double *buf_recv);

    void Parallel_flush_recv_tensor(Parallel_Boundary *bound, Phy_Tensor &Faker_Bnd, int depth, double *buf_recv);

    // 将目标块的scalar数组传递到Faker区域
    void Inner_flush_scalar(Inner_Boundary *bound, double3D &Faker_Bnd, int depth, double3D &scalar_tar);
    void Inner_flush_tensor(Inner_Boundary *bound, Phy_Tensor &Faker_Bnd, int depth, Phy_Tensor &ten_tar);

    //=============================================================================================

    // Inner 角区的传值
    void Inner_corner_scalar(Edge *edge, double3D &scalar, int depth, double3D &scalar_tar, int32_t dimension);

    /*
     * @brief Parallel_corner_send_scalar 是将包含角区的数据取出按顺序存入buf_send
     */
    void Parallel_corner_send_scalar(Parallel_Boundary *bound, std::vector<Edge> &bound_corner_send, int depth, double3D &scalar, double *buf_send, int32_t dimension);
    /*
     * @brief Parallel_corner_send_scalar 是将包含角区的数据从buf_recv中取出
     */
    void Parallel_corner_recv_scalar(Parallel_Boundary *bound, std::vector<Edge> &bound_corner, int depth, double3D &scalar, double *buf_recv, int32_t dimension);

    //=============================================================================================

    void Parallel_corner3D_send_scalar(Parallel_Boundary *bound, std::vector<Edge> &bound_corner_send, int depth, double3D &scalar, double *buf_send, int32_t dimension);

    void Parallel_corner3D_recv_scalar(Parallel_Boundary *bound, std::vector<Edge> &bound_corner, int depth, double3D &scalar, double *buf_recv, int32_t dimension);

    void Inner_corner_scalar_corner3D(Edge *edge, double3D &scalar, int depth, double3D &scalar_tar, int32_t dimension);
}