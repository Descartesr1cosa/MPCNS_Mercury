#include "1_grid/Grid_Trans.h"

GRID_TRANS::Para_data_transfer::Para_data_transfer()
{
    whether_allocate = false;
};

GRID_TRANS::Para_data_transfer::~Para_data_transfer()
{
    Deallocate();
};

/**
 * @brief   释放所分配的空间
 */
void GRID_TRANS::Para_data_transfer::Deallocate()
{
    if (whether_allocate)
    {
        for (int index = 0; index < num_parallel_face; index++)
        {
            // 释放空间
            delete[] buf_send[index];
            delete[] buf_recv[index];
        }
        delete[] buf_send;
        delete[] buf_recv;
        delete[] length;
        delete[] request_s;
        delete[] request_r;
        delete[] status_s;
        delete[] status_r;
        whether_allocate = false;
    }
}

/**
 * @brief 分配用于所有界面虚网格传递的缓存空间，需要网格指针、所需要传递变量的个数
 * @remark 周期边界条件跳过，不开辟空间
 * @remark 不同物理块面耦合界面也会分配缓存空间，但是网格预处理不会传递
 */
void GRID_TRANS::Para_data_transfer::Initial_Allocate_Grid(Grid *grd, int32_t var_num)
{
    //------------------------------------------------------------------
    // 预处理
    grids_info = &(grd->grids);
    depth = grd->ngg + 1;
    if (whether_allocate)
    {
        std::cout << "DATATRANS::Para_data_transfer Has Been Allocated! ! !\n";
        exit(-1);
    }
    whether_allocate = true; // 分配了空间，后期需要释放
    //------------------------------------------------------------------

    //------------------------------------------------------------------
    // 创建临时变量
    int32_t nblock = grd->nblock;
    int32_t index = 0;
    int32_t num_parallel; // 临时变量，表示每个网格块的并行界面个数
    int32_t sub[3], sup[3];
    //------------------------------------------------------------------

    //------------------------------------------------------------------
    // 开始统计一共有多少面
    num_parallel_face = 0;
    //------------------------------
    // 处理统计网格，注意周期边界条件不需要传递
    for (int i = 0; i < nblock; i++)
    {
        num_parallel = (*grids_info)(i).parallel_bc.size();
        for (int j = 0; j < num_parallel; j++)
        {
            // 周期边界条件send recv flag均用奇数表示,故而只计算偶数
            if (fmod((*grids_info)(i).parallel_bc[j].send_flag, 2) == 0)
                num_parallel_face++;
        }
    }
    //------------------------------------------------------------------

    //------------------------------------------------------------------
    // 开缓冲数组空间1
    buf_send = new double *[num_parallel_face]; // 一共 num_parallel_face个指针数组
    buf_recv = new double *[num_parallel_face]; // 一共 num_parallel_face个指针数组
    length = new int32_t[num_parallel_face];
    request_s = new MPI_Request[num_parallel_face];
    request_r = new MPI_Request[num_parallel_face];
    status_s = new MPI_Status[num_parallel_face];
    status_r = new MPI_Status[num_parallel_face];
    //------------------------------------------------------------------

    //------------------------------------------------------------------
    // 开缓冲数组空间2，分配空间
    index = 0;
    for (int i = 0; i < nblock; i++)
    {
        num_parallel = (*grids_info)(i).parallel_bc.size();
        for (int j = 0; j < num_parallel; j++)
        {
            // 周期边界条件send recv flag均用奇数表示,故而只计算偶数
            // 物理量传输，周期边界条件也需要传递；若是网格传输，周期边界跳过
            if (fmod((*grids_info)(i).parallel_bc[j].send_flag, 2) != 0)
                continue;
            for (int k = 0; k < 3; k++)
            {
                sub[k] = abs((*grids_info)(i).parallel_bc[j].sub[k]);
                sup[k] = abs((*grids_info)(i).parallel_bc[j].sup[k]);
            }
            // 求出传输数组的长度
            length[index] = (abs(sub[0] - sup[0]) + 1) * (abs(sub[1] - sup[1]) + 1) * (abs(sub[2] - sup[2]) + 1) * (depth + 1) * var_num;
            buf_send[index] = new double[length[index]];
            buf_recv[index] = new double[length[index]];
            index++;
        }
    }
    //------------------------------------------------------------------
};

/**
 * @brief 分配用于不同物理块【面耦合界面】虚网格几何信息传递的缓存空间，需要网格指针、所需要传递变量的个数
 * @remark Faker.is_multi_phys在判断中被使用了，需要注意构造Faker区域后调用
 */
void GRID_TRANS::Para_data_transfer::Initial_Allocate_Couple_Grid(Grid *grd, int32_t var_num)
{
    //------------------------------------------------------------------
    // 预处理
    grids_info = &(grd->grids);
    depth = grd->ngg;
    if (whether_allocate)
    {
        std::cout << "DATATRANS::Para_data_transfer Has Been Allocated! ! !\n";
        exit(-1);
    }
    whether_allocate = true; // 分配了空间，后期需要释放
    //------------------------------------------------------------------

    //------------------------------------------------------------------
    int32_t nblock = grd->nblock;
    int32_t index = 0;
    int32_t num_parallel; // 临时变量，表示每个网格块的并行界面个数
    int32_t sub[3], sup[3];
    //------------------------------------------------------------------

    //------------------------------------------------------------------
    // 开始统计一共有多少面
    num_parallel_face = 0;
    //------------------------------
    // 处理统计网格，注意只传递不同物理块界面
    for (int i = 0; i < nblock; i++)
    {
        num_parallel = (*grids_info)(i).parallel_bc.size();
        for (int j = 0; j < num_parallel; j++)
        {
            // 只传递不同物理块界面
            if ((*grids_info)(i).parallel_bc[j].Faker.is_multi_phys)
                num_parallel_face++;
        }
    }
    //------------------------------------------------------------------

    //------------------------------------------------------------------
    // 开缓冲数组空间1
    buf_send = new double *[num_parallel_face]; // 一共 num_parallel_face个指针数组
    buf_recv = new double *[num_parallel_face]; // 一共 num_parallel_face个指针数组
    length = new int32_t[num_parallel_face];
    request_s = new MPI_Request[num_parallel_face];
    request_r = new MPI_Request[num_parallel_face];
    status_s = new MPI_Status[num_parallel_face];
    status_r = new MPI_Status[num_parallel_face];
    //------------------------------------------------------------------

    //------------------------------------------------------------------
    // 分配空间2
    index = 0;
    for (int i = 0; i < nblock; i++)
    {
        num_parallel = (*grids_info)(i).parallel_bc.size();
        for (int j = 0; j < num_parallel; j++)
        {
            // 非耦合边界跳过
            if (!(*grids_info)(i).parallel_bc[j].Faker.is_multi_phys)
                continue;
            for (int k = 0; k < 3; k++)
            {
                sub[k] = abs((*grids_info)(i).parallel_bc[j].sub[k]);
                sup[k] = abs((*grids_info)(i).parallel_bc[j].sup[k]);
            }
            // 求出传输数组的长度
            length[index] = (abs(sub[0] - sup[0]) + 1) * (abs(sub[1] - sup[1]) + 1) * (abs(sub[2] - sup[2]) + 1) * (depth + 1) * var_num;
            buf_send[index] = new double[length[index]];
            buf_recv[index] = new double[length[index]];
            index++;
        }
    }
    //------------------------------------------------------------------
}

//=============================================================================================

void GRID_TRANS::Parallel_send_scalar(Parallel_Boundary *bound, int depth, double3D &scalar, double *buf_send)
{
    int32_t index = 0;
    int32_t sub[3], sup[3], inver_Transform[3], dir[3], cycle[3];
    // 注：Transform为ijk到中间状态的映射 inver_Transform为中间状态到ijk的映射
    // 中间状态：边界中定义的以1表示坐标相同的方向（法），以2表示为负的方向， 以0表示剩下的方向（记为中间状态）
    //  inp信息就是利用【中间状态】找到对应关系
    //-------------------------------------------------------------------------
    for (int i = 0; i < 3; i++)
    {
        sub[i] = abs(bound->sub[i]);
        sup[i] = abs(bound->sup[i]);
    }
    //-------------------------------------------------------------------------
    // cycle是用来表征哪一个方向为法向
    // 这里的正负号是为了保证取数据的时候选的是内部的点，direction>0大号面，取负从而获得内部点
    for (int i = 0; i < 3; i++)
        cycle[i] = -bound->cycle[i] * depth;
    // dir表示从sub到sup是增大还是减小，便于循环统一书写
    for (int i = 0; i < 3; i++)
        dir[i] = (sub[i] < sup[i] + cycle[i]) ? 1 : -1;
    // 求出inver_Transform
    for (int i = 0; i < 3; i++)
        inver_Transform[bound->Transform[i]] = i;
    //-------------------------------------------------------------------------
    // 本程序中并行传输数据是以中间状态顺序实现sub到sup的循环
    // 以1,法向方向为最外层循环，从而能够使得边界与虚网格处理分开,默认0为中间层，2为内层循环
    int32_t out, mid, inner;
    out = inver_Transform[1];
    mid = inver_Transform[0];
    inner = inver_Transform[2];
    //--------------------------------------------------------------------------
    // 开始取出数据
    int32_t ijk[3];
    index = 0;
    for (ijk[out] = sub[out]; (ijk[out] - (sup[out] + cycle[out])) * dir[out] <= 0; ijk[out] += dir[out])
    {
        for (ijk[mid] = sub[mid]; (ijk[mid] - (sup[mid] + cycle[mid])) * dir[mid] <= 0; ijk[mid] += dir[mid])
        {
            for (ijk[inner] = sub[inner]; (ijk[inner] - (sup[inner] + cycle[inner])) * dir[inner] <= 0; ijk[inner] += dir[inner])
            {
                buf_send[index] = scalar(ijk[0], ijk[1], ijk[2]);
                index++;
            }
        }
    }
}

void GRID_TRANS::Parallel_recv_scalar(Parallel_Boundary *bound, int depth, double3D &scalar, double *buf_recv)
{
    int32_t index = 0;
    int32_t sub[3], sup[3], inver_Transform[3], dir[3], cycle[3];
    // 注：Transform为ijk到中间状态的映射  inver_Transform为中间状态到ijk的映射
    // 中间状态：边界中定义的以1表示坐标相同的方向（法），以2表示为负的方向， 以0表示剩下的方向（记为中间状态）
    //  inp信息就是利用【中间状态】找到对应关系
    //-------------------------------------------------------------------------
    for (int i = 0; i < 3; i++)
    {
        sub[i] = abs(bound->sub[i]);
        sup[i] = abs(bound->sup[i]);
    }
    //-------------------------------------------------------------------------
    // cycle是用来表征哪一个方向为法向
    // 这里的正负号是为了保证取数据的时候选的是外部虚网格的点，direction>0大号面，取正从而获得外部点
    for (int i = 0; i < 3; i++)
        cycle[i] = bound->cycle[i] * depth;
    // dir表示从sub到sup是增大还是减小，便于循环统一书写
    for (int i = 0; i < 3; i++)
        dir[i] = (sub[i] < sup[i] + cycle[i]) ? 1 : -1;
    // 求出inver_Transform
    for (int i = 0; i < 3; i++)
        inver_Transform[bound->Transform[i]] = i;
    //-------------------------------------------------------------------------
    // 本程序中并行传输数据是以中间状态顺序实现sub到sup的循环
    // 以1,法向方向为最外层循环，从而能够使得边界与虚网格处理分开,默认0为中间层，2为内层循环
    int32_t out, mid, inner;
    out = inver_Transform[1];
    mid = inver_Transform[0];
    inner = inver_Transform[2];
    //-------------------------------------------------------------------------
    // 开始取出数据更新scalar
    int32_t ijk[3];
    index = 0;
    // 两者相加除2 边界为平均
    for (ijk[out] = sub[out]; (ijk[out] - sup[out]) * dir[out] <= 0; ijk[out] += dir[out])
    {
        for (ijk[mid] = sub[mid]; (ijk[mid] - sup[mid]) * dir[mid] <= 0; ijk[mid] += dir[mid])
        {
            for (ijk[inner] = sub[inner]; (ijk[inner] - sup[inner]) * dir[inner] <= 0; ijk[inner] += dir[inner])
            {
                scalar(ijk[0], ijk[1], ijk[2]) = 0.5 * (scalar(ijk[0], ijk[1], ijk[2]) + buf_recv[index]);
                index++;
            }
        }
    }
    // 虚网格直接赋值
    for (ijk[out] = sub[out] + cycle[out] / depth; (ijk[out] - (sup[out] + cycle[out])) * dir[out] <= 0; ijk[out] += dir[out])
    {
        for (ijk[mid] = sub[mid] + cycle[mid] / depth; (ijk[mid] - (sup[mid] + cycle[mid])) * dir[mid] <= 0; ijk[mid] += dir[mid])
        {
            for (ijk[inner] = sub[inner] + cycle[inner] / depth; (ijk[inner] - (sup[inner] + cycle[inner])) * dir[inner] <= 0; ijk[inner] += dir[inner])
            {
                scalar(ijk[0], ijk[1], ijk[2]) = buf_recv[index];
                index++;
            }
        }
    }
}

void GRID_TRANS::Inner_trans_scalar(Inner_Boundary *bound, int depth, double3D &scalar, double3D &scalar_tar)
{
    int32_t sub[3], sup[3], tar_sub[3], tar_sup[3], dir[3], cycle[3], tar_dir[3], tar_cycle[3];
    int32_t Trans[3], tar_Trans[3], my_to_tar[3], tar_inver_Transform[3];
    // 注：Transform为目标块ijk到中间状态的映射
    // 中间状态：边界中定义的以1表示坐标相同的方向（法），以2表示为负的方向， 以0表示剩下的方向（记为中间状态）
    //  inp信息就是利用【中间状态】找到对应关系
    //============================================================================
    for (int i = 0; i < 3; i++)
    {
        Trans[i] = bound->Transform[i];
        tar_Trans[i] = bound->tar_Transform[i];
        sub[i] = abs(bound->sub[i]);
        sup[i] = abs(bound->sup[i]);
        tar_sub[i] = abs(bound->tar_sub[i]);
        tar_sup[i] = abs(bound->tar_sup[i]);

        // cycle是用来表征哪一个方向为法向, 这里的正负号是为了保证取数据的时候选的是内部的点，direction>0大号面，取负从而获得内部点
        cycle[i] = -bound->cycle[i] * depth;
        tar_cycle[i] = bound->tar_cycle[i] * depth;
    }
    for (int i = 0; i < 3; i++)
    {
        tar_inver_Transform[tar_Trans[i]] = i;
    }
    for (int i = 0; i < 3; i++)
        my_to_tar[i] = tar_inver_Transform[Trans[i]];
    //============================================================================
    //-------------------------------------------------------------------------
    // dir表示从sub到sup是增大还是减小，便于循环统一书写
    for (int i = 0; i < 3; i++)
        dir[i] = (sub[i] < sup[i] + cycle[i]) ? 1 : -1;
    //---------------------------------
    for (int i = 0; i < 3; i++)
        tar_dir[i] = (tar_sub[i] < tar_sup[i] + tar_cycle[i]) ? 1 : -1;
    //-------------------------------------------------------------------------
    // 首先处理将本块的内部网格传输到目标块的虚网格
    int32_t ijk[3], tar_ijk[3];
    for (ijk[0] = sub[0] + cycle[0] / depth; (ijk[0] - (sup[0] + cycle[0])) * dir[0] <= 0; ijk[0] += dir[0])
    {
        for (ijk[1] = sub[1] + cycle[1] / depth; (ijk[1] - (sup[1] + cycle[1])) * dir[1] <= 0; ijk[1] += dir[1])
        {
            for (ijk[2] = sub[2] + cycle[2] / depth; (ijk[2] - (sup[2] + cycle[2])) * dir[2] <= 0; ijk[2] += dir[2])
            {
                for (int32_t l = 0; l < 3; l++)
                    tar_ijk[my_to_tar[l]] = tar_sub[my_to_tar[l]] + (ijk[l] - sub[l]) * dir[l] * tar_dir[my_to_tar[l]];
                scalar_tar(tar_ijk[0], tar_ijk[1], tar_ijk[2]) = scalar(ijk[0], ijk[1], ijk[2]);
            }
        }
    }
    //============================================================================
    //-------------------------------------------------------------------------
    // 然后处理将目标块的内部网格传输到本块的虚网格
    for (int i = 0; i < 3; i++)
    {
        cycle[i] = -cycle[i];
        tar_cycle[i] = -tar_cycle[i];
    }
    for (int i = 0; i < 3; i++)
    {
        dir[i] = (sub[i] < sup[i] + cycle[i]) ? 1 : -1;
        tar_dir[i] = (tar_sub[i] < tar_sup[i] + tar_cycle[i]) ? 1 : -1;
    }
    //-------------------------------------------------------------------------
    for (ijk[0] = sub[0] + cycle[0] / depth; (ijk[0] - (sup[0] + cycle[0])) * dir[0] <= 0; ijk[0] += dir[0])
    {
        for (ijk[1] = sub[1] + cycle[1] / depth; (ijk[1] - (sup[1] + cycle[1])) * dir[1] <= 0; ijk[1] += dir[1])
        {
            for (ijk[2] = sub[2] + cycle[2] / depth; (ijk[2] - (sup[2] + cycle[2])) * dir[2] <= 0; ijk[2] += dir[2])
            {
                for (int32_t l = 0; l < 3; l++)
                    tar_ijk[my_to_tar[l]] = tar_sub[my_to_tar[l]] + (ijk[l] - sub[l]) * dir[l] * tar_dir[my_to_tar[l]];
                scalar(ijk[0], ijk[1], ijk[2]) = scalar_tar(tar_ijk[0], tar_ijk[1], tar_ijk[2]);
            }
        }
    }
    //============================================================================
    //-------------------------------------------------------------------------
    // 然后处理将目标块的边界网格与本块的边界网格平均
    double temp_data;
    for (ijk[0] = sub[0]; (ijk[0] - sup[0]) * dir[0] <= 0; ijk[0] += dir[0])
    {
        for (ijk[1] = sub[1]; (ijk[1] - sup[1]) * dir[1] <= 0; ijk[1] += dir[1])
        {
            for (ijk[2] = sub[2]; (ijk[2] - sup[2]) * dir[2] <= 0; ijk[2] += dir[2])
            {
                for (int32_t l = 0; l < 3; l++)
                    tar_ijk[my_to_tar[l]] = tar_sub[my_to_tar[l]] + (ijk[l] - sub[l]) * dir[l] * tar_dir[my_to_tar[l]];
                temp_data = 0.5 * (scalar_tar(tar_ijk[0], tar_ijk[1], tar_ijk[2]) + scalar(ijk[0], ijk[1], ijk[2]));
                scalar_tar(tar_ijk[0], tar_ijk[1], tar_ijk[2]) = temp_data;
                scalar(ijk[0], ijk[1], ijk[2]) = temp_data;
            }
        }
    }
}

//=============================================================================================

void GRID_TRANS::Parallel_send_tensor(Parallel_Boundary *bound, int depth, Phy_Tensor &ten, double *buf_send)
{
    int32_t index = 0;
    int32_t sub[3], sup[3], inver_Transform[3], dir[3], cycle[3];
    int32_t vec_size1 = ten.Getsizeten1(), vec_size2 = ten.Getsizeten2();
    // 注：Transform为ijk到中间状态的映射 inver_Transform为中间状态到ijk的映射
    // 中间状态：边界中定义的以1表示坐标相同的方向（法），以2表示为负的方向， 以0表示剩下的方向（记为中间状态）
    //  inp信息就是利用【中间状态】找到对应关系
    //-------------------------------------------------------------------------
    for (int i = 0; i < 3; i++)
    {
        sub[i] = abs(bound->sub[i]);
        sup[i] = abs(bound->sup[i]);
    }
    //-------------------------------------------------------------------------
    // cycle是用来表征哪一个方向为法向
    // 这里的正负号是为了保证取数据的时候选的是内部的点，direction>0大号面，取负从而获得内部点
    for (int i = 0; i < 3; i++)
        cycle[i] = -bound->cycle[i] * depth;
    // dir表示从sub到sup是增大还是减小，便于循环统一书写
    for (int i = 0; i < 3; i++)
        dir[i] = (sub[i] < sup[i] + cycle[i]) ? 1 : -1;
    // 求出inver_Transform
    for (int i = 0; i < 3; i++)
        inver_Transform[bound->Transform[i]] = i;
    //-------------------------------------------------------------------------
    // 本程序中并行传输数据是以中间状态顺序实现sub到sup的循环
    // 以1,法向方向为最外层循环，从而能够使得边界与虚网格处理分开,默认0为中间层，2为内层循环
    int32_t out, mid, inner;
    out = inver_Transform[1];
    mid = inver_Transform[0];
    inner = inver_Transform[2];
    //--------------------------------------------------------------------------
    // 开始取出数据
    int32_t ijk[3];
    index = 0;
    for (ijk[out] = sub[out]; (ijk[out] - (sup[out] + cycle[out])) * dir[out] <= 0; ijk[out] += dir[out])
        for (ijk[mid] = sub[mid]; (ijk[mid] - (sup[mid] + cycle[mid])) * dir[mid] <= 0; ijk[mid] += dir[mid])
            for (ijk[inner] = sub[inner]; (ijk[inner] - (sup[inner] + cycle[inner])) * dir[inner] <= 0; ijk[inner] += dir[inner])
                for (int ll = 0; ll < vec_size1; ll++)
                    for (int kk = 0; kk < vec_size2; kk++)
                    {
                        buf_send[index] = ten(ijk[0], ijk[1], ijk[2], ll, kk);
                        index++;
                    }
}

void GRID_TRANS::Parallel_flush_recv_scalar(Parallel_Boundary *bound, double3D &Faker_Bnd, int depth, double *buf_recv)
{
    // 将目标块的MPI传过来的缓冲数组传递到Faker区域
    int32_t index = 0;
    int32_t sub[3], sup[3], inver_Transform[3], dir[3], cycle[3];
    // 注：Transform为ijk到中间状态的映射  inver_Transform为中间状态到ijk的映射
    // 中间状态：边界中定义的以1表示坐标相同的方向（法），以2表示为负的方向， 以0表示剩下的方向（记为中间状态）
    //  inp信息就是利用【中间状态】找到对应关系
    //-------------------------------------------------------------------------
    for (int i = 0; i < 3; i++)
    {
        sub[i] = abs(bound->sub[i]);
        sup[i] = abs(bound->sup[i]);
    }
    //-------------------------------------------------------------------------
    // cycle是用来表征哪一个方向为法向
    // 这里的正负号是为了保证取数据的时候选的是外部虚网格的点，direction>0大号面，取正从而获得外部点
    for (int i = 0; i < 3; i++)
        cycle[i] = bound->cycle[i] * depth;
    // cycle[abs(bound->direction) - 1] = (bound->direction > 0) ? depth : -depth;
    // dir表示从sub到sup是增大还是减小，便于循环统一书写
    for (int i = 0; i < 3; i++)
        dir[i] = (sub[i] < sup[i] + cycle[i]) ? 1 : -1;
    // 求出inver_Transform
    for (int i = 0; i < 3; i++)
        inver_Transform[bound->Transform[i]] = i;
    //-------------------------------------------------------------------------
    // 本程序中并行传输数据是以中间状态顺序实现sub到sup的循环
    // 以1,法向方向为最外层循环，从而能够使得边界与虚网格处理分开,默认0为中间层，2为内层循环
    int32_t out, mid, inner;
    out = inver_Transform[1];
    mid = inver_Transform[0];
    inner = inver_Transform[2];
    //-------------------------------------------------------------------------
    // 开始取出数据更新scalar
    int32_t ijk[3];
    int32_t min_sub[3];
    for (int i = 0; i < 3; i++)
        min_sub[i] = fmin(sub[i], sup[i] + cycle[i]);
    index = 0;
    // 虚网格直接赋值
    for (ijk[out] = sub[out]; (ijk[out] - (sup[out] + cycle[out])) * dir[out] <= 0; ijk[out] += dir[out])
    {
        for (ijk[mid] = sub[mid]; (ijk[mid] - (sup[mid] + cycle[mid])) * dir[mid] <= 0; ijk[mid] += dir[mid])
        {
            for (ijk[inner] = sub[inner]; (ijk[inner] - (sup[inner] + cycle[inner])) * dir[inner] <= 0; ijk[inner] += dir[inner])
            {
                Faker_Bnd(ijk[0] - min_sub[0], ijk[1] - min_sub[1], ijk[2] - min_sub[2]) = buf_recv[index];
                index++;
            }
        }
    }
}

void GRID_TRANS::Parallel_flush_recv_tensor(Parallel_Boundary *bound, Phy_Tensor &Faker_Bnd, int depth, double *buf_recv)
{
    // 将目标块的MPI传过来的缓冲数组传递到Faker区域
    int32_t index = 0;
    int32_t sub[3], sup[3], inver_Transform[3], dir[3], cycle[3];
    int32_t vec_size1 = Faker_Bnd.Getsizeten1(), vec_size2 = Faker_Bnd.Getsizeten2();
    // 注：Transform为ijk到中间状态的映射  inver_Transform为中间状态到ijk的映射
    // 中间状态：边界中定义的以1表示坐标相同的方向（法），以2表示为负的方向， 以0表示剩下的方向（记为中间状态）
    //  inp信息就是利用【中间状态】找到对应关系
    //-------------------------------------------------------------------------
    for (int i = 0; i < 3; i++)
    {
        sub[i] = abs(bound->sub[i]);
        sup[i] = abs(bound->sup[i]);
    }
    //-------------------------------------------------------------------------
    // cycle是用来表征哪一个方向为法向
    // 这里的正负号是为了保证取数据的时候选的是外部虚网格的点，direction>0大号面，取正从而获得外部点
    for (int i = 0; i < 3; i++)
        cycle[i] = bound->cycle[i] * depth;
    // cycle[abs(bound->direction) - 1] = (bound->direction > 0) ? depth : -depth;
    // dir表示从sub到sup是增大还是减小，便于循环统一书写
    for (int i = 0; i < 3; i++)
        dir[i] = (sub[i] < sup[i] + cycle[i]) ? 1 : -1;
    // 求出inver_Transform
    for (int i = 0; i < 3; i++)
        inver_Transform[bound->Transform[i]] = i;
    //-------------------------------------------------------------------------
    // 本程序中并行传输数据是以中间状态顺序实现sub到sup的循环
    // 以1,法向方向为最外层循环，从而能够使得边界与虚网格处理分开,默认0为中间层，2为内层循环
    int32_t out, mid, inner;
    out = inver_Transform[1];
    mid = inver_Transform[0];
    inner = inver_Transform[2];
    //-------------------------------------------------------------------------
    // 开始取出数据更新scalar
    int32_t ijk[3];
    int32_t min_sub[3];
    for (int i = 0; i < 3; i++)
        min_sub[i] = fmin(sub[i], sup[i] + cycle[i]);
    index = 0;
    // 虚网格直接赋值
    for (ijk[out] = sub[out]; (ijk[out] - (sup[out] + cycle[out])) * dir[out] <= 0; ijk[out] += dir[out])
        for (ijk[mid] = sub[mid]; (ijk[mid] - (sup[mid] + cycle[mid])) * dir[mid] <= 0; ijk[mid] += dir[mid])
            for (ijk[inner] = sub[inner]; (ijk[inner] - (sup[inner] + cycle[inner])) * dir[inner] <= 0; ijk[inner] += dir[inner])
                for (int ll = 0; ll < vec_size1; ll++)
                    for (int kk = 0; kk < vec_size1; kk++)
                    {
                        Faker_Bnd(ijk[0] - min_sub[0], ijk[1] - min_sub[1], ijk[2] - min_sub[2], ll, kk) = buf_recv[index];
                        index++;
                    }
}

void GRID_TRANS::Inner_flush_scalar(Inner_Boundary *bound, double3D &Faker_Bnd, int depth, double3D &scalar_tar)
{
    // 将scalar_tar的bound范围内的数据拷入Faker_Bnd，注意要在Inner_Boundary建立好后调用
    if (bound->Faker.is_multi_phys == false)
    {
        std::cout << "#Fatal Error: The Inner_Boundary is not a Couple Boundary:" << bound->target_block_name << "<-\t<-" << bound->this_block_name << "\n";
        exit(-1);
    }
    int32_t sub[3], sup[3], tar_sub[3], tar_sup[3], dir[3], cycle[3], tar_dir[3], tar_cycle[3];
    int32_t Trans[3], tar_Trans[3], my_to_tar[3], tar_inver_Transform[3];
    // 注：Transform为目标块ijk到中间状态的映射
    // 中间状态：边界中定义的以1表示坐标相同的方向（法），以2表示为负的方向， 以0表示剩下的方向（记为中间状态）
    //  inp信息就是利用【中间状态】找到对应关系
    //============================================================================
    for (int i = 0; i < 3; i++)
    {
        Trans[i] = bound->Transform[i];
        tar_Trans[i] = bound->tar_Transform[i];
        sub[i] = abs(bound->sub[i]);
        sup[i] = abs(bound->sup[i]);
        tar_sub[i] = abs(bound->tar_sub[i]);
        tar_sup[i] = abs(bound->tar_sup[i]);

        // cycle是用来表征哪一个方向为法向, 这里的正负号是为了保证取数据的时候选的是内部的点，direction>0大号面，取负从而获得内部点
        cycle[i] = bound->cycle[i] * depth;
        tar_cycle[i] = -bound->tar_cycle[i] * depth;
    }
    for (int i = 0; i < 3; i++)
    {
        tar_inver_Transform[tar_Trans[i]] = i;
    }
    for (int i = 0; i < 3; i++)
        my_to_tar[i] = tar_inver_Transform[Trans[i]];
    //============================================================================
    //-------------------------------------------------------------------------
    // dir表示从sub到sup是增大还是减小，便于循环统一书写
    for (int i = 0; i < 3; i++)
        dir[i] = (sub[i] < sup[i] + cycle[i]) ? 1 : -1;
    //---------------------------------
    for (int i = 0; i < 3; i++)
        tar_dir[i] = (tar_sub[i] < tar_sup[i] + tar_cycle[i]) ? 1 : -1;
    //-------------------------------------------------------------------------
    //============================================================================
    //-------------------------------------------------------------------------
    // 将目标块的内部网格传输到本块Faker的虚网格
    int32_t ijk[3], tar_ijk[3], my_sub[3];
    for (int i = 0; i < 3; i++)
        my_sub[i] = fmin(abs(bound->sub[i]), abs(bound->sup[i]));
    if (bound->direction < 0)
        my_sub[abs(bound->direction) - 1] -= depth;
    for (ijk[0] = sub[0]; (ijk[0] - (sup[0] + cycle[0])) * dir[0] <= 0; ijk[0] += dir[0])
    {
        for (ijk[1] = sub[1]; (ijk[1] - (sup[1] + cycle[1])) * dir[1] <= 0; ijk[1] += dir[1])
        {
            for (ijk[2] = sub[2]; (ijk[2] - (sup[2] + cycle[2])) * dir[2] <= 0; ijk[2] += dir[2])
            {
                for (int32_t l = 0; l < 3; l++)
                    tar_ijk[my_to_tar[l]] = tar_sub[my_to_tar[l]] + (ijk[l] - sub[l]) * dir[l] * tar_dir[my_to_tar[l]];
                Faker_Bnd(ijk[0] - my_sub[0], ijk[1] - my_sub[1], ijk[2] - my_sub[2]) = scalar_tar(tar_ijk[0], tar_ijk[1], tar_ijk[2]);
            }
        }
    }
}

void GRID_TRANS::Inner_flush_tensor(Inner_Boundary *bound, Phy_Tensor &Faker_Bnd, int depth, Phy_Tensor &ten_tar)
{
    // 将ten_tar的bound范围内的数据拷入Faker_Bnd，注意要在Inner_Boundary建立好后调用
    if (bound->Faker.is_multi_phys == false)
    {
        std::cout << "#Fatal Error: The Inner_Boundary is not a Couple Boundary:" << bound->target_block_name << "<-\t<-" << bound->this_block_name << "\n";
        exit(-1);
    }
    int32_t sub[3], sup[3], tar_sub[3], tar_sup[3], dir[3], cycle[3], tar_dir[3], tar_cycle[3];
    int32_t Trans[3], tar_Trans[3], my_to_tar[3], tar_inver_Transform[3];
    // 注：Transform为目标块ijk到中间状态的映射
    // 中间状态：边界中定义的以1表示坐标相同的方向（法），以2表示为负的方向， 以0表示剩下的方向（记为中间状态）
    //  inp信息就是利用【中间状态】找到对应关系
    //============================================================================
    for (int i = 0; i < 3; i++)
    {
        Trans[i] = bound->Transform[i];
        tar_Trans[i] = bound->tar_Transform[i];
        sub[i] = abs(bound->sub[i]);
        sup[i] = abs(bound->sup[i]);
        tar_sub[i] = abs(bound->tar_sub[i]);
        tar_sup[i] = abs(bound->tar_sup[i]);

        // cycle是用来表征哪一个方向为法向, 这里的正负号是为了保证取数据的时候选的是内部的点，direction>0大号面，取负从而获得内部点
        cycle[i] = bound->cycle[i] * depth;
        tar_cycle[i] = -bound->tar_cycle[i] * depth;
    }
    for (int i = 0; i < 3; i++)
    {
        tar_inver_Transform[tar_Trans[i]] = i;
    }
    for (int i = 0; i < 3; i++)
        my_to_tar[i] = tar_inver_Transform[Trans[i]];
    //============================================================================
    //-------------------------------------------------------------------------
    // dir表示从sub到sup是增大还是减小，便于循环统一书写
    for (int i = 0; i < 3; i++)
        dir[i] = (sub[i] < sup[i] + cycle[i]) ? 1 : -1;
    //---------------------------------
    for (int i = 0; i < 3; i++)
        tar_dir[i] = (tar_sub[i] < tar_sup[i] + tar_cycle[i]) ? 1 : -1;
    //-------------------------------------------------------------------------
    //============================================================================
    //-------------------------------------------------------------------------
    // 将目标块的内部网格传输到本块Faker的虚网格
    int32_t ijk[3], tar_ijk[3], my_sub[3];
    int32_t vec_num1, vec_num2;
    vec_num1 = ten_tar.Getsizeten1();
    vec_num2 = ten_tar.Getsizeten2();

    for (int i = 0; i < 3; i++)
        my_sub[i] = fmin(abs(bound->sub[i]), abs(bound->sup[i]));
    if (bound->direction < 0)
        my_sub[abs(bound->direction) - 1] -= depth;

    for (ijk[0] = sub[0]; (ijk[0] - (sup[0] + cycle[0])) * dir[0] <= 0; ijk[0] += dir[0])
        for (ijk[1] = sub[1]; (ijk[1] - (sup[1] + cycle[1])) * dir[1] <= 0; ijk[1] += dir[1])
            for (ijk[2] = sub[2]; (ijk[2] - (sup[2] + cycle[2])) * dir[2] <= 0; ijk[2] += dir[2])
            {
                for (int32_t l = 0; l < 3; l++)
                    tar_ijk[my_to_tar[l]] = tar_sub[my_to_tar[l]] + (ijk[l] - sub[l]) * dir[l] * tar_dir[my_to_tar[l]];
                for (int32_t l = 0; l < vec_num1; l++)
                    for (int32_t m = 0; m < vec_num2; m++)
                        Faker_Bnd(ijk[0] - my_sub[0], ijk[1] - my_sub[1], ijk[2] - my_sub[2], l, m) = ten_tar(tar_ijk[0], tar_ijk[1], tar_ijk[2], l, m);
            }
}
//=============================================================================================