#include "1_grid/1_MPCNS_Grid.h"
#include "1_grid/Grid_Trans.h"
#include <fstream>
#include <math.h>
#include <algorithm>

/**
 * @brief 针对包含内部边界条件的角区进行传值处理
 */
void Grid::MeshTrans_Corner_Inner()
{
    int32_t block, tar_block;

    for (int ii = 0; ii < grids.Getsize1(); ii++)
    {
        Block &iblock = grids(ii);
        // Inner Boundary Corner Process
        for (auto &e : iblock.edge)
        {
            if (e.index != 0)
                continue;

            Inner_Boundary &inner = *(e.inner_bound);
            // if (inner.target_block_name != inner.this_block_name || e.inner_bound->is_period || e.is_singular)
            // 修改2026/01/14 16：04 耦合面也需要传值
            // if (inner.target_block_name != inner.this_block_name || e.inner_bound->is_period)
            if (e.inner_bound->is_period)
                continue;
            block = inner.this_block_num;
            tar_block = inner.tar_block_num;
            GRID_TRANS::Inner_corner_scalar(&e, grids(block).x, ngg + 1, grids(tar_block).x, dimension);
            GRID_TRANS::Inner_corner_scalar(&e, grids(block).y, ngg + 1, grids(tar_block).y, dimension);
            GRID_TRANS::Inner_corner_scalar(&e, grids(block).z, ngg + 1, grids(tar_block).z, dimension);
        }
    }
}

/**
 * @brief 针对包含并行边界条件的角区进行传值处理
 */
void Grid::MeshTrans_Corner_Parallel()
{
    // 确认哪些Para边需要给对方传角区信息，即确定is_corner_send
    //=============================================================================================
    PARALLEL::mpi_barrier();
    GRID_TRANS::Para_data_transfer DATA;
    DATA.Initial_Allocate_Corner_Grid(this, 1);
    PARALLEL::mpi_barrier();

    //===============================================================================
    int32_t index, num_parallel;
    //===============================================================================
    //----------------------------------------------------------------------
    // 收集信息
    index = 0;
    for (int i = 0; i < nblock; i++)
    {
        num_parallel = grids(i).parallel_bc.size();
        for (int j = 0; j < num_parallel; j++)
        {
            // 耦合面不传
            // 修改2026/01/14 16：04 耦合面也需要传值
            // if (grids(i).parallel_bc[j].target_block_name != grids(i).parallel_bc[j].this_block_name)
            //     continue;
            // 周期不传
            if (fmod(grids(i).parallel_bc[j].send_flag, 2) != 0)
                continue;
            // 是否为corner的信息发送方
            GRID_TRANS::Parallel_corner_send_scalar(&grids(i).parallel_bc[j], grids(i).corner_send[j], ngg + 1, grids(i).x, DATA.buf_send[index], dimension);
            index++;
        }
    }
    PARALLEL::mpi_barrier();
    //----------------------------------------------------------------------
    // 发送接受所有数据
    index = 0;
    for (int i = 0; i < nblock; i++)
    {
        num_parallel = grids(i).parallel_bc.size();
        for (int j = 0; j < num_parallel; j++)
        {
            // 修改2026/01/14 16：04 耦合面也需要传值
            // if (grids(i).parallel_bc[j].target_block_name != grids(i).parallel_bc[j].this_block_name || fmod(grids(i).parallel_bc[j].send_flag, 2) != 0)
            if (fmod(grids(i).parallel_bc[j].send_flag, 2) != 0)
                continue;

            PARALLEL::mpi_data_send(grids(i).parallel_bc[j].tar_myid, grids(i).parallel_bc[j].send_flag, DATA.buf_send[index], DATA.length[index], &DATA.request_s[index]);
            // else
            // 	DATA.request_s[index] = MPI_REQUEST_NULL;

            PARALLEL::mpi_data_recv(grids(i).parallel_bc[j].tar_myid, grids(i).parallel_bc[j].rece_flag, DATA.buf_recv[index], DATA.length[index + DATA.num_parallel_face], &DATA.request_r[index]);
            // else
            // 	DATA.request_r[index] = MPI_REQUEST_NULL;

            index++;
        }
    }
    //----------------------------------------------------------------------
    // 等待完成
    PARALLEL::mpi_wait(DATA.num_parallel_face, DATA.request_s, DATA.status_s);
    PARALLEL::mpi_wait(DATA.num_parallel_face, DATA.request_r, DATA.status_r);
    PARALLEL::mpi_barrier();
    //----------------------------------------------------------------------
    // 取出所有数据
    index = 0;
    for (int i = 0; i < nblock; i++)
    {
        num_parallel = grids(i).parallel_bc.size();
        for (int j = 0; j < num_parallel; j++)
        {
            // 修改2026/01/14 16：04 耦合面也需要传值
            // if (grids(i).parallel_bc[j].target_block_name != grids(i).parallel_bc[j].this_block_name || fmod(grids(i).parallel_bc[j].send_flag, 2) != 0)
            if (fmod(grids(i).parallel_bc[j].send_flag, 2) != 0)
                continue;

            GRID_TRANS::Parallel_corner_recv_scalar(&grids(i).parallel_bc[j], grids(i).corner[j], ngg + 1, grids(i).x, DATA.buf_recv[index], dimension);
            index++;
        }
    }
    PARALLEL::mpi_barrier();
    //===============================================================================
    //===============================================================================
    //----------------------------------------------------------------------
    // 收集信息
    index = 0;
    for (int i = 0; i < nblock; i++)
    {
        num_parallel = grids(i).parallel_bc.size();
        for (int j = 0; j < num_parallel; j++)
        {
            // 耦合面不传
            // 修改2026/01/14 16：04 耦合面也需要传值
            // if (grids(i).parallel_bc[j].target_block_name != grids(i).parallel_bc[j].this_block_name)
            //     continue;
            // 周期不传
            if (fmod(grids(i).parallel_bc[j].send_flag, 2) != 0)
                continue;
            // 是否为corner的信息发送方
            GRID_TRANS::Parallel_corner_send_scalar(&grids(i).parallel_bc[j], grids(i).corner_send[j], ngg + 1, grids(i).y, DATA.buf_send[index], dimension);
            index++;
        }
    }
    PARALLEL::mpi_barrier();
    //----------------------------------------------------------------------
    // 发送接受所有数据
    index = 0;
    for (int i = 0; i < nblock; i++)
    {
        num_parallel = grids(i).parallel_bc.size();
        for (int j = 0; j < num_parallel; j++)
        {
            // 修改2026/01/14 16：04 耦合面也需要传值
            // if (grids(i).parallel_bc[j].target_block_name != grids(i).parallel_bc[j].this_block_name || fmod(grids(i).parallel_bc[j].send_flag, 2) != 0)
            if (fmod(grids(i).parallel_bc[j].send_flag, 2) != 0)
                continue;

            PARALLEL::mpi_data_send(grids(i).parallel_bc[j].tar_myid, grids(i).parallel_bc[j].send_flag, DATA.buf_send[index], DATA.length[index], &DATA.request_s[index]);
            // else
            // 	DATA.request_s[index] = MPI_REQUEST_NULL;

            PARALLEL::mpi_data_recv(grids(i).parallel_bc[j].tar_myid, grids(i).parallel_bc[j].rece_flag, DATA.buf_recv[index], DATA.length[index + DATA.num_parallel_face], &DATA.request_r[index]);
            // else
            // 	DATA.request_r[index] = MPI_REQUEST_NULL;

            index++;
        }
    }
    //----------------------------------------------------------------------
    // 等待完成
    PARALLEL::mpi_wait(DATA.num_parallel_face, DATA.request_s, DATA.status_s);
    PARALLEL::mpi_wait(DATA.num_parallel_face, DATA.request_r, DATA.status_r);
    PARALLEL::mpi_barrier();
    //----------------------------------------------------------------------
    // 取出所有数据
    index = 0;
    for (int i = 0; i < nblock; i++)
    {
        num_parallel = grids(i).parallel_bc.size();
        for (int j = 0; j < num_parallel; j++)
        {
            // 修改2026/01/14 16：04 耦合面也需要传值
            // if (grids(i).parallel_bc[j].target_block_name != grids(i).parallel_bc[j].this_block_name || fmod(grids(i).parallel_bc[j].send_flag, 2) != 0)
            if (fmod(grids(i).parallel_bc[j].send_flag, 2) != 0)
                continue;

            GRID_TRANS::Parallel_corner_recv_scalar(&grids(i).parallel_bc[j], grids(i).corner[j], ngg + 1, grids(i).y, DATA.buf_recv[index], dimension);
            index++;
        }
    }
    PARALLEL::mpi_barrier();
    //===============================================================================
    //===============================================================================
    //----------------------------------------------------------------------
    // 收集信息
    index = 0;
    for (int i = 0; i < nblock; i++)
    {
        num_parallel = grids(i).parallel_bc.size();
        for (int j = 0; j < num_parallel; j++)
        {
            // 耦合面不传
            // 修改2026/01/14 16：04 耦合面也需要传值
            // if (grids(i).parallel_bc[j].target_block_name != grids(i).parallel_bc[j].this_block_name)
            //     continue;
            // 周期不传
            if (fmod(grids(i).parallel_bc[j].send_flag, 2) != 0)
                continue;
            // 是否为corner的信息发送方
            GRID_TRANS::Parallel_corner_send_scalar(&grids(i).parallel_bc[j], grids(i).corner_send[j], ngg + 1, grids(i).z, DATA.buf_send[index], dimension);
            index++;
        }
    }
    PARALLEL::mpi_barrier();
    //----------------------------------------------------------------------
    // 发送接受所有数据
    index = 0;
    for (int i = 0; i < nblock; i++)
    {
        num_parallel = grids(i).parallel_bc.size();
        for (int j = 0; j < num_parallel; j++)
        {
            // 修改2026/01/14 16：04 耦合面也需要传值
            // if (grids(i).parallel_bc[j].target_block_name != grids(i).parallel_bc[j].this_block_name || fmod(grids(i).parallel_bc[j].send_flag, 2) != 0)
            if (fmod(grids(i).parallel_bc[j].send_flag, 2) != 0)
                continue;

            PARALLEL::mpi_data_send(grids(i).parallel_bc[j].tar_myid, grids(i).parallel_bc[j].send_flag, DATA.buf_send[index], DATA.length[index], &DATA.request_s[index]);
            // else
            // 	DATA.request_s[index] = MPI_REQUEST_NULL;

            PARALLEL::mpi_data_recv(grids(i).parallel_bc[j].tar_myid, grids(i).parallel_bc[j].rece_flag, DATA.buf_recv[index], DATA.length[index + DATA.num_parallel_face], &DATA.request_r[index]);
            // else
            // 	DATA.request_r[index] = MPI_REQUEST_NULL;

            index++;
        }
    }
    //----------------------------------------------------------------------
    // 等待完成
    PARALLEL::mpi_wait(DATA.num_parallel_face, DATA.request_s, DATA.status_s);
    PARALLEL::mpi_wait(DATA.num_parallel_face, DATA.request_r, DATA.status_r);
    PARALLEL::mpi_barrier();
    //----------------------------------------------------------------------
    // 取出所有数据
    index = 0;
    for (int i = 0; i < nblock; i++)
    {
        num_parallel = grids(i).parallel_bc.size();
        for (int j = 0; j < num_parallel; j++)
        {
            // 修改2026/01/14 16：04 耦合面也需要传值
            // if (grids(i).parallel_bc[j].target_block_name != grids(i).parallel_bc[j].this_block_name || fmod(grids(i).parallel_bc[j].send_flag, 2) != 0)
            if (fmod(grids(i).parallel_bc[j].send_flag, 2) != 0)
                continue;

            GRID_TRANS::Parallel_corner_recv_scalar(&grids(i).parallel_bc[j], grids(i).corner[j], ngg + 1, grids(i).z, DATA.buf_recv[index], dimension);
            index++;
        }
    }
    PARALLEL::mpi_barrier();
    //===============================================================================
    // ===============================================================================
    // 释放空间
    DATA.Deallocate();
    PARALLEL::mpi_barrier();
    //=============================================================================================
}

void Grid::MeshTrans_Corner3D_Inner()
{
    int32_t block, tar_block;

    for (int ii = 0; ii < grids.Getsize1(); ii++)
    {
        Block &iblock = grids(ii);
        // Inner Boundary Corner Process
        for (auto &e : iblock.edge)
        {
            if (e.index != 0)
                continue;

            Inner_Boundary &inner = *(e.inner_bound);
            // if (inner.target_block_name != inner.this_block_name || e.inner_bound->is_period || e.is_singular)
            // 修改2026/01/14 16：04 耦合面也需要传值
            // if (inner.target_block_name != inner.this_block_name || e.inner_bound->is_period)
            if (e.inner_bound->is_period)
                continue;
            block = inner.this_block_num;
            tar_block = inner.tar_block_num;
            GRID_TRANS::Inner_corner_scalar_corner3D(&e, grids(block).x, ngg + 1, grids(tar_block).x, dimension);
            GRID_TRANS::Inner_corner_scalar_corner3D(&e, grids(block).y, ngg + 1, grids(tar_block).y, dimension);
            GRID_TRANS::Inner_corner_scalar_corner3D(&e, grids(block).z, ngg + 1, grids(tar_block).z, dimension);
        }
    }
}

void Grid::MeshTrans_Corner3D_Parallel()
{
    // 确认哪些Para边需要给对方传角区信息，即确定is_corner_send
    //=============================================================================================
    PARALLEL::mpi_barrier();
    GRID_TRANS::Para_data_transfer DATA;
    DATA.Initial_Allocate_Corner3D_Grid(this, 1);
    PARALLEL::mpi_barrier();

    //===============================================================================
    int32_t index, num_parallel;
    //===============================================================================
    //----------------------------------------------------------------------
    // 收集信息
    index = 0;
    for (int i = 0; i < nblock; i++)
    {
        num_parallel = grids(i).parallel_bc.size();
        for (int j = 0; j < num_parallel; j++)
        {
            // 耦合面不传
            // 修改2026/01/14 16：04 耦合面也需要传值
            // if (grids(i).parallel_bc[j].target_block_name != grids(i).parallel_bc[j].this_block_name)
            //     continue;
            // 周期不传
            if (fmod(grids(i).parallel_bc[j].send_flag, 2) != 0)
                continue;
            // 是否为corner的信息发送方
            GRID_TRANS::Parallel_corner3D_send_scalar(&grids(i).parallel_bc[j], grids(i).corner_send[j], ngg + 1, grids(i).x, DATA.buf_send[index], dimension);
            index++;
        }
    }
    PARALLEL::mpi_barrier();
    //----------------------------------------------------------------------
    // 发送接受所有数据
    index = 0;
    for (int i = 0; i < nblock; i++)
    {
        num_parallel = grids(i).parallel_bc.size();
        for (int j = 0; j < num_parallel; j++)
        {
            // 修改2026/01/14 16：04 耦合面也需要传值
            // if (grids(i).parallel_bc[j].target_block_name != grids(i).parallel_bc[j].this_block_name || fmod(grids(i).parallel_bc[j].send_flag, 2) != 0)
            if (fmod(grids(i).parallel_bc[j].send_flag, 2) != 0)
                continue;

            PARALLEL::mpi_data_send(grids(i).parallel_bc[j].tar_myid, grids(i).parallel_bc[j].send_flag, DATA.buf_send[index], DATA.length[index], &DATA.request_s[index]);
            // else
            // 	DATA.request_s[index] = MPI_REQUEST_NULL;

            PARALLEL::mpi_data_recv(grids(i).parallel_bc[j].tar_myid, grids(i).parallel_bc[j].rece_flag, DATA.buf_recv[index], DATA.length[index + DATA.num_parallel_face], &DATA.request_r[index]);
            // else
            // 	DATA.request_r[index] = MPI_REQUEST_NULL;

            index++;
        }
    }
    //----------------------------------------------------------------------
    // 等待完成
    PARALLEL::mpi_wait(DATA.num_parallel_face, DATA.request_s, DATA.status_s);
    PARALLEL::mpi_wait(DATA.num_parallel_face, DATA.request_r, DATA.status_r);
    PARALLEL::mpi_barrier();
    //----------------------------------------------------------------------
    // 取出所有数据
    index = 0;
    for (int i = 0; i < nblock; i++)
    {
        num_parallel = grids(i).parallel_bc.size();
        for (int j = 0; j < num_parallel; j++)
        {
            // 修改2026/01/14 16：04 耦合面也需要传值
            // if (grids(i).parallel_bc[j].target_block_name != grids(i).parallel_bc[j].this_block_name || fmod(grids(i).parallel_bc[j].send_flag, 2) != 0)
            if (fmod(grids(i).parallel_bc[j].send_flag, 2) != 0)
                continue;

            GRID_TRANS::Parallel_corner3D_recv_scalar(&grids(i).parallel_bc[j], grids(i).corner[j], ngg + 1, grids(i).x, DATA.buf_recv[index], dimension);
            index++;
        }
    }
    PARALLEL::mpi_barrier();
    //===============================================================================
    //===============================================================================
    //----------------------------------------------------------------------
    // 收集信息
    index = 0;
    for (int i = 0; i < nblock; i++)
    {
        num_parallel = grids(i).parallel_bc.size();
        for (int j = 0; j < num_parallel; j++)
        {
            // 耦合面不传
            // 修改2026/01/14 16：04 耦合面也需要传值
            // if (grids(i).parallel_bc[j].target_block_name != grids(i).parallel_bc[j].this_block_name)
            //     continue;
            // 周期不传
            if (fmod(grids(i).parallel_bc[j].send_flag, 2) != 0)
                continue;
            // 是否为corner的信息发送方
            GRID_TRANS::Parallel_corner3D_send_scalar(&grids(i).parallel_bc[j], grids(i).corner_send[j], ngg + 1, grids(i).y, DATA.buf_send[index], dimension);
            index++;
        }
    }
    PARALLEL::mpi_barrier();
    //----------------------------------------------------------------------
    // 发送接受所有数据
    index = 0;
    for (int i = 0; i < nblock; i++)
    {
        num_parallel = grids(i).parallel_bc.size();
        for (int j = 0; j < num_parallel; j++)
        {
            // 修改2026/01/14 16：04 耦合面也需要传值
            // if (grids(i).parallel_bc[j].target_block_name != grids(i).parallel_bc[j].this_block_name || fmod(grids(i).parallel_bc[j].send_flag, 2) != 0)
            if (fmod(grids(i).parallel_bc[j].send_flag, 2) != 0)
                continue;

            PARALLEL::mpi_data_send(grids(i).parallel_bc[j].tar_myid, grids(i).parallel_bc[j].send_flag, DATA.buf_send[index], DATA.length[index], &DATA.request_s[index]);
            // else
            // 	DATA.request_s[index] = MPI_REQUEST_NULL;

            PARALLEL::mpi_data_recv(grids(i).parallel_bc[j].tar_myid, grids(i).parallel_bc[j].rece_flag, DATA.buf_recv[index], DATA.length[index + DATA.num_parallel_face], &DATA.request_r[index]);
            // else
            // 	DATA.request_r[index] = MPI_REQUEST_NULL;

            index++;
        }
    }
    //----------------------------------------------------------------------
    // 等待完成
    PARALLEL::mpi_wait(DATA.num_parallel_face, DATA.request_s, DATA.status_s);
    PARALLEL::mpi_wait(DATA.num_parallel_face, DATA.request_r, DATA.status_r);
    PARALLEL::mpi_barrier();
    //----------------------------------------------------------------------
    // 取出所有数据
    index = 0;
    for (int i = 0; i < nblock; i++)
    {
        num_parallel = grids(i).parallel_bc.size();
        for (int j = 0; j < num_parallel; j++)
        {
            // 修改2026/01/14 16：04 耦合面也需要传值
            // if (grids(i).parallel_bc[j].target_block_name != grids(i).parallel_bc[j].this_block_name || fmod(grids(i).parallel_bc[j].send_flag, 2) != 0)
            if (fmod(grids(i).parallel_bc[j].send_flag, 2) != 0)
                continue;

            GRID_TRANS::Parallel_corner3D_recv_scalar(&grids(i).parallel_bc[j], grids(i).corner[j], ngg + 1, grids(i).y, DATA.buf_recv[index], dimension);
            index++;
        }
    }
    PARALLEL::mpi_barrier();
    //===============================================================================
    //===============================================================================
    //----------------------------------------------------------------------
    // 收集信息
    index = 0;
    for (int i = 0; i < nblock; i++)
    {
        num_parallel = grids(i).parallel_bc.size();
        for (int j = 0; j < num_parallel; j++)
        {
            // 修改2026/01/14 16：04 耦合面也需要传值
            // 耦合面不传
            // if (grids(i).parallel_bc[j].target_block_name != grids(i).parallel_bc[j].this_block_name)
            //     continue;
            // 周期不传
            if (fmod(grids(i).parallel_bc[j].send_flag, 2) != 0)
                continue;
            // 是否为corner的信息发送方
            GRID_TRANS::Parallel_corner3D_send_scalar(&grids(i).parallel_bc[j], grids(i).corner_send[j], ngg + 1, grids(i).z, DATA.buf_send[index], dimension);
            index++;
        }
    }
    PARALLEL::mpi_barrier();
    //----------------------------------------------------------------------
    // 发送接受所有数据
    index = 0;
    for (int i = 0; i < nblock; i++)
    {
        num_parallel = grids(i).parallel_bc.size();
        for (int j = 0; j < num_parallel; j++)
        {
            // 修改2026/01/14 16：04 耦合面也需要传值
            // if (grids(i).parallel_bc[j].target_block_name != grids(i).parallel_bc[j].this_block_name || fmod(grids(i).parallel_bc[j].send_flag, 2) != 0)
            if (fmod(grids(i).parallel_bc[j].send_flag, 2) != 0)
                continue;

            PARALLEL::mpi_data_send(grids(i).parallel_bc[j].tar_myid, grids(i).parallel_bc[j].send_flag, DATA.buf_send[index], DATA.length[index], &DATA.request_s[index]);
            // else
            // 	DATA.request_s[index] = MPI_REQUEST_NULL;

            PARALLEL::mpi_data_recv(grids(i).parallel_bc[j].tar_myid, grids(i).parallel_bc[j].rece_flag, DATA.buf_recv[index], DATA.length[index + DATA.num_parallel_face], &DATA.request_r[index]);
            // else
            // 	DATA.request_r[index] = MPI_REQUEST_NULL;

            index++;
        }
    }
    //----------------------------------------------------------------------
    // 等待完成
    PARALLEL::mpi_wait(DATA.num_parallel_face, DATA.request_s, DATA.status_s);
    PARALLEL::mpi_wait(DATA.num_parallel_face, DATA.request_r, DATA.status_r);
    PARALLEL::mpi_barrier();
    //----------------------------------------------------------------------
    // 取出所有数据
    index = 0;
    for (int i = 0; i < nblock; i++)
    {
        num_parallel = grids(i).parallel_bc.size();
        for (int j = 0; j < num_parallel; j++)
        {
            // 修改2026/01/14 16：04 耦合面也需要传值
            // if (grids(i).parallel_bc[j].target_block_name != grids(i).parallel_bc[j].this_block_name || fmod(grids(i).parallel_bc[j].send_flag, 2) != 0)
            if (fmod(grids(i).parallel_bc[j].send_flag, 2) != 0)
                continue;

            GRID_TRANS::Parallel_corner3D_recv_scalar(&grids(i).parallel_bc[j], grids(i).corner[j], ngg + 1, grids(i).z, DATA.buf_recv[index], dimension);
            index++;
        }
    }
    PARALLEL::mpi_barrier();
    //===============================================================================
    // ===============================================================================
    // 释放空间
    DATA.Deallocate();
    PARALLEL::mpi_barrier();
    //=============================================================================================
}
