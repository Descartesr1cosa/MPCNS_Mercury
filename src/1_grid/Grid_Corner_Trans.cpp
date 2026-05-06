#include "1_grid/Grid_Trans.h"

namespace GRID_TRANS
{
    void Inner_corner_scalar(Edge *edge, double3D &scalar, int depth, double3D &scalar_tar, int32_t dimension)
    {
        int32_t sub[3], sup[3], tar_sub[3], tar_sup[3], dir[3], cycle[3], tar_dir[3], tar_cycle[3];
        int32_t Trans[3], tar_Trans[3], my_to_tar[3], tar_inver_Transform[3];
        bool is_3D = (dimension == 3);
        Inner_Boundary *bound = edge->inner_bound;
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
        //============================================================================
        //-------------------------------------------------------------------------
        // 处理将目标块的内部网格传输到本块的虚网格
        int32_t ijk[3], tar_ijk[3];
        int32_t range_sub[3], range_sup[3], cycle_dir[3];
        for (int i = 0; i < 3; i++)
        {
            range_sub[i] = abs(edge->sub[i]);
            range_sup[i] = abs(edge->sup[i]);
        }
        if (edge->direction1 > 0)
        {
            range_sup[abs(edge->direction1) - 1] += depth;
            range_sub[abs(edge->direction1) - 1] += 1;
        }
        else
        {
            range_sup[abs(edge->direction1) - 1] -= depth;
            range_sub[abs(edge->direction1) - 1] -= 1;
        }
        if (edge->inner_bound->direction > 0)
        {
            range_sup[abs(edge->inner_bound->direction) - 1] += depth;
            range_sub[abs(edge->inner_bound->direction) - 1] += 1;
        }
        else
        {
            range_sup[abs(edge->inner_bound->direction) - 1] -= depth;
            range_sub[abs(edge->inner_bound->direction) - 1] -= 1;
        }

        for (int i = 0; i < 3; i++)
            cycle_dir[i] = (range_sub[i] <= range_sup[i]) ? 1 : -1;

        //-------------------------------------------------------------------------
        for (ijk[0] = range_sub[0]; (ijk[0] - range_sup[0]) * cycle_dir[0] <= 0; ijk[0] += cycle_dir[0])
            for (ijk[1] = range_sub[1]; (ijk[1] - range_sup[1]) * cycle_dir[1] <= 0; ijk[1] += cycle_dir[1])
                for (ijk[2] = range_sub[2]; (ijk[2] - range_sup[2]) * cycle_dir[2] <= 0; ijk[2] += cycle_dir[2])
                {
                    for (int32_t l = 0; l < 3; l++)
                        tar_ijk[my_to_tar[l]] = tar_sub[my_to_tar[l]] + (ijk[l] - sub[l]) * dir[l] * tar_dir[my_to_tar[l]];
                    scalar(ijk[0], ijk[1], ijk[2]) = scalar_tar(tar_ijk[0], tar_ijk[1], tar_ijk[2]);
                }
    }

    void Inner_corner_scalar_corner3D(Edge *edge, double3D &scalar, int depth, double3D &scalar_tar, int32_t dimension)
    {
        int32_t sub[3], sup[3], tar_sub[3], tar_sup[3];
        int32_t dir[3], cycle[3], tar_dir[3], tar_cycle[3];
        int32_t Trans[3], tar_Trans[3], my_to_tar[3], tar_inver_Transform[3];
        bool is_3D = (dimension == 3);
        if (!is_3D)
            return;
        Inner_Boundary *bound = edge->inner_bound;

        //============================================================================
        // Transform 等初始化，与原函数一致
        for (int i = 0; i < 3; i++)
        {
            Trans[i] = bound->Transform[i];
            tar_Trans[i] = bound->tar_Transform[i];
            sub[i] = abs(bound->sub[i]);
            sup[i] = abs(bound->sup[i]);
            tar_sub[i] = abs(bound->tar_sub[i]);
            tar_sup[i] = abs(bound->tar_sup[i]);

            cycle[i] = bound->cycle[i] * depth;
            tar_cycle[i] = -bound->tar_cycle[i] * depth;
        }

        for (int i = 0; i < 3; i++)
            tar_inver_Transform[tar_Trans[i]] = i;
        for (int i = 0; i < 3; i++)
            my_to_tar[i] = tar_inver_Transform[Trans[i]];

        //============================================================================
        // 正反方向
        for (int i = 0; i < 3; i++)
            dir[i] = (sub[i] < sup[i] + cycle[i]) ? 1 : -1;
        for (int i = 0; i < 3; i++)
            tar_dir[i] = (tar_sub[i] < tar_sup[i] + tar_cycle[i]) ? 1 : -1;

        //============================================================================
        // 计算范围
        int32_t ijk[3], tar_ijk[3];
        int32_t range_sub[3], range_sup[3], cycle_dir[3];

        for (int i = 0; i < 3; i++)
        {
            range_sub[i] = abs(edge->sub[i]);
            range_sup[i] = abs(edge->sup[i]);
        }

        // 第一个法向面扩展
        if (edge->direction1 > 0)
        {
            range_sup[abs(edge->direction1) - 1] += depth;
            range_sub[abs(edge->direction1) - 1] += 1;
        }
        else
        {
            range_sup[abs(edge->direction1) - 1] -= depth;
            range_sub[abs(edge->direction1) - 1] -= 1;
        }

        // 第二个法向面扩展
        if (edge->inner_bound->direction > 0)
        {
            range_sup[abs(edge->inner_bound->direction) - 1] += depth;
            range_sub[abs(edge->inner_bound->direction) - 1] += 1;
        }
        else
        {
            range_sup[abs(edge->inner_bound->direction) - 1] -= depth;
            range_sub[abs(edge->inner_bound->direction) - 1] -= 1;
        }

        //=============================== 新增：第三方向扩展 ============================
        // 找到不在前两个法向中的第三方向
        int32_t dir1 = abs(edge->direction1) - 1;
        int32_t dir2 = abs(edge->inner_bound->direction) - 1;
        int32_t dir3 = 3 - dir1 - dir2;

        // 沿第三方向扩展 depth
        if (range_sup[dir3] > range_sub[dir3])
        {
            range_sup[dir3] += depth;
            range_sub[dir3] -= depth;
        }
        else
        {
            range_sup[dir3] -= depth;
            range_sub[dir3] += depth;
        }

        //==========================================================================

        for (int i = 0; i < 3; i++)
            cycle_dir[i] = (range_sub[i] <= range_sup[i]) ? 1 : -1;

        //============================================================================
        // 传递数据
        for (ijk[0] = range_sub[0]; (ijk[0] - range_sup[0]) * cycle_dir[0] <= 0; ijk[0] += cycle_dir[0])
            for (ijk[1] = range_sub[1]; (ijk[1] - range_sup[1]) * cycle_dir[1] <= 0; ijk[1] += cycle_dir[1])
                for (ijk[2] = range_sub[2]; (ijk[2] - range_sup[2]) * cycle_dir[2] <= 0; ijk[2] += cycle_dir[2])
                {
                    for (int32_t l = 0; l < 3; l++)
                        tar_ijk[my_to_tar[l]] = tar_sub[my_to_tar[l]] + (ijk[l] - sub[l]) * dir[l] * tar_dir[my_to_tar[l]];
                    scalar(ijk[0], ijk[1], ijk[2]) = scalar_tar(tar_ijk[0], tar_ijk[1], tar_ijk[2]);
                }
    }

    /**
     * @brief 分配用于所有界面网格并行通信传递的缓存空间，需要网格指针
     */
    void Para_data_transfer::Initial_Allocate_Corner_Grid(Grid *grd, int32_t var_num)
    {
        //------------------------------------------------------------------
        // 预处理
        grids_info = &(grd->grids);
        depth = grd->ngg + 1;

        if (whether_allocate)
        {
            std::cout << "PARALLEL::Para_data_transfer Has Been Allocated! ! !\n";
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
        // 设置空间，清空对应para_bound的corner corner_send
        for (int i = 0; i < nblock; i++)
        {
            num_parallel = (*grids_info)(i).parallel_bc.size();
            (*grids_info)(i).corner_send.resize(num_parallel);
            (*grids_info)(i).corner.resize(num_parallel);

            for (int j = 0; j < num_parallel; j++)
            {
                // 周期边界条件send recv flag均用奇数表示
                // 同一物理场才传
                // 修改2026/01/14 16：04 耦合面也需要传值
                // if (fmod((*grids_info)(i).parallel_bc[j].send_flag, 2) == 0 && (*grids_info)(i).parallel_bc[j].this_block_name == (*grids_info)(i).parallel_bc[j].target_block_name)
                if (fmod((*grids_info)(i).parallel_bc[j].send_flag, 2) == 0)
                {
                    (*grids_info)(i).corner_send[j].clear();
                    (*grids_info)(i).corner[j].clear();
                }
            }
        }
        // 处理统计每个面的corner，注意周期边界条件不需要传递网格
        //--------------------------------------------------
        // 寻找para对应的index
        auto find_para_index = [&](Block *blk, Parallel_Boundary *p) -> int
        {
            std::vector<Parallel_Boundary> &para = blk->parallel_bc;
            for (int i = 0; i < para.size(); i++)
            {
                if (para[i].sub[0] == p->sub[0] && para[i].sub[1] == p->sub[1] && para[i].sub[2] == p->sub[2] && para[i].sup[0] == p->sup[0] && para[i].sup[1] == p->sup[1] && para[i].sup[2] == p->sup[2])
                    return i;
            }
            std::cout << "Fata error, can not find the para_boundary ! ! !\n";
            return -1;
        };
        for (int i = 0; i < nblock; i++)
        {
            for (Edge &e : (*grids_info)(i).edge)
            {
                Edge e_temp(e);
                if (e.index != 1)
                    continue;
                // 修改2026/01/14 16：04 耦合面也需要传值
                // if (e.para_bound->target_block_name != e.para_bound->this_block_name)
                //     continue;
                int j = find_para_index(&(*grids_info)(i), e.para_bound);
                (*grids_info)(i).corner[j].push_back(e_temp);
            }
        }
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
                // 周期边界条件send recv flag均用奇数表示
                // 同一物理场才传
                // 修改2026/01/14 16：04 耦合面也需要传值
                // if (fmod((*grids_info)(i).parallel_bc[j].send_flag, 2) == 0 && (*grids_info)(i).parallel_bc[j].this_block_name == (*grids_info)(i).parallel_bc[j].target_block_name)
                if (fmod((*grids_info)(i).parallel_bc[j].send_flag, 2) == 0)
                    num_parallel_face++;
            }
        }
        //------------------------------
        //------------------------------------------------------------------
        // 开缓冲数组空间1
        buf_send = new double *[num_parallel_face];  // 一共 num_parallel_face个指针数组
        buf_recv = new double *[num_parallel_face];  // 一共 num_parallel_face个指针数组
        length = new int32_t[num_parallel_face * 2]; // 0~num_parallel_face-1 send num_parallel_face~2*num_parallel_face-1 recv
        request_s = new MPI_Request[num_parallel_face];
        request_r = new MPI_Request[num_parallel_face];
        status_s = new MPI_Status[num_parallel_face];
        status_r = new MPI_Status[num_parallel_face];
        //------------------------------------------------------------------
        //=============================================================================================================

        //===============================================================================
        //----------------------------------------------------------------------
        // 收集信息
        index = 0;
        double require_send_num_edge[num_parallel_face], require_recv_num_edge[num_parallel_face];
        for (int i = 0; i < nblock; i++)
        {
            num_parallel = (*grids_info)(i).parallel_bc.size();
            for (int j = 0; j < num_parallel; j++)
            {
                // 耦合面不传
                // 修改2026/01/14 16：04 耦合面也需要传值
                // if ((*grids_info)(i).parallel_bc[j].target_block_name != (*grids_info)(i).parallel_bc[j].this_block_name)
                //     continue;
                // 周期不传
                if (fmod((*grids_info)(i).parallel_bc[j].send_flag, 2) != 0)
                    continue;
                // 作为为corner的信息接受方，发送所需要接受的Edge数量
                int edge_num = (*grids_info)(i).corner[j].size();
                require_recv_num_edge[index] = (double)edge_num;
                index++;
            }
        }
        PARALLEL::mpi_barrier();
        //----------------------------------------------------------------------
        // 发送接受所有数据
        index = 0;
        for (int i = 0; i < nblock; i++)
        {
            num_parallel = (*grids_info)(i).parallel_bc.size();
            for (int j = 0; j < num_parallel; j++)
            {
                // 修改2026/01/14 16：04 耦合面也需要传值
                // if ((*grids_info)(i).parallel_bc[j].target_block_name != (*grids_info)(i).parallel_bc[j].this_block_name || fmod((*grids_info)(i).parallel_bc[j].send_flag, 2) != 0)
                if (fmod((*grids_info)(i).parallel_bc[j].send_flag, 2) != 0)
                    continue;
                PARALLEL::mpi_data_send((*grids_info)(i).parallel_bc[j].tar_myid, (*grids_info)(i).parallel_bc[j].send_flag, &(require_recv_num_edge[index]), 1, &request_s[index]);
                PARALLEL::mpi_data_recv((*grids_info)(i).parallel_bc[j].tar_myid, (*grids_info)(i).parallel_bc[j].rece_flag, &(require_send_num_edge[index]), 1, &request_r[index]);
                index++;
            }
        }
        //----------------------------------------------------------------------
        // 等待完成
        PARALLEL::mpi_wait(num_parallel_face, request_s, status_s);
        PARALLEL::mpi_wait(num_parallel_face, request_r, status_r);
        PARALLEL::mpi_barrier();
        //----------------------------------------------------------------------
        // 取出所有数据
        index = 0;
        for (int i = 0; i < nblock; i++)
        {
            num_parallel = (*grids_info)(i).parallel_bc.size();
            for (int j = 0; j < num_parallel; j++)
            {
                // 修改2026/01/14 16：04 耦合面也需要传值
                // if ((*grids_info)(i).parallel_bc[j].target_block_name != (*grids_info)(i).parallel_bc[j].this_block_name || fmod((*grids_info)(i).parallel_bc[j].send_flag, 2) != 0)
                if (fmod((*grids_info)(i).parallel_bc[j].send_flag, 2) != 0)
                    continue;
                double temp_num = require_send_num_edge[index];
                (*grids_info)(i).corner_send[j].resize((int)temp_num);
                index++;
            }
        }
        PARALLEL::mpi_barrier();
        //===============================================================================
        // 发送接收Edge的起始点和位移范围，一个Edge有6个数，按面来处理
        //----------------------------------------------------------------------
        // 收集信息
        index = 0;
        double **require_send_edge;
        double **require_recv_edge;
        require_send_edge = new double *[num_parallel_face]; // 一共 num_parallel_face个指针数组
        require_recv_edge = new double *[num_parallel_face]; // 一共 num_parallel_face个指针数组
        for (int i = 0; i < nblock; i++)
        {
            num_parallel = (*grids_info)(i).parallel_bc.size();
            for (int j = 0; j < num_parallel; j++)
            {
                // 耦合面不传,周期不传
                // 修改2026/01/14 16：04 耦合面也需要传值
                // if ((*grids_info)(i).parallel_bc[j].target_block_name != (*grids_info)(i).parallel_bc[j].this_block_name || fmod((*grids_info)(i).parallel_bc[j].send_flag, 2) != 0)
                if (fmod((*grids_info)(i).parallel_bc[j].send_flag, 2) != 0)
                {
                    continue;
                }
                // 作为为corner的信息接受方
                int edge_num = (*grids_info)(i).corner[j].size();
                require_recv_edge[index] = new double[edge_num * 6 + 1];
                // 作为为corner的信息发送方
                edge_num = (*grids_info)(i).corner_send[j].size();
                require_send_edge[index] = new double[edge_num * 6 + 1];
                //-----------------------------------------------------------------------
                // 从corner取出放入require_recv_edge
                // ###########################################################################
                Parallel_Boundary &para = (*grids_info)(i).parallel_bc[j];
                int dir[3], cycle[3], direction_face, sub[3], sup[3];
                for (int ii = 0; ii < 3; ii++)
                {
                    sub[ii] = abs(para.sub[ii]);
                    sup[ii] = abs(para.sup[ii]);
                }
                for (int ii = 0; ii < 3; ii++)
                    cycle[ii] = para.cycle[ii];
                // dir表示从sub到sup是增大还是减小，便于循环统一书写
                for (int ii = 0; ii < 3; ii++)
                    dir[ii] = (sub[ii] <= sup[ii] + cycle[ii]) ? 1 : -1;
                direction_face = abs(para.direction) - 1;
                for (int ii = 0; ii < (*grids_info)(i).corner[j].size(); ii++)
                {
                    Edge &e = (*grids_info)(i).corner[j][ii];
                    int direct_eside = abs(e.direction1) - 1;
                    int direct_edge = 3 - direct_eside - direction_face;
                    int esub[3], esup[3], range_sub[3], range_sup[3], mid_sub[3], mid_sup[3];
                    for (int jj = 0; jj < 3; jj++)
                    {
                        esub[jj] = abs(e.sub[jj]);
                        esup[jj] = abs(e.sup[jj]);
                    }
                    // edge方向
                    range_sub[direct_edge] = (esub[direct_edge] - sub[direct_edge]) * dir[direct_edge];
                    range_sup[direct_edge] = (esup[direct_edge] - sub[direct_edge]) * dir[direct_edge];
                    // face的side方向
                    range_sub[direct_eside] = (esub[direct_eside] + e.cycle1[direct_eside] - sub[direct_eside]) * dir[direct_eside];
                    range_sup[direct_eside] = (esub[direct_eside] + e.cycle1[direct_eside] * depth - sub[direct_eside]) * dir[direct_eside];
                    // face的norm方向
                    range_sub[direction_face] = (esub[direction_face] + cycle[direction_face] - sub[direction_face]) * dir[direction_face];
                    range_sup[direction_face] = (esub[direction_face] + cycle[direction_face] * depth - sub[direction_face]) * dir[direction_face];
                    // Transform
                    for (int jj = 0; jj < 3; jj++)
                    {
                        mid_sub[para.Transform[jj]] = range_sub[jj];
                        mid_sup[para.Transform[jj]] = range_sup[jj];
                    }
                    for (int jj = 0; jj < 3; jj++)
                    {
                        e.sub_mid[jj] = mid_sub[jj];
                        e.sup_mid[jj] = mid_sup[jj];
                    }
                    require_recv_edge[index][6 * ii + 0] = (mid_sub[0] + 0.0);
                    require_recv_edge[index][6 * ii + 1] = (mid_sup[0] + 0.0);
                    require_recv_edge[index][6 * ii + 2] = (mid_sub[1] + 0.0);
                    require_recv_edge[index][6 * ii + 3] = (mid_sup[1] + 0.0);
                    require_recv_edge[index][6 * ii + 4] = (mid_sub[2] + 0.0);
                    require_recv_edge[index][6 * ii + 5] = (mid_sup[2] + 0.0);

                    // if(grd->grids(0).par->GetInt("myid") == 1 && para.tar_myid == 0)
                    // {
                    // 	std::cout << esub[0] << "\t" << esup[0] << "\t" << esub[1] << "\t" << esup[1] << "\t" << esub[2] << "\t" << esup[2] << "\n";
                    // 	std::cout << sub[0] << "\t\t" << sub[1] << "\t\t" << sub[2] << "\n";
                    // 	std::cout << para.Transform[0] << "\t\t" << para.Transform[1] << "\t\t" << para.Transform[2] << "\n";
                    // 	std::cout << mid_sub[0] << "\t" << mid_sup[0] << "\t" << mid_sub[1] << "\t" << mid_sup[1] << "\t" << mid_sub[2] << "\t" << mid_sup[2]
                    // 	<< "\n\n----------------------------\n" << std::flush;
                    // }
                }
                // ###########################################################################
                //-----------------------------------------------------------------------
                index++;
            }
        }
        PARALLEL::mpi_barrier();
        //----------------------------------------------------------------------
        // 发送接受所有数据
        index = 0;
        for (int i = 0; i < nblock; i++)
        {
            num_parallel = (*grids_info)(i).parallel_bc.size();
            for (int j = 0; j < num_parallel; j++)
            {
                // 修改2026/01/14 16：04 耦合面也需要传值
                // if ((*grids_info)(i).parallel_bc[j].target_block_name != (*grids_info)(i).parallel_bc[j].this_block_name || fmod((*grids_info)(i).parallel_bc[j].send_flag, 2) != 0)
                if (fmod((*grids_info)(i).parallel_bc[j].send_flag, 2) != 0)
                    continue;

                int edge_num = (*grids_info)(i).corner[j].size();
                PARALLEL::mpi_data_send((*grids_info)(i).parallel_bc[j].tar_myid, (*grids_info)(i).parallel_bc[j].send_flag, require_recv_edge[index], 6 * edge_num + 1, &request_s[index]);
                edge_num = (*grids_info)(i).corner_send[j].size();
                PARALLEL::mpi_data_recv((*grids_info)(i).parallel_bc[j].tar_myid, (*grids_info)(i).parallel_bc[j].rece_flag, require_send_edge[index], 6 * edge_num + 1, &request_r[index]);
                index++;
            }
        }
        //----------------------------------------------------------------------
        // 等待完成
        PARALLEL::mpi_wait(num_parallel_face, request_s, status_s);
        PARALLEL::mpi_wait(num_parallel_face, request_r, status_r);
        PARALLEL::mpi_barrier();
        //----------------------------------------------------------------------
        // 取出所有数据
        index = 0;
        for (int i = 0; i < nblock; i++)
        {
            num_parallel = (*grids_info)(i).parallel_bc.size();
            for (int j = 0; j < num_parallel; j++)
            {
                // 修改2026/01/14 16：04 耦合面也需要传值
                // if ((*grids_info)(i).parallel_bc[j].target_block_name != (*grids_info)(i).parallel_bc[j].this_block_name || fmod((*grids_info)(i).parallel_bc[j].send_flag, 2) != 0)
                if (fmod((*grids_info)(i).parallel_bc[j].send_flag, 2) != 0)
                    continue;
                std::vector<Edge> &e = (*grids_info)(i).corner_send[j];
                for (int ii = 0; ii < e.size(); ii++)
                {
                    e[ii].sub_mid[0] = require_send_edge[index][6 * ii + 0];
                    e[ii].sup_mid[0] = require_send_edge[index][6 * ii + 1];
                    e[ii].sub_mid[1] = require_send_edge[index][6 * ii + 2];
                    e[ii].sup_mid[1] = require_send_edge[index][6 * ii + 3];
                    e[ii].sub_mid[2] = require_send_edge[index][6 * ii + 4];
                    e[ii].sup_mid[2] = require_send_edge[index][6 * ii + 5];
                }
                index++;
            }
        }
        for (int index = 0; index < num_parallel_face; index++)
        {
            // 释放空间
            delete[] require_send_edge[index];
            delete[] require_recv_edge[index];
        }
        delete[] require_send_edge;
        delete[] require_recv_edge;
        PARALLEL::mpi_barrier();
        //===============================================================================

        //=============================================================================================================
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
                // 修改2026/01/14 16：04 耦合面也需要传值
                // if ((*grids_info)(i).parallel_bc[j].this_block_name != (*grids_info)(i).parallel_bc[j].target_block_name)
                //     continue;

                int direct_face = abs((*grids_info)(i).parallel_bc[j].direction) - 1, send_length = 1, recv_length = 1;
                // send
                for (int k = 0; k < (*grids_info)(i).corner_send[j].size(); k++)
                {
                    Edge &e = (*grids_info)(i).corner_send[j][k];
                    send_length += (abs(abs(e.sub_mid[0]) - abs(e.sup_mid[0])) + 1) * (abs(abs(e.sub_mid[1]) - abs(e.sup_mid[1])) + 1) * (abs(abs(e.sub_mid[2]) - abs(e.sup_mid[2])) + 1) * var_num;
                }
                // recv
                for (int k = 0; k < (*grids_info)(i).corner[j].size(); k++)
                {
                    Edge &e = (*grids_info)(i).corner[j][k];
                    int dir = 3 - (abs(e.direction1) - 1) - direct_face;

                    recv_length += (abs(abs(e.sub[dir]) - abs(e.sup[dir])) + 1) * depth * depth * var_num;
                }
                // 存储传输数组的长度
                length[index] = send_length;
                length[index + num_parallel_face] = recv_length;
                buf_send[index] = new double[length[index]];
                buf_recv[index] = new double[length[index + num_parallel_face]];
                index++;
            }
        }
        //------------------------------------------------------------------
    };

    void Parallel_corner_send_scalar(Parallel_Boundary *bound, std::vector<Edge> &bound_corner_send, int depth, double3D &scalar, double *buf_send, int32_t dimension)
    {
        int32_t inver_Transform[3], sub[3], sup[3], dir[3], cycle[3];
        int32_t ijk[3];
        // 面的范围
        for (int i = 0; i < 3; i++)
        {
            sub[i] = abs(bound->sub[i]);
            sup[i] = abs(bound->sup[i]);
        }
        //-------------------------------------------------------------------------
        // cycle是用来表征哪一个方向为法向
        // 这里的正负号是为了保证取数据的时候选的是内部的点，direction>0大号面，取负从而获得内部点
        for (int i = 0; i < 3; i++)
            cycle[i] = -bound->cycle[i];
        // dir表示从sub到sup是增大还是减小，便于循环统一书写
        for (int i = 0; i < 3; i++)
            dir[i] = (sub[i] <= sup[i] + cycle[i]) ? 1 : -1;
        //-------------------------------------------------------------------------
        // 求出inver_Transform
        for (int i = 0; i < 3; i++)
            inver_Transform[bound->Transform[i]] = i;
        // 以1,法向方向为最外层循环，从而能够使得边界与虚网格处理分开,默认0为中间层，2为内层循环
        int32_t out, mid, inner;
        out = inver_Transform[1];
        mid = inver_Transform[0];
        inner = inver_Transform[2];

        // Check depth
        if (bound_corner_send.size() != 0)
        {
            Edge &e = bound_corner_send[0];
            int32_t sub_[3], sup_[3], sub_mid[3], sup_mid[3];
            for (int i = 0; i < 3; i++)
            {
                sub_mid[i] = e.sub_mid[bound->Transform[i]];
                sup_mid[i] = e.sup_mid[bound->Transform[i]];
            }
            for (int i = 0; i < 3; i++)
            {
                sub_[i] = sub[i] + sub_mid[i] * dir[i];
                sup_[i] = sub[i] + sup_mid[i] * dir[i];
            }
            if (abs(sub_[abs(bound->direction) - 1] - sup_[abs(bound->direction) - 1]) + 1 != depth)
            {
                std::cout << "#Fatal Error: PARALLEL Corner Send 'depth' is not match with corner_send ! !\n";
                std::cout << "depth is " << depth << ", the corresponding corner_send size is " << abs(sub_[abs(bound->direction) - 1] - sup_[abs(bound->direction) - 1]) + 1 << std::endl;
                exit(-1);
            }
        }

        // 循环该并行通信面将所有的需要发送的棱边信息提取出来放入buf_send
        int32_t index = 0;
        for (int i_edge = 0; i_edge < bound_corner_send.size(); i_edge++)
        {
            Edge &e = bound_corner_send[i_edge];
            int32_t sub_[3], sup_[3], dir_[3], sub_mid[3], sup_mid[3];

            for (int i = 0; i < 3; i++)
            {
                sub_mid[i] = e.sub_mid[bound->Transform[i]];
                sup_mid[i] = e.sup_mid[bound->Transform[i]];
            }

            for (int i = 0; i < 3; i++)
            {
                sub_[i] = sub[i] + sub_mid[i] * dir[i];
                sup_[i] = sub[i] + sup_mid[i] * dir[i];
                dir_[i] = (sub_[i] <= sup_[i]) ? 1 : -1;
            }

            for (ijk[out] = sub_[out]; (ijk[out] - sup_[out]) * dir_[out] <= 0; ijk[out] += dir_[out])
                for (ijk[mid] = sub_[mid]; (ijk[mid] - sup_[mid]) * dir_[mid] <= 0; ijk[mid] += dir_[mid])
                    for (ijk[inner] = sub_[inner]; (ijk[inner] - sup_[inner]) * dir_[inner] <= 0; ijk[inner] += dir_[inner])
                    {
                        buf_send[index] = scalar(ijk[0], ijk[1], ijk[2]);
                        index++;
                    }
        }
    }

    void Parallel_corner_recv_scalar(Parallel_Boundary *bound, std::vector<Edge> &bound_corner, int depth, double3D &scalar, double *buf_recv, int32_t dimension)
    {
        int32_t inver_Transform[3], sub[3], sup[3], dir[3], cycle[3];
        int32_t ijk[3];
        // 面的范围
        for (int i = 0; i < 3; i++)
        {
            sub[i] = abs(bound->sub[i]);
            sup[i] = abs(bound->sup[i]);
        }
        //-------------------------------------------------------------------------
        // cycle是用来表征哪一个方向为法向
        // 这里的正负号是为了保证取数据的时候选的是内部的点，direction>0大号面，取负从而获得内部点
        for (int i = 0; i < 3; i++)
            cycle[i] = bound->cycle[i];
        // dir表示从sub到sup是增大还是减小，便于循环统一书写
        for (int i = 0; i < 3; i++)
            dir[i] = (sub[i] <= sup[i] + cycle[i]) ? 1 : -1;
        //-------------------------------------------------------------------------
        // 求出inver_Transform
        for (int i = 0; i < 3; i++)
            inver_Transform[bound->Transform[i]] = i;
        // 以1,法向方向为最外层循环，从而能够使得边界与虚网格处理分开,默认0为中间层，2为内层循环
        int32_t out, mid, inner;
        out = inver_Transform[1];
        mid = inver_Transform[0];
        inner = inver_Transform[2];

        // Check depth
        if (bound_corner.size() != 0)
        {
            Edge &e = bound_corner[0];
            int32_t sub_[3], sup_[3], sub_mid[3], sup_mid[3];
            for (int i = 0; i < 3; i++)
            {
                sub_mid[i] = e.sub_mid[bound->Transform[i]];
                sup_mid[i] = e.sup_mid[bound->Transform[i]];
            }
            for (int i = 0; i < 3; i++)
            {
                sub_[i] = sub[i] + sub_mid[i] * dir[i];
                sup_[i] = sub[i] + sup_mid[i] * dir[i];
            }
            if (abs(sub_[abs(bound->direction) - 1] - sup_[abs(bound->direction) - 1]) + 1 != depth)
            {
                std::cout << "#Fatal Error: PARALLEL Corner Recv 'depth' is not match with corner ! !\n";
                std::cout << "depth is " << depth << ", the corresponding corner_recv size is " << abs(sub_[abs(bound->direction) - 1] - sup_[abs(bound->direction) - 1]) + 1 << std::endl;
                exit(-1);
            }
        }

        // 循环该并行通信面将所有的需要发送的棱边信息提取出来放入buf_send
        int32_t index = 0;
        for (int i_edge = 0; i_edge < bound_corner.size(); i_edge++)
        {
            Edge &e = bound_corner[i_edge];
            int32_t sub_[3], sup_[3], dir_[3], sub_mid[3], sup_mid[3];

            for (int i = 0; i < 3; i++)
            {
                sub_mid[i] = e.sub_mid[bound->Transform[i]];
                sup_mid[i] = e.sup_mid[bound->Transform[i]];
            }

            for (int i = 0; i < 3; i++)
            {
                sub_[i] = sub[i] + sub_mid[i] * dir[i];
                sup_[i] = sub[i] + sup_mid[i] * dir[i];
                dir_[i] = (sub_[i] <= sup_[i]) ? 1 : -1;
            }

            for (ijk[out] = sub_[out]; (ijk[out] - sup_[out]) * dir_[out] <= 0; ijk[out] += dir_[out])
                for (ijk[mid] = sub_[mid]; (ijk[mid] - sup_[mid]) * dir_[mid] <= 0; ijk[mid] += dir_[mid])
                    for (ijk[inner] = sub_[inner]; (ijk[inner] - sup_[inner]) * dir_[inner] <= 0; ijk[inner] += dir_[inner])
                    {
                        scalar(ijk[0], ijk[1], ijk[2]) = buf_recv[index];
                        index++;
                    }
        }
    }

    /**
     * @brief 分配用于所有界面网格并行通信传递的缓存空间，需要网格指针
     */
    void Para_data_transfer::Initial_Allocate_Corner3D_Grid(Grid *grd, int32_t var_num)
    {
        //------------------------------------------------------------------
        // 预处理
        grids_info = &(grd->grids);
        depth = grd->ngg + 1;

        if (whether_allocate)
        {
            std::cout << "PARALLEL::Para_data_transfer Has Been Allocated! ! !\n";
            exit(-1);
        }
        whether_allocate = true; // 分配了空间，后期需要释放
        //------------------------------------------------------------------

        //------------------------------------------------------------------
        // 创建临时变量
        int32_t nblock = grd->nblock;
        int32_t index = 0;
        int32_t num_parallel; // 临时变量，表示每个网格块的并行界面个数
        //------------------------------------------------------------------

        //------------------------------------------------------------------
        // 设置空间，清空对应para_bound的corner corner_send
        for (int i = 0; i < nblock; i++)
        {
            num_parallel = (*grids_info)(i).parallel_bc.size();
            (*grids_info)(i).corner_send.resize(num_parallel);
            (*grids_info)(i).corner.resize(num_parallel);

            for (int j = 0; j < num_parallel; j++)
            {
                // 周期边界条件send recv flag均用奇数表示
                // 同一物理场才传
                // 修改2026/01/14 16：04 耦合面也需要传值
                // if (fmod((*grids_info)(i).parallel_bc[j].send_flag, 2) == 0 && (*grids_info)(i).parallel_bc[j].this_block_name == (*grids_info)(i).parallel_bc[j].target_block_name)
                if (fmod((*grids_info)(i).parallel_bc[j].send_flag, 2) == 0)
                {
                    (*grids_info)(i).corner_send[j].clear();
                    (*grids_info)(i).corner[j].clear();
                }
            }
        }
        //------------------------------------------------------------------
        // 寻找para对应的index
        auto find_para_index = [&](Block *blk, Parallel_Boundary *p) -> int
        {
            std::vector<Parallel_Boundary> &para = blk->parallel_bc;
            for (int i = 0; i < para.size(); i++)
            {
                if (para[i].sub[0] == p->sub[0] && para[i].sub[1] == p->sub[1] && para[i].sub[2] == p->sub[2] && para[i].sup[0] == p->sup[0] && para[i].sup[1] == p->sup[1] && para[i].sup[2] == p->sup[2])
                    return i;
            }
            std::cout << "Fata error, can not find the para_boundary ! ! !\n";
            return -1;
        };
        for (int i = 0; i < nblock; i++)
        {
            for (Edge &e : (*grids_info)(i).edge)
            {
                Edge e_temp(e);
                if (e.index != 1)
                    continue;
                // 修改2026/01/14 16：04 耦合面也需要传值
                // if (e.para_bound->target_block_name != e.para_bound->this_block_name)
                //     continue;
                int j = find_para_index(&(*grids_info)(i), e.para_bound);
                (*grids_info)(i).corner[j].push_back(e_temp);
            }
        }

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
                // 周期边界条件send recv flag均用奇数表示
                // 同一物理场才传
                // 修改2026/01/14 16：04 耦合面也需要传值
                // if (fmod((*grids_info)(i).parallel_bc[j].send_flag, 2) == 0 && (*grids_info)(i).parallel_bc[j].this_block_name == (*grids_info)(i).parallel_bc[j].target_block_name)
                if (fmod((*grids_info)(i).parallel_bc[j].send_flag, 2) == 0)
                    num_parallel_face++;
            }
        }
        //------------------------------
        //------------------------------------------------------------------
        // 开缓冲数组空间1
        buf_send = new double *[num_parallel_face];  // 一共 num_parallel_face个指针数组
        buf_recv = new double *[num_parallel_face];  // 一共 num_parallel_face个指针数组
        length = new int32_t[num_parallel_face * 2]; // 0~num_parallel_face-1 send num_parallel_face~2*num_parallel_face-1 recv
        request_s = new MPI_Request[num_parallel_face];
        request_r = new MPI_Request[num_parallel_face];
        status_s = new MPI_Status[num_parallel_face];
        status_r = new MPI_Status[num_parallel_face];
        //------------------------------------------------------------------
        //=============================================================================================================

        //===============================================================================
        //----------------------------------------------------------------------
        // 收集信息
        index = 0;
        double require_send_num_edge[num_parallel_face], require_recv_num_edge[num_parallel_face];
        for (int i = 0; i < nblock; i++)
        {
            num_parallel = (*grids_info)(i).parallel_bc.size();
            for (int j = 0; j < num_parallel; j++)
            {
                // 耦合面不传
                // 修改2026/01/14 16：04 耦合面也需要传值
                // if ((*grids_info)(i).parallel_bc[j].target_block_name != (*grids_info)(i).parallel_bc[j].this_block_name)
                //     continue;
                // 周期不传
                if (fmod((*grids_info)(i).parallel_bc[j].send_flag, 2) != 0)
                    continue;
                // 作为为corner的信息接受方，发送所需要接受的Edge数量
                int edge_num = (*grids_info)(i).corner[j].size();
                require_recv_num_edge[index] = (double)edge_num;
                index++;
            }
        }
        PARALLEL::mpi_barrier();
        //----------------------------------------------------------------------
        // 发送接受所有数据
        index = 0;
        for (int i = 0; i < nblock; i++)
        {
            num_parallel = (*grids_info)(i).parallel_bc.size();
            for (int j = 0; j < num_parallel; j++)
            {
                // 修改2026/01/14 16：04 耦合面也需要传值
                // if ((*grids_info)(i).parallel_bc[j].target_block_name != (*grids_info)(i).parallel_bc[j].this_block_name || fmod((*grids_info)(i).parallel_bc[j].send_flag, 2) != 0)
                if (fmod((*grids_info)(i).parallel_bc[j].send_flag, 2) != 0)
                    continue;
                PARALLEL::mpi_data_send((*grids_info)(i).parallel_bc[j].tar_myid, (*grids_info)(i).parallel_bc[j].send_flag, &(require_recv_num_edge[index]), 1, &request_s[index]);
                PARALLEL::mpi_data_recv((*grids_info)(i).parallel_bc[j].tar_myid, (*grids_info)(i).parallel_bc[j].rece_flag, &(require_send_num_edge[index]), 1, &request_r[index]);
                index++;
            }
        }
        //----------------------------------------------------------------------
        // 等待完成
        PARALLEL::mpi_wait(num_parallel_face, request_s, status_s);
        PARALLEL::mpi_wait(num_parallel_face, request_r, status_r);
        PARALLEL::mpi_barrier();
        //----------------------------------------------------------------------
        // 取出所有数据
        index = 0;
        for (int i = 0; i < nblock; i++)
        {
            num_parallel = (*grids_info)(i).parallel_bc.size();
            for (int j = 0; j < num_parallel; j++)
            {
                // 修改2026/01/14 16：04 耦合面也需要传值
                // if ((*grids_info)(i).parallel_bc[j].target_block_name != (*grids_info)(i).parallel_bc[j].this_block_name || fmod((*grids_info)(i).parallel_bc[j].send_flag, 2) != 0)
                if (fmod((*grids_info)(i).parallel_bc[j].send_flag, 2) != 0)
                    continue;
                double temp_num = require_send_num_edge[index];
                (*grids_info)(i).corner_send[j].resize((int)temp_num);
                index++;
            }
        }
        PARALLEL::mpi_barrier();
        //===============================================================================
        // 发送接收Edge的起始点和位移范围，一个Edge有6个数，按面来处理
        //----------------------------------------------------------------------
        // 收集信息
        index = 0;
        double **require_send_edge;
        double **require_recv_edge;
        require_send_edge = new double *[num_parallel_face]; // 一共 num_parallel_face个指针数组
        require_recv_edge = new double *[num_parallel_face]; // 一共 num_parallel_face个指针数组
        for (int i = 0; i < nblock; i++)
        {
            num_parallel = (*grids_info)(i).parallel_bc.size();
            for (int j = 0; j < num_parallel; j++)
            {
                // 耦合面不传,周期不传
                // 修改2026/01/14 16：04 耦合面也需要传值
                // if ((*grids_info)(i).parallel_bc[j].target_block_name != (*grids_info)(i).parallel_bc[j].this_block_name || fmod((*grids_info)(i).parallel_bc[j].send_flag, 2) != 0)
                if (fmod((*grids_info)(i).parallel_bc[j].send_flag, 2) != 0)
                {
                    continue;
                }
                // 作为为corner的信息接受方
                int edge_num = (*grids_info)(i).corner[j].size();
                require_recv_edge[index] = new double[edge_num * 6 + 1];
                // 作为为corner的信息发送方
                edge_num = (*grids_info)(i).corner_send[j].size();
                require_send_edge[index] = new double[edge_num * 6 + 1];
                //-----------------------------------------------------------------------
                // 从corner取出放入require_recv_edge
                // ###########################################################################
                Parallel_Boundary &para = (*grids_info)(i).parallel_bc[j];
                int dir[3], cycle[3], direction_face, sub[3], sup[3];
                for (int ii = 0; ii < 3; ii++)
                {
                    sub[ii] = abs(para.sub[ii]);
                    sup[ii] = abs(para.sup[ii]);
                }
                for (int ii = 0; ii < 3; ii++)
                    cycle[ii] = para.cycle[ii];
                // dir表示从sub到sup是增大还是减小，便于循环统一书写
                for (int ii = 0; ii < 3; ii++)
                    dir[ii] = (sub[ii] <= sup[ii] + cycle[ii]) ? 1 : -1;
                direction_face = abs(para.direction) - 1;
                for (int ii = 0; ii < (*grids_info)(i).corner[j].size(); ii++)
                {
                    Edge &e = (*grids_info)(i).corner[j][ii];
                    int direct_eside = abs(e.direction1) - 1;
                    int direct_edge = 3 - direct_eside - direction_face;
                    int esub[3], esup[3], range_sub[3], range_sup[3], mid_sub[3], mid_sup[3];
                    for (int jj = 0; jj < 3; jj++)
                    {
                        esub[jj] = abs(e.sub[jj]);
                        esup[jj] = abs(e.sup[jj]);
                    }
                    // edge方向
                    range_sub[direct_edge] = (esub[direct_edge] - sub[direct_edge]) * dir[direct_edge];
                    range_sup[direct_edge] = (esup[direct_edge] - sub[direct_edge]) * dir[direct_edge];
                    // face的side方向
                    range_sub[direct_eside] = (esub[direct_eside] + e.cycle1[direct_eside] - sub[direct_eside]) * dir[direct_eside];
                    range_sup[direct_eside] = (esub[direct_eside] + e.cycle1[direct_eside] * depth - sub[direct_eside]) * dir[direct_eside];
                    // face的norm方向
                    range_sub[direction_face] = (esub[direction_face] + cycle[direction_face] - sub[direction_face]) * dir[direction_face];
                    range_sup[direction_face] = (esub[direction_face] + cycle[direction_face] * depth - sub[direction_face]) * dir[direction_face];
                    // 在第三方向 direct_edge 上再延长 depth 层，形成角点立方体
                    if (depth > 0)
                    {
                        if (dir[direct_edge] > 0)
                            range_sup[direct_edge] += depth; // 向正方向延长
                        else
                            range_sub[direct_edge] -= depth; // 向负方向延长
                    }
                    // Transform
                    for (int jj = 0; jj < 3; jj++)
                    {
                        mid_sub[para.Transform[jj]] = range_sub[jj];
                        mid_sup[para.Transform[jj]] = range_sup[jj];
                    }
                    for (int jj = 0; jj < 3; jj++)
                    {
                        e.sub_mid[jj] = mid_sub[jj];
                        e.sup_mid[jj] = mid_sup[jj];
                    }
                    require_recv_edge[index][6 * ii + 0] = (mid_sub[0] + 0.0);
                    require_recv_edge[index][6 * ii + 1] = (mid_sup[0] + 0.0);
                    require_recv_edge[index][6 * ii + 2] = (mid_sub[1] + 0.0);
                    require_recv_edge[index][6 * ii + 3] = (mid_sup[1] + 0.0);
                    require_recv_edge[index][6 * ii + 4] = (mid_sub[2] + 0.0);
                    require_recv_edge[index][6 * ii + 5] = (mid_sup[2] + 0.0);

                    // if(grd->grids(0).par->GetInt("myid") == 1 && para.tar_myid == 0)
                    // {
                    // 	std::cout << esub[0] << "\t" << esup[0] << "\t" << esub[1] << "\t" << esup[1] << "\t" << esub[2] << "\t" << esup[2] << "\n";
                    // 	std::cout << sub[0] << "\t\t" << sub[1] << "\t\t" << sub[2] << "\n";
                    // 	std::cout << para.Transform[0] << "\t\t" << para.Transform[1] << "\t\t" << para.Transform[2] << "\n";
                    // 	std::cout << mid_sub[0] << "\t" << mid_sup[0] << "\t" << mid_sub[1] << "\t" << mid_sup[1] << "\t" << mid_sub[2] << "\t" << mid_sup[2]
                    // 	<< "\n\n----------------------------\n" << std::flush;
                    // }
                }
                // ###########################################################################
                //-----------------------------------------------------------------------
                index++;
            }
        }
        PARALLEL::mpi_barrier();
        //----------------------------------------------------------------------
        // 发送接受所有数据
        index = 0;
        for (int i = 0; i < nblock; i++)
        {
            num_parallel = (*grids_info)(i).parallel_bc.size();
            for (int j = 0; j < num_parallel; j++)
            {
                // 修改2026/01/14 16：04 耦合面也需要传值
                // if ((*grids_info)(i).parallel_bc[j].target_block_name != (*grids_info)(i).parallel_bc[j].this_block_name || fmod((*grids_info)(i).parallel_bc[j].send_flag, 2) != 0)
                if (fmod((*grids_info)(i).parallel_bc[j].send_flag, 2) != 0)
                    continue;

                int edge_num = (*grids_info)(i).corner[j].size();
                PARALLEL::mpi_data_send((*grids_info)(i).parallel_bc[j].tar_myid, (*grids_info)(i).parallel_bc[j].send_flag, require_recv_edge[index], 6 * edge_num + 1, &request_s[index]);
                edge_num = (*grids_info)(i).corner_send[j].size();
                PARALLEL::mpi_data_recv((*grids_info)(i).parallel_bc[j].tar_myid, (*grids_info)(i).parallel_bc[j].rece_flag, require_send_edge[index], 6 * edge_num + 1, &request_r[index]);
                index++;
            }
        }
        //----------------------------------------------------------------------
        // 等待完成
        PARALLEL::mpi_wait(num_parallel_face, request_s, status_s);
        PARALLEL::mpi_wait(num_parallel_face, request_r, status_r);
        PARALLEL::mpi_barrier();
        //----------------------------------------------------------------------
        // 取出所有数据
        index = 0;
        for (int i = 0; i < nblock; i++)
        {
            num_parallel = (*grids_info)(i).parallel_bc.size();
            for (int j = 0; j < num_parallel; j++)
            {
                // 修改2026/01/14 16：04 耦合面也需要传值
                // if ((*grids_info)(i).parallel_bc[j].target_block_name != (*grids_info)(i).parallel_bc[j].this_block_name || fmod((*grids_info)(i).parallel_bc[j].send_flag, 2) != 0)
                if (fmod((*grids_info)(i).parallel_bc[j].send_flag, 2) != 0)
                    continue;
                std::vector<Edge> &e = (*grids_info)(i).corner_send[j];
                for (int ii = 0; ii < e.size(); ii++)
                {
                    e[ii].sub_mid[0] = require_send_edge[index][6 * ii + 0];
                    e[ii].sup_mid[0] = require_send_edge[index][6 * ii + 1];
                    e[ii].sub_mid[1] = require_send_edge[index][6 * ii + 2];
                    e[ii].sup_mid[1] = require_send_edge[index][6 * ii + 3];
                    e[ii].sub_mid[2] = require_send_edge[index][6 * ii + 4];
                    e[ii].sup_mid[2] = require_send_edge[index][6 * ii + 5];
                }
                index++;
            }
        }
        for (int index = 0; index < num_parallel_face; index++)
        {
            // 释放空间
            delete[] require_send_edge[index];
            delete[] require_recv_edge[index];
        }
        delete[] require_send_edge;
        delete[] require_recv_edge;
        PARALLEL::mpi_barrier();
        //===============================================================================

        //=============================================================================================================
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
                // 修改2026/01/14 16：04 耦合面也需要传值
                // if ((*grids_info)(i).parallel_bc[j].this_block_name != (*grids_info)(i).parallel_bc[j].target_block_name)
                //     continue;

                int direct_face = abs((*grids_info)(i).parallel_bc[j].direction) - 1, send_length = 1, recv_length = 1;
                // send
                for (int k = 0; k < (*grids_info)(i).corner_send[j].size(); k++)
                {
                    Edge &e = (*grids_info)(i).corner_send[j][k];
                    // send_length += (abs(abs(e.sub_mid[0]) - abs(e.sup_mid[0])) + 1) * (abs(abs(e.sub_mid[1]) - abs(e.sup_mid[1])) + 1) * (abs(abs(e.sub_mid[2]) - abs(e.sup_mid[2])) + 1) * var_num;
                    send_length += (abs(e.sub_mid[0] - e.sup_mid[0]) + 1) * (abs(e.sub_mid[1] - e.sup_mid[1]) + 1) * (abs(e.sub_mid[2] - e.sup_mid[2]) + 1) * var_num;
                }
                // recv
                for (int k = 0; k < (*grids_info)(i).corner[j].size(); k++)
                {
                    Edge &e = (*grids_info)(i).corner[j][k];
                    int dir = 3 - (abs(e.direction1) - 1) - direct_face;

                    // recv_length += (abs(abs(e.sub[dir]) - abs(e.sup[dir])) + 1) * depth * depth * var_num;
                    // Corner3D +depth
                    recv_length += (abs(e.sub[dir] - e.sup[dir]) + 1 + depth) * depth * depth * var_num;
                }
                // 存储传输数组的长度
                length[index] = send_length;
                length[index + num_parallel_face] = recv_length;
                buf_send[index] = new double[length[index]];
                buf_recv[index] = new double[length[index + num_parallel_face]];
                index++;
            }
        }
        //------------------------------------------------------------------
    };

    void Parallel_corner3D_send_scalar(Parallel_Boundary *bound, std::vector<Edge> &bound_corner_send, int depth, double3D &scalar, double *buf_send, int32_t dimension)
    {
        int32_t inver_Transform[3], sub[3], sup[3], dir[3], cycle[3];
        int32_t ijk[3];
        // 面的范围
        for (int i = 0; i < 3; i++)
        {
            sub[i] = abs(bound->sub[i]);
            sup[i] = abs(bound->sup[i]);
        }
        //-------------------------------------------------------------------------
        // cycle是用来表征哪一个方向为法向
        // 这里的正负号是为了保证取数据的时候选的是内部的点，direction>0大号面，取负从而获得内部点
        for (int i = 0; i < 3; i++)
            cycle[i] = -bound->cycle[i];
        // dir表示从sub到sup是增大还是减小，便于循环统一书写
        for (int i = 0; i < 3; i++)
            dir[i] = (sub[i] <= sup[i] + cycle[i]) ? 1 : -1;
        //-------------------------------------------------------------------------
        // 求出inver_Transform
        for (int i = 0; i < 3; i++)
            inver_Transform[bound->Transform[i]] = i;
        // 以1,法向方向为最外层循环，从而能够使得边界与虚网格处理分开,默认0为中间层，2为内层循环
        int32_t out, mid, inner;
        out = inver_Transform[1];
        mid = inver_Transform[0];
        inner = inver_Transform[2];

        // Check depth
        if (bound_corner_send.size() != 0)
        {
            Edge &e = bound_corner_send[0];
            int32_t sub_[3], sup_[3], sub_mid[3], sup_mid[3];
            for (int i = 0; i < 3; i++)
            {
                sub_mid[i] = e.sub_mid[bound->Transform[i]];
                sup_mid[i] = e.sup_mid[bound->Transform[i]];
            }
            for (int i = 0; i < 3; i++)
            {
                sub_[i] = sub[i] + sub_mid[i] * dir[i];
                sup_[i] = sub[i] + sup_mid[i] * dir[i];
            }
            if (abs(sub_[abs(bound->direction) - 1] - sup_[abs(bound->direction) - 1]) + 1 != depth)
            {
                std::cout << "#Fatal Error: PARALLEL Corner Send 'depth' is not match with corner_send ! !\n";
                std::cout << "depth is " << depth << ", the corresponding corner_send size is " << abs(sub_[abs(bound->direction) - 1] - sup_[abs(bound->direction) - 1]) + 1 << std::endl;
                exit(-1);
            }
        }

        // 循环该并行通信面将所有的需要发送的棱边信息提取出来放入buf_send
        int32_t index = 0;
        for (int i_edge = 0; i_edge < bound_corner_send.size(); i_edge++)
        {
            Edge &e = bound_corner_send[i_edge];
            int32_t sub_[3], sup_[3], dir_[3], sub_mid[3], sup_mid[3];

            for (int i = 0; i < 3; i++)
            {
                sub_mid[i] = e.sub_mid[bound->Transform[i]];
                sup_mid[i] = e.sup_mid[bound->Transform[i]];
            }

            for (int i = 0; i < 3; i++)
            {
                sub_[i] = sub[i] + sub_mid[i] * dir[i];
                sup_[i] = sub[i] + sup_mid[i] * dir[i];
                dir_[i] = (sub_[i] <= sup_[i]) ? 1 : -1;
            }

            for (ijk[out] = sub_[out]; (ijk[out] - sup_[out]) * dir_[out] <= 0; ijk[out] += dir_[out])
                for (ijk[mid] = sub_[mid]; (ijk[mid] - sup_[mid]) * dir_[mid] <= 0; ijk[mid] += dir_[mid])
                    for (ijk[inner] = sub_[inner]; (ijk[inner] - sup_[inner]) * dir_[inner] <= 0; ijk[inner] += dir_[inner])
                    {
                        buf_send[index] = scalar(ijk[0], ijk[1], ijk[2]);
                        index++;
                    }
        }
    }

    void Parallel_corner3D_recv_scalar(Parallel_Boundary *bound, std::vector<Edge> &bound_corner, int depth, double3D &scalar, double *buf_recv, int32_t dimension)
    {
        int32_t inver_Transform[3], sub[3], sup[3], dir[3], cycle[3];
        int32_t ijk[3];
        // 面的范围
        for (int i = 0; i < 3; i++)
        {
            sub[i] = abs(bound->sub[i]);
            sup[i] = abs(bound->sup[i]);
        }
        //-------------------------------------------------------------------------
        // cycle是用来表征哪一个方向为法向
        // 这里的正负号是为了保证取数据的时候选的是内部的点，direction>0大号面，取负从而获得内部点
        for (int i = 0; i < 3; i++)
            cycle[i] = bound->cycle[i];
        // dir表示从sub到sup是增大还是减小，便于循环统一书写
        for (int i = 0; i < 3; i++)
            dir[i] = (sub[i] <= sup[i] + cycle[i]) ? 1 : -1;
        //-------------------------------------------------------------------------
        // 求出inver_Transform
        for (int i = 0; i < 3; i++)
            inver_Transform[bound->Transform[i]] = i;
        // 以1,法向方向为最外层循环，从而能够使得边界与虚网格处理分开,默认0为中间层，2为内层循环
        int32_t out, mid, inner;
        out = inver_Transform[1];
        mid = inver_Transform[0];
        inner = inver_Transform[2];

        // Check depth
        if (bound_corner.size() != 0)
        {
            Edge &e = bound_corner[0];
            int32_t sub_[3], sup_[3], sub_mid[3], sup_mid[3];
            for (int i = 0; i < 3; i++)
            {
                sub_mid[i] = e.sub_mid[bound->Transform[i]];
                sup_mid[i] = e.sup_mid[bound->Transform[i]];
            }
            for (int i = 0; i < 3; i++)
            {
                sub_[i] = sub[i] + sub_mid[i] * dir[i];
                sup_[i] = sub[i] + sup_mid[i] * dir[i];
            }
            if (abs(sub_[abs(bound->direction) - 1] - sup_[abs(bound->direction) - 1]) + 1 != depth)
            {
                std::cout << "#Fatal Error: PARALLEL Corner Recv 'depth' is not match with corner ! !\n";
                std::cout << "depth is " << depth << ", the corresponding corner_recv size is " << abs(sub_[abs(bound->direction) - 1] - sup_[abs(bound->direction) - 1]) + 1 << std::endl;
                exit(-1);
            }
        }

        // 循环该并行通信面将所有的需要发送的棱边信息提取出来放入buf_send
        int32_t index = 0;
        for (int i_edge = 0; i_edge < bound_corner.size(); i_edge++)
        {
            Edge &e = bound_corner[i_edge];
            int32_t sub_[3], sup_[3], dir_[3], sub_mid[3], sup_mid[3];

            for (int i = 0; i < 3; i++)
            {
                sub_mid[i] = e.sub_mid[bound->Transform[i]];
                sup_mid[i] = e.sup_mid[bound->Transform[i]];
            }

            for (int i = 0; i < 3; i++)
            {
                sub_[i] = sub[i] + sub_mid[i] * dir[i];
                sup_[i] = sub[i] + sup_mid[i] * dir[i];
                dir_[i] = (sub_[i] <= sup_[i]) ? 1 : -1;
            }

            for (ijk[out] = sub_[out]; (ijk[out] - sup_[out]) * dir_[out] <= 0; ijk[out] += dir_[out])
                for (ijk[mid] = sub_[mid]; (ijk[mid] - sup_[mid]) * dir_[mid] <= 0; ijk[mid] += dir_[mid])
                    for (ijk[inner] = sub_[inner]; (ijk[inner] - sup_[inner]) * dir_[inner] <= 0; ijk[inner] += dir_[inner])
                    {
                        scalar(ijk[0], ijk[1], ijk[2]) = buf_recv[index];
                        index++;
                    }
        }
    }
}