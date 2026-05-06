#include "0_basic/DEFINE.h"
#include "1_grid/1_MPCNS_Grid.h"
#include "1_grid/Grid_Trans.h"
#include <fstream>

/**
 * @brief 读取grd、inp、fvbnd文件
 * @remark 本程序依赖ngg参数
 */
void Grid::ReadGeometry(int32_t my_id, int32_t ghostmesh, int32_t Method)
{
    ngg = ghostmesh;

    std::string _my_id_s;
    if (my_id < 10)
    {
        _my_id_s = "   " + std::to_string(my_id);
    }
    else if (my_id < 100)
    {
        _my_id_s = "  " + std::to_string(my_id);
    }
    else if (my_id < 1000)
    {
        _my_id_s = " " + std::to_string(my_id);
    }
    else // 这说明并行进程数不得超过9999
    {
        _my_id_s = std::to_string(my_id);
    }
    if (Method == 0)
        Read_Grid(_my_id_s, my_id);
    else
        Read_Grid_Binary(_my_id_s, my_id);
    Read_Phy_Boundary(_my_id_s);
    Read_Phy_Name(_my_id_s, my_id);
    Adjust_Phy_Boundary();
    Read_Inner_Boundary(_my_id_s, my_id);
    Read_Parallel_Boundary(_my_id_s, my_id);
}

void Grid::Read_Grid(std::string _my_id_s, int32_t my_id)
{
    /*****************读取网格******************/
    std::ifstream grd("./CASE/geometry/grid" + _my_id_s + ".grd", std::ios::in);
    std::string input;
    /****************读取网格块数*****************/
    int32_t numofblock = 0;
    grd >> numofblock;
    nblock = numofblock;
    grids.SetSize(numofblock);
    /***************读取网格名称******************/
    std::string phy_name;
    for (int iblock = 0; iblock < nblock; iblock++)
    {
        grd >> phy_name;
        grids(iblock).block_name = phy_name; // 每一块都需要对应一个名称，用来选择物理求解器
    }
    /***************读取网格点数****************/
    int32_t *mx = new int32_t[numofblock];
    int32_t *my = new int32_t[numofblock];
    int32_t *mz = new int32_t[numofblock];
    dimension = 3;
    for (int iblock = 0; iblock < nblock; iblock++)
        grd >> mx[iblock] >> my[iblock] >> mz[iblock];
    // 注意这里的mx my mz是临时的网格点数，且是从1开始的

    if (mz[0] == 1)
        dimension = 2;

    /*****************读取网格******************/
    /**************读取网格点坐标***************/
    for (int izone = 0; izone < nblock; izone++)
    {
        Block &iblock = grids(izone);

        iblock.imax = mx[izone];
        iblock.jmax = my[izone];
        iblock.kmax = mz[izone];
        iblock.mx = mx[izone] - 1;
        iblock.my = my[izone] - 1;
        iblock.mz = mz[izone] - 1;
        iblock.dimension = dimension;

        double3D &x = iblock.x;
        double3D &y = iblock.y;
        double3D &z = iblock.z;
        double3D &jacobi = iblock.jacobi;
        Phy_Tensor &metric = iblock.metric;
        Phy_Tensor &inverse_metric = iblock.inverse_metric;

        x.SetSize(mx[izone] + 2 * (ngg + 1), my[izone] + 2 * (ngg + 1), mz[izone] + 2 * (ngg + 1), ngg + 1);
        y.SetSize(mx[izone] + 2 * (ngg + 1), my[izone] + 2 * (ngg + 1), mz[izone] + 2 * (ngg + 1), ngg + 1);
        z.SetSize(mx[izone] + 2 * (ngg + 1), my[izone] + 2 * (ngg + 1), mz[izone] + 2 * (ngg + 1), ngg + 1);
        jacobi.SetSize(mx[izone] + 2 * ngg, my[izone] + 2 * ngg, mz[izone] + 2 * ngg, ngg);
        metric.SetSize(mx[izone] + 2 * ngg, my[izone] + 2 * ngg, mz[izone] + 2 * ngg, ngg, 3, 3);
        inverse_metric.SetSize(mx[izone] + 2 * ngg, my[izone] + 2 * ngg, mz[izone] + 2 * ngg, ngg, 3, 3);
        // 由于需要计算度量系数，因此需要额外添加一层虚网格，注意对于正常网格点100个的mx=100
        // 就坐标而言需要0~99号的正常网格和-ngg-1~0以及100~100+ngg虚网格共100+ngg+ngg+1+1个点
        // 一般表示式即为mx+2*(ngg+1)
        // 设置虚网格为ngg+1表示数组可以从-ngg-1开始索引

        /*读取grd坐标*/
        for (int k = 0; k < mz[izone]; ++k)
            for (int j = 0; j < my[izone]; ++j)
                for (int i = 0; i < mx[izone]; ++i)
                    grd >> x(i, j, k);
        for (int k = 0; k < mz[izone]; ++k)
            for (int j = 0; j < my[izone]; ++j)
                for (int i = 0; i < mx[izone]; ++i)
                    grd >> y(i, j, k);
        for (int k = 0; k < mz[izone]; ++k)
            for (int j = 0; j < my[izone]; ++j)
                for (int i = 0; i < mx[izone]; ++i)
                    grd >> z(i, j, k);
    }

    if (my_id == 0)
        std::cout << "\t--> Read grid successfully!" << std::endl;

    delete[] mx, my, mz;
    grd.close();
}

void Grid::Read_Grid_Binary(std::string _my_id_s, int32_t my_id)
{
    /*****************读取网格******************/
    std::ifstream grd("./CASE/geometry/grid" + _my_id_s + ".grd", std::ios::in | std::ios::binary);
    int32_t temp_int, doublelength = sizeof(double);
    double temp_double;
    /****************读取网格块数*****************/
    int32_t numofblock = 0;
    // grd >> numofblock;
    grd.read((char *)&temp_int, sizeof(int32_t));
    numofblock = temp_int;
    nblock = numofblock;
    grids.SetSize(numofblock);
    /***************读取网格名称******************/
    char *temp_char; // 临时存储string
    for (int32_t iblock = 0; iblock < nblock; iblock++)
    {
        grd.read((char *)&temp_int, sizeof(int32_t)); // 字符串长度

        temp_char = new char[temp_int];
        grd.read(temp_char, temp_int); // 读取字符串
        grids(iblock).block_name.reserve(temp_int);
        grids(iblock).block_name = temp_char; // 每一块都需要对应一个名称，用来选择物理求解器
        delete[] temp_char;
    }
    /***************读取网格点数****************/
    int32_t *mx = new int32_t[numofblock];
    int32_t *my = new int32_t[numofblock];
    int32_t *mz = new int32_t[numofblock];
    dimension = 3;
    for (int iblock = 0; iblock < nblock; iblock++)
    {
        grd.read((char *)&temp_int, sizeof(int32_t)); // imax
        mx[iblock] = temp_int;
        grd.read((char *)&temp_int, sizeof(int32_t)); // jmax
        my[iblock] = temp_int;
        grd.read((char *)&temp_int, sizeof(int32_t)); // kmax
        mz[iblock] = temp_int;
    }
    // 注意这里的mx my mz是临时的网格点数，且是从1开始的

    if (mz[0] == 1)
        dimension = 2;

    /*****************读取网格******************/
    /**************读取网格点坐标***************/
    for (int izone = 0; izone < nblock; izone++)
    {
        Block &iblock = grids(izone);

        iblock.imax = mx[izone];
        iblock.jmax = my[izone];
        iblock.kmax = mz[izone];
        iblock.mx = mx[izone] - 1;
        iblock.my = my[izone] - 1;
        iblock.mz = mz[izone] - 1;
        iblock.dimension = dimension;

        double3D &x = iblock.x;
        double3D &y = iblock.y;
        double3D &z = iblock.z;
        double3D &jacobi = iblock.jacobi;
        Phy_Tensor &metric = iblock.metric;
        Phy_Tensor &inverse_metric = iblock.inverse_metric;

        x.SetSize(mx[izone] + 2 * (ngg + 1), my[izone] + 2 * (ngg + 1), mz[izone] + 2 * (ngg + 1), ngg + 1);
        y.SetSize(mx[izone] + 2 * (ngg + 1), my[izone] + 2 * (ngg + 1), mz[izone] + 2 * (ngg + 1), ngg + 1);
        z.SetSize(mx[izone] + 2 * (ngg + 1), my[izone] + 2 * (ngg + 1), mz[izone] + 2 * (ngg + 1), ngg + 1);
        jacobi.SetSize(mx[izone] + 2 * ngg, my[izone] + 2 * ngg, mz[izone] + 2 * ngg, ngg);
        metric.SetSize(mx[izone] + 2 * ngg, my[izone] + 2 * ngg, mz[izone] + 2 * ngg, ngg, 3, 3);
        inverse_metric.SetSize(mx[izone] + 2 * ngg, my[izone] + 2 * ngg, mz[izone] + 2 * ngg, ngg, 3, 3);
        // 由于需要计算度量系数，因此需要额外添加一层虚网格，注意对于正常网格点100个的mx=100
        // 就坐标而言需要0~99号的正常网格和-ngg-1~0以及100~100+ngg虚网格共100+ngg+ngg+1+1个点
        // 一般表示式即为mx+2*(ngg+1)
        // 设置虚网格为ngg+1表示数组可以从-ngg-1开始索引

        /*读取grd坐标*/
        for (int k = 0; k < mz[izone]; ++k)
            for (int j = 0; j < my[izone]; ++j)
                for (int i = 0; i < mx[izone]; ++i)
                {
                    grd.read((char *)&temp_double, doublelength);
                    x(i, j, k) = temp_double;
                }
        for (int k = 0; k < mz[izone]; ++k)
            for (int j = 0; j < my[izone]; ++j)
                for (int i = 0; i < mx[izone]; ++i)
                {
                    grd.read((char *)&temp_double, doublelength);
                    y(i, j, k) = temp_double;
                }
        for (int k = 0; k < mz[izone]; ++k)
            for (int j = 0; j < my[izone]; ++j)
                for (int i = 0; i < mx[izone]; ++i)
                {
                    grd.read((char *)&temp_double, doublelength);
                    z(i, j, k) = temp_double;
                }
    }

    if (my_id == 0)
        std::cout << "\t--> Read grid successfully!" << std::endl;

    delete[] mx, my, mz;
    grd.close();
}

void Grid::Read_Phy_Boundary(std::string _my_id_s)
{
    /*****************边界条件******************/
    /***********读取boundary外边界条件**********/
    std::ifstream grd;
    grd.open("./CASE/geometry/boundary_condition/boundary" + _my_id_s + ".txt", std::ios::in);
    for (int izone = 0; izone < nblock; izone++)
    {
        Block &iblock = grids(izone);

        int nbcface;
        grd >> nbcface; // 读取外边界面个数

        for (int l = 0; l < nbcface; l++)
        {
            Physical_Boundary bcl;
            int32_t bctype;
            int32_t pointst[3], pointed[3], pointst_[3], pointed_[3];

            grd >> pointst_[0] >> pointed_[0] >> pointst_[1] >> pointed_[1] >> pointst_[2] >> pointed_[2];
            grd >> bctype;
            // Remark:边界条件文件是按照以0为起点

            for (int p = 0; p < 3; p++)
            {
                pointst[p] = (int)fmin(pointst_[p], pointed_[p]);
                pointed[p] = (int)fmax(pointst_[p], pointed_[p]);
                bcl.cycle[p] = 0;
            }

            for (int p = 0; p < 3; p++)
            {
                bcl.sub[p] = pointst[p];
                bcl.sup[p] = pointed[p];
                if (p == 2 && dimension == 2)
                    continue;

                if (bcl.sub[p] == bcl.sup[p])
                {
                    if (bcl.sub[p] == 0)
                    {
                        bcl.direction = -(p + 1);
                        bcl.cycle[p] = -1;
                    }
                    else
                    {
                        bcl.direction = (p + 1);
                        bcl.cycle[p] = 1;
                    }
                }
            }

            bcl.boundary_num = bctype;

            bcl.this_block_num = izone;

            bcl.this_block_name = iblock.block_name;

            iblock.physical_bc.push_back(bcl);
        }
    }
    grd.close();
}

void Grid::Read_Phy_Name(std::string _my_id_s, int32_t my_id)
{
    /***********读取fvbnd物理边界条件以更新bcname**********/
    std::ifstream grd;
    std::string input;
    grd.open("./CASE/geometry/boundary_condition/boundary_name.fvbnd", std::ios::in);
    std::getline(grd, input);                                              // FVBND 1 3
    std::getline(grd, input);                                              //(0)No_boundary_Conditions,注意没有为1的类型
    boundary_name.push_back(input);                                        //(0)No_boundary_Conditions
    boundary_name.push_back("");                                           //(1)""
    boundary_name_index.insert(std::pair<std::string, int32_t>(input, 0)); //(0)No_boundary_Conditions
    int bc_number = 1;
    while (true)
    {
        bc_number = bc_number + 1;
        grd >> input; // 读入名称
        // std::getline(grd, input);		//读入名称
        // input.resize(input.size() - 1); //末尾有换行
        if (input == "BOUNDARIES")
        {
            break;
        }
        boundary_name_index.insert(std::pair<std::string, int>(input, bc_number));
        boundary_name.push_back(input);
    }
    grd.close();
    // 读取fvbnd物理边界条件名称成功,上面boundary_name由号索引字符串，boundary_name_index由字符串索引号码
    for (int32_t izone = 0; izone < nblock; izone++)
    {
        Block &iblock = grids(izone);

        for (int32_t l = 0; l < iblock.physical_bc.size(); l++)
        {
            iblock.physical_bc[l].boundary_name = boundary_name[iblock.physical_bc[l].boundary_num];
        }
    }
    // 更新bcname完成

    if (my_id == 0)
        std::cout << "\t--> Read physical boundary successfully!" << std::endl;
}

void Grid::Adjust_Phy_Boundary()
{
    List<int> Priority = par->GetInt_List("Boundary_Priority");
    List<int> Priority_order;
    Priority_order.num = 0;
    int32_t Priority_size = Priority.num;

    std::map<std::string, int32_t>::iterator iter_memo;
    do
    {
        int index_max = -1;
        for (std::map<std::string, int32_t>::iterator iter = Priority.data.begin(); iter != Priority.data.end(); iter++)
        {
            if (iter->second >= index_max)
            {
                index_max = iter->second;
                iter_memo = iter;
            }
        }
        if (Priority.data.size() > 0)
        {
            Priority_order.data[iter_memo->first] = iter_memo->second;
            Priority_order.num++;
            Priority.data.erase(iter_memo);
        }
        else
            break;
    } while (true);

    //=============================================================================================
    // Adjust Boundary Priority
    // Proposed by Dch, Added by Descartes 2023-12-26
    for (int izone = 0; izone < nblock; izone++)
    {
        std::vector<Physical_Boundary> &iphysical_bc = grids(izone).physical_bc;
        if (iphysical_bc.size() <= 1)
            continue;
        // Adjust
        for (std::map<std::string, int32_t>::iterator iter = Priority_order.data.begin(); iter != Priority_order.data.end(); iter++)
        {
            std::string name = iter->first;
            bool if_break;
            for (int32_t i = iphysical_bc.size() - 1; i >= 0; i--)
            {
                if (iphysical_bc[i].boundary_name == name)
                {
                    Physical_Boundary temp = iphysical_bc[i];
                    iphysical_bc.push_back(temp);
                    std::vector<Physical_Boundary>::iterator index = iphysical_bc.begin();
                    iphysical_bc.erase(index + i);
                }
            }
        }
        //  //==========================================================================================================
        // //Test only
        // for(std::vector<Physical_Boundary>::iterator iter = iphysical_bc.begin(); iter != iphysical_bc.end(); iter++)
        //     std::cout<<iter->boundary_name<<"\t";
        // std::cout<<"\n";
        //  //==========================================================================================================
    }
    //=============================================================================================
}

void Grid::Read_Inner_Boundary(std::string _my_id_s, int32_t my_id)
{
    /*****************边界条件******************/
    /**********读取inner内边界对应关系**********/
    std::ifstream grd;
    std::string input;
    grd.open("./CASE/geometry/boundary_condition/inner" + _my_id_s + ".txt", std::ios::in);

    // 先读内边界面总数
    std::getline(grd, input);
    int num_inner_face = std::stoi(input);
    int *faces = new int[nblock];
    // 然后读各个block的内边界数
    for (int iblock = 0; iblock < nblock; iblock++)
        grd >> faces[iblock];
    std::getline(grd, input);
    std::getline(grd, input);

    // 读各个内边界的点范围以及接触面id
    for (int izone = 0; izone < nblock; izone++)
    {
        Block &iblock = grids(izone);

        int inner_iblock = faces[izone];

        for (int iface = 0; iface < inner_iblock; iface++)
        {
            Inner_Boundary innerbc;
            int target_block, source_block; // 内部传值目标块，本块的标号
            int pointst[3], pointed[3], targetst[3], targeted[3];
            std::string namest, nameed;

            // 注意：读取的边界点起始坐标从1开始，需要在后面转换成从0开始
            grd >> pointst[0] >> pointed[0] >> pointst[1] >> pointed[1] >> pointst[2] >> pointed[2];
            grd >> source_block;
            grd >> namest;

            grd >> targetst[0] >> targeted[0] >> targetst[1] >> targeted[1] >> targetst[2] >> targeted[2];
            grd >> target_block;
            grd >> nameed;

            for (int i = 0; i < 3; i++)
            {
                innerbc.sub[i] = pointst[i];
                innerbc.sup[i] = pointed[i];
                innerbc.tar_sub[i] = targetst[i];
                innerbc.tar_sup[i] = targeted[i];
            }
            innerbc.this_block_name = namest;
            innerbc.target_block_name = nameed;
            // 注意，这里的本块在进程中的编号是从1开始的，需要转化为从0开始故而为source_block-1
            ////// innerbc.this_block_num = source_block - 1;
            // Modify(2022/10/03/18:50)该数可正可负，不再直接减1，使用时注意取绝对值再减1
            // innerbc.this_block_num = source_block;
            // Modify(2024/04/26/17:13)该数可正可负，从1开始，读入后直接取绝对值再减1，并将是否为周期信息存入is_period
            innerbc.is_period = (source_block < 0);
            innerbc.this_block_num = abs(source_block) - 1;

            // 注意，这里的目标块在进程中的编号也是从1开始的，需要转化为从0开始故而为target_block-1
            ////// innerbc.tar_block_num = target_block - 1;
            // Modify(2022/10/03/18:50)该数可正可负，不再直接减1，使用时注意取绝对值再减1
            // innerbc.tar_block_num = target_block;
            // Modify(2024/04/26/17:13)该数可正可负，从1开始，读入后直接取绝对值再减1，并将是否为周期信息存入is_period
            innerbc.tar_block_num = abs(target_block) - 1;

            // 注意：读取的边界点起始坐标从1开始，换成从0开始
            innerbc.Pre_process(dimension);

            iblock.inner_bc.push_back(innerbc);
        }
    }
    if (my_id == 0)
        std::cout << "\t--> Read inner boundary successfully!" << std::endl;

    delete[] faces;
    grd.close();
    // 处理this_block_index tar_block_index主要在内部传值时使用
    for (int izone = 0; izone < nblock; izone++)
    {
        Block &iblock = grids(izone);

        int inner_iblock = iblock.inner_bc.size();

        for (int iface = 0; iface < inner_iblock; iface++)
        {
            // 若此inner面的index不是初始给的-1,表明已经处理，直接跳过
            if (iblock.inner_bc[iface].this_block_index != -1)
                continue;
            //-------------------------------------------------
            Inner_Boundary &this_inner = iblock.inner_bc[iface];
            // Modify(2024/04/26/17:13)从0开始且总为正
            Block &tar_block = grids(iblock.inner_bc[iface].tar_block_num);
            // Block &tar_block = grids(abs(iblock.inner_bc[iface].tar_block_num) - 1);

            // 若此inner面的index是初始的-1 ，修正
            this_inner.this_block_index = iface;

            // 获取目标块的总inner面数，依次循环搜索
            int tar_iblock_innernum = tar_block.inner_bc.size();
            for (int jface = 0; jface < tar_iblock_innernum; jface++)
            {
                Inner_Boundary &tar_inner = tar_block.inner_bc[jface];
                // 若此inner面的index不是初始给的-1,表明已经处理，直接跳过
                if (tar_inner.this_block_index != -1)
                    continue;
                if (tar_inner.sub[0] == this_inner.tar_sub[0] && tar_inner.sub[1] == this_inner.tar_sub[1] && tar_inner.sub[2] == this_inner.tar_sub[2] && tar_inner.sup[0] == this_inner.tar_sup[0] && tar_inner.sup[1] == this_inner.tar_sup[1] && tar_inner.sup[2] == this_inner.tar_sup[2])
                {
                    // 找到目标块的index，存入this_inner中，同时更新目标块的信息
                    this_inner.tar_block_index = jface;
                    tar_inner.tar_block_index = iface;
                    tar_inner.this_block_index = jface;
                    break;
                }
            }
        }
    }
}

void Grid::Read_Parallel_Boundary(std::string _my_id_s, int32_t my_id)
{
    /*****************边界条件******************/
    /********读取parallel 进程间对接关系********/
    std::ifstream grd;
    std::string input;
    grd.open("./CASE/geometry/boundary_condition/parallel" + _my_id_s + ".txt", std::ios::in);
    int blk_num1 = 0;
    grd >> blk_num1; // 首先读入该进程的块数，看看与网格坐标是否匹配
    if (blk_num1 != nblock)
    {
        std::cout << "Read parallel error! #block in parallel file and parallel file is different! The problem is in process " << my_id << std::endl;
    }

    int num_parallel_face;
    grd >> num_parallel_face; // 读入需要进程间通信的面数

    int *nface = new int[nblock];
    for (int izone = 0; izone < nblock; izone++)
        grd >> nface[izone]; // 获得每块网格的通信面数
    std::getline(grd, input);
    std::getline(grd, input); // 文字注释行

    for (int izone = 0; izone < nblock; izone++)
    {
        Block &iblock = grids(izone);
        int iface_num = nface[izone];

        for (int iiface = 0; iiface < iface_num; iiface++)
        {

            Parallel_Boundary iParallelBC;
            int pointst[3], pointed[3];
            int srid, sflag, rflag;
            std::string nameed;
            // 读取imin、imax、jmin、jmax、kmin、kmax, 读取s_r_id、s_flag、r_flag
            grd >> pointst[0] >> pointed[0] >> pointst[1] >> pointed[1] >> pointst[2] >> pointed[2] >> srid >> sflag >> rflag;
            grd >> nameed;
            // 与srid进程中的块通信，这里进程号从0开始故而不用改变
            iParallelBC.tar_myid = srid;
            iParallelBC.this_myid = my_id;
            iParallelBC.this_block_num = izone;

            for (int i = 0; i < 3; i++)
            {
                iParallelBC.sub[i] = pointst[i];
                iParallelBC.sup[i] = pointed[i];
            }

            iParallelBC.send_flag = sflag;
            iParallelBC.rece_flag = rflag;

            iParallelBC.this_block_name = iblock.block_name;
            iParallelBC.target_block_name = nameed;

            // 注意：读取的边界点起始坐标从1开始，在SetPoint中转换成从0开始
            iParallelBC.Pre_process(dimension);

            iblock.parallel_bc.push_back(iParallelBC);
        }
    }
    if (my_id == 0)
        std::cout << "\t--> Read parallel information successfully!\n";

    delete[] nface;
    grd.close();
}

/**
 * @brief 添加虚网格
 * @remark 本程序依赖ngg参数
 */
void Block::Add_ghostmesh(int32_t &ngg)
{
    int grmin[3] = {-ngg - 1, -ngg - 1, -ngg - 1};
    int grmax[3] = {mx + ngg + 1, my + ngg + 1, mz + ngg + 1}; // mx my mz 分别为ijk方向上从0开始的网格点（100个点则为99）
    int vrmin[3] = {0, 0, 0};
    int vrmax[3] = {mx, my, mz};

    int grange[3]; // 算上虚网格后，各个方向点数
    for (int i = 0; i < 3; i++)
        grange[i] = grmax[i] - grmin[i] + 1;

    // 加了虚网格之后的范围 从-ngg-1到mx+ngg+1，
    //[-ngg-1, -1]，[mx+1，mx+ngg+1]为虚网格,
    //[0, mx]为原网格

    // 外推虚网格之前先判断右手系z方向朝上还是朝下
    double xi = x(1, 0, 0) - x(0, 0, 0);
    double yi = y(1, 0, 0) - y(0, 0, 0);
    double xj = x(0, 1, 0) - x(0, 0, 0);
    double yj = y(0, 1, 0) - y(0, 0, 0);
    double dir = xi * yj - xj * yi;
    double ghostdz = 1.0;
    if (dir < 0)
        ghostdz = -1.0;

    // 开始推虚网格
    // 两部分：（1）非拐角点，（2）拐角点
    for (int gh = 1; gh <= ngg + 1; gh++)
    {
        // 每一步都只往外推1层虚网格

        int gimin = 0 - gh;
        int gimax = mx + gh;
        int gjmin = 0 - gh;
        int gjmax = my + gh;
        int gkmin = 0 - gh;
        int gkmax = mz + gh;

        // （1）非角点
        /*i方向*/
        for (int j = gjmin + 1; j < gjmax; j++)
        {
            for (int k = gkmin + 1; k < gkmax; k++)
            {
                x(gimin, j, k) = 2 * x(0, j, k) - x(-gimin, j, k);
                y(gimin, j, k) = 2 * y(0, j, k) - y(-gimin, j, k);
                z(gimin, j, k) = 2 * z(0, j, k) - z(-gimin, j, k);

                x(gimax, j, k) = 2 * x(mx, j, k) - x(2 * mx - gimax, j, k);
                y(gimax, j, k) = 2 * y(mx, j, k) - y(2 * mx - gimax, j, k);
                z(gimax, j, k) = 2 * z(mx, j, k) - z(2 * mx - gimax, j, k);
            }
        }
        /*j方向*/
        for (int i = gimin + 1; i < gimax; i++)
        {
            for (int k = gkmin + 1; k < gkmax; k++)
            {
                x(i, gjmin, k) = 2 * x(i, 0, k) - x(i, -gjmin, k);
                y(i, gjmin, k) = 2 * y(i, 0, k) - y(i, -gjmin, k);
                z(i, gjmin, k) = 2 * z(i, 0, k) - z(i, -gjmin, k);

                x(i, gjmax, k) = 2 * x(i, my, k) - x(i, 2 * my - gjmax, k);
                y(i, gjmax, k) = 2 * y(i, my, k) - y(i, 2 * my - gjmax, k);
                z(i, gjmax, k) = 2 * z(i, my, k) - z(i, 2 * my - gjmax, k);
            }
        }
        /*k方向*/
        if (dimension == 3)
        {
            for (int i = gimin + 1; i < gimax; i++)
            {
                for (int j = gjmin + 1; j < gjmax; j++)
                {
                    x(i, j, gkmin) = 2 * x(i, j, 0) - x(i, j, -gkmin);
                    y(i, j, gkmin) = 2 * y(i, j, 0) - y(i, j, -gkmin);
                    z(i, j, gkmin) = 2 * z(i, j, 0) - z(i, j, -gkmin);

                    x(i, j, gkmax) = 2 * x(i, j, mz) - x(i, j, 2 * mz - gkmax);
                    y(i, j, gkmax) = 2 * y(i, j, mz) - y(i, j, 2 * mz - gkmax);
                    z(i, j, gkmax) = 2 * z(i, j, mz) - z(i, j, 2 * mz - gkmax);
                }
            }
        }
        else
        {
            for (int i = gimin + 1; i < gimax; i++)
            {
                for (int j = gjmin + 1; j < gjmax; j++)
                {
                    x(i, j, gkmin) = x(i, j, gkmin + 1);
                    y(i, j, gkmin) = y(i, j, gkmin + 1);
                    z(i, j, gkmin) = z(i, j, gkmin + 1) - ghostdz;

                    x(i, j, gkmax) = x(i, j, gkmax - 1);
                    y(i, j, gkmax) = y(i, j, gkmax - 1);
                    z(i, j, gkmax) = z(i, j, gkmax - 1) + ghostdz;
                }
            }
        }

        // （2）边区和角点：二维4个角点再在k方向两侧各推ngg层，三维8个角点8条边区
        /*k方向，4个点/4条边区*/
        for (int k = gkmin + 1; k < gkmax; k++)
        {
            /*imin，jmin，k*/
            x(gimin, gjmin, k) = x(gimin + 1, gjmin, k) + x(gimin, gjmin + 1, k) - x(gimin + 1, gjmin + 1, k);
            y(gimin, gjmin, k) = y(gimin + 1, gjmin, k) + y(gimin, gjmin + 1, k) - y(gimin + 1, gjmin + 1, k);
            z(gimin, gjmin, k) = z(gimin + 1, gjmin, k) + z(gimin, gjmin + 1, k) - z(gimin + 1, gjmin + 1, k);

            /*imax，jmin，k*/
            x(gimax, gjmin, k) = x(gimax - 1, gjmin, k) + x(gimax, gjmin + 1, k) - x(gimax - 1, gjmin + 1, k);
            y(gimax, gjmin, k) = y(gimax - 1, gjmin, k) + y(gimax, gjmin + 1, k) - y(gimax - 1, gjmin + 1, k);
            z(gimax, gjmin, k) = z(gimax - 1, gjmin, k) + z(gimax, gjmin + 1, k) - z(gimax - 1, gjmin + 1, k);

            /*imin，jmax，k*/
            x(gimin, gjmax, k) = x(gimin + 1, gjmax, k) + x(gimin, gjmax - 1, k) - x(gimin + 1, gjmax - 1, k);
            y(gimin, gjmax, k) = y(gimin + 1, gjmax, k) + y(gimin, gjmax - 1, k) - y(gimin + 1, gjmax - 1, k);
            z(gimin, gjmax, k) = z(gimin + 1, gjmax, k) + z(gimin, gjmax - 1, k) - z(gimin + 1, gjmax - 1, k);

            /*imax，jmax，k*/
            x(gimax, gjmax, k) = x(gimax - 1, gjmax, k) + x(gimax, gjmax - 1, k) - x(gimax - 1, gjmax - 1, k);
            y(gimax, gjmax, k) = y(gimax - 1, gjmax, k) + y(gimax, gjmax - 1, k) - y(gimax - 1, gjmax - 1, k);
            z(gimax, gjmax, k) = z(gimax - 1, gjmax, k) + z(gimax, gjmax - 1, k) - z(gimax - 1, gjmax - 1, k);
        }
        if (dimension == 2)
        {
            /*j方向，4条边区*/
            for (int j = gjmin + 1; j < gjmax; j++)
            {
                /*imin，j，kmin*/
                x(gimin, j, gkmin) = x(gimin, j, gkmin + 1);
                y(gimin, j, gkmin) = y(gimin, j, gkmin + 1);
                z(gimin, j, gkmin) = z(gimin, j, gkmin + 1) - ghostdz;

                /*imax，j，kmin*/
                x(gimax, j, gkmin) = x(gimax, j, gkmin + 1);
                y(gimax, j, gkmin) = y(gimax, j, gkmin + 1);
                z(gimax, j, gkmin) = z(gimax, j, gkmin + 1) - ghostdz;

                /*imin，j，kmax*/
                x(gimin, j, gkmax) = x(gimin, j, gkmax - 1);
                y(gimin, j, gkmax) = y(gimin, j, gkmax - 1);
                z(gimin, j, gkmax) = z(gimin, j, gkmax - 1) + ghostdz;

                /*imax，j，kmax*/
                x(gimax, j, gkmax) = x(gimax, j, gkmax - 1);
                y(gimax, j, gkmax) = y(gimax, j, gkmax - 1);
                z(gimax, j, gkmax) = z(gimax, j, gkmax - 1) + ghostdz;
            }
            /*i方向，4条边区+8个角点*/
            for (int i = gimin; i < gimax + 1; i++)
            {
                /*i,gjmin,gkmin*/
                x(i, gjmin, gkmin) = x(i, gjmin, gkmin + 1);
                y(i, gjmin, gkmin) = y(i, gjmin, gkmin + 1);
                z(i, gjmin, gkmin) = z(i, gjmin, gkmin + 1) - ghostdz;

                /*i,gjmax,gkmin*/
                x(i, gjmax, gkmin) = x(i, gjmax, gkmin + 1);
                y(i, gjmax, gkmin) = y(i, gjmax, gkmin + 1);
                z(i, gjmax, gkmin) = z(i, gjmax, gkmin + 1) - ghostdz;

                /*i,gjmin,gkmax*/
                x(i, gjmin, gkmax) = x(i, gjmin, gkmax - 1);
                y(i, gjmin, gkmax) = y(i, gjmin, gkmax - 1);
                z(i, gjmin, gkmax) = z(i, gjmin, gkmax - 1) + ghostdz;

                /*i,gjmax.gkmax*/
                x(i, gjmax, gkmax) = x(i, gjmax, gkmax - 1);
                y(i, gjmax, gkmax) = y(i, gjmax, gkmax - 1);
                z(i, gjmax, gkmax) = z(i, gjmax, gkmax - 1) + ghostdz;
            }
        }
        else if (dimension == 3)
        {
            /*j方向，4条边区*/
            for (int j = gjmin + 1; j < gjmax; j++)
            {
                /*imin，j，kmin*/
                x(gimin, j, gkmin) = x(gimin + 1, j, gkmin) + x(gimin, j, gkmin + 1) - x(gimin + 1, j, gkmin + 1);
                y(gimin, j, gkmin) = y(gimin + 1, j, gkmin) + y(gimin, j, gkmin + 1) - y(gimin + 1, j, gkmin + 1);
                z(gimin, j, gkmin) = z(gimin + 1, j, gkmin) + z(gimin, j, gkmin + 1) - z(gimin + 1, j, gkmin + 1);

                /*imax，j，kmin*/
                x(gimax, j, gkmin) = x(gimax - 1, j, gkmin) + x(gimax, j, gkmin + 1) - x(gimax - 1, j, gkmin + 1);
                y(gimax, j, gkmin) = y(gimax - 1, j, gkmin) + y(gimax, j, gkmin + 1) - y(gimax - 1, j, gkmin + 1);
                z(gimax, j, gkmin) = z(gimax - 1, j, gkmin) + z(gimax, j, gkmin + 1) - z(gimax - 1, j, gkmin + 1);

                /*imin，j，kmax*/
                x(gimin, j, gkmax) = x(gimin + 1, j, gkmax) + x(gimin, j, gkmax - 1) - x(gimin + 1, j, gkmax - 1);
                y(gimin, j, gkmax) = y(gimin + 1, j, gkmax) + y(gimin, j, gkmax - 1) - y(gimin + 1, j, gkmax - 1);
                z(gimin, j, gkmax) = z(gimin + 1, j, gkmax) + z(gimin, j, gkmax - 1) - z(gimin + 1, j, gkmax - 1);

                /*imax，j，kmax*/
                x(gimax, j, gkmax) = x(gimax - 1, j, gkmax) + x(gimax, j, gkmax - 1) - x(gimax - 1, j, gkmax - 1);
                y(gimax, j, gkmax) = y(gimax - 1, j, gkmax) + y(gimax, j, gkmax - 1) - y(gimax - 1, j, gkmax - 1);
                z(gimax, j, gkmax) = z(gimax - 1, j, gkmax) + z(gimax, j, gkmax - 1) - z(gimax - 1, j, gkmax - 1);
            }
            /*i方向，4条边区+8个角点*/
            for (int i = gimin; i < gimax + 1; i++)
            {
                /*i,gjmin,gkmin*/
                x(i, gjmin, gkmin) = x(i, gjmin + 1, gkmin) + x(i, gjmin, gkmin + 1) - x(i, gjmin + 1, gkmin + 1);
                y(i, gjmin, gkmin) = y(i, gjmin + 1, gkmin) + y(i, gjmin, gkmin + 1) - y(i, gjmin + 1, gkmin + 1);
                z(i, gjmin, gkmin) = z(i, gjmin + 1, gkmin) + z(i, gjmin, gkmin + 1) - z(i, gjmin + 1, gkmin + 1);

                /*i,gjmax,gkmin*/
                x(i, gjmax, gkmin) = x(i, gjmax - 1, gkmin) + x(i, gjmax, gkmin + 1) - x(i, gjmax - 1, gkmin + 1);
                y(i, gjmax, gkmin) = y(i, gjmax - 1, gkmin) + y(i, gjmax, gkmin + 1) - y(i, gjmax - 1, gkmin + 1);
                z(i, gjmax, gkmin) = z(i, gjmax - 1, gkmin) + z(i, gjmax, gkmin + 1) - z(i, gjmax - 1, gkmin + 1);

                /*i,gjmin,gkmax*/
                x(i, gjmin, gkmax) = x(i, gjmin + 1, gkmax) + x(i, gjmin, gkmax - 1) - x(i, gjmin + 1, gkmax - 1);
                y(i, gjmin, gkmax) = y(i, gjmin + 1, gkmax) + y(i, gjmin, gkmax - 1) - y(i, gjmin + 1, gkmax - 1);
                z(i, gjmin, gkmax) = z(i, gjmin + 1, gkmax) + z(i, gjmin, gkmax - 1) - z(i, gjmin + 1, gkmax - 1);

                /*i,gjmax.gkmax*/
                x(i, gjmax, gkmax) = x(i, gjmax - 1, gkmax) + x(i, gjmax, gkmax - 1) - x(i, gjmax - 1, gkmax - 1);
                y(i, gjmax, gkmax) = y(i, gjmax - 1, gkmax) + y(i, gjmax, gkmax - 1) - y(i, gjmax - 1, gkmax - 1);
                z(i, gjmax, gkmax) = z(i, gjmax - 1, gkmax) + z(i, gjmax, gkmax - 1) - z(i, gjmax - 1, gkmax - 1);
            }
        }
    }

    /*输出查看*/
#ifdef test_ghost_mesh

    FILE *p;
    int32_t myid;
    DATATRANS::mpi_rank(&myid);
    std::string filename = "./ghostmesh" + std::to_string(myid) + ".dat";
    p = fopen(filename.c_str(), "w");
    fprintf(p, " variables=\"x\",\"y\",\"z\" \n");
    fprintf(p, "zone, i = %d, j = %d, k = %d, f = point \n", imax + 2 * (ngg + 1), jmax + 2 * (ngg + 1), kmax + 2 * (ngg + 1));
    for (int k = grmin[2]; k <= grmax[2]; k++)
    {
        for (int j = grmin[1]; j <= grmax[1]; j++)
        {
            for (int i = grmin[0]; i <= grmax[0]; i++)
            {
                fprintf(p, "%f	%f	%f \n", x(i, j, k), y(i, j, k), z(i, j, k));
            }
        }
    }
    fclose(p);
#endif // test_ghost_mesh
}

/**
 * @brief 计算Jacob行列式、度量系数
 * @remark 本程序依赖ngg参数
 */
void Block::calc_Jacobi_and_Metrics(int32_t &ngg)
{
    int rmin[3], rmax[3]; // 内部点坐标范围
    rmin[0] = 0 - ngg;
    rmin[1] = 0 - ngg;
    rmin[2] = 0 - ngg;
    rmax[0] = mx + ngg + 1;
    rmax[1] = my + ngg + 1;
    rmax[2] = mz + ngg + 1;
    // 有度量系数和Jacob的网格点从-ngg开始的到mx+nggj即imax+ngg,这里多加一个便于for循环中用<号

    jacobi.SetSize(1 + mx + 2 * ngg, 1 + my + 2 * ngg, 1 + mz + 2 * ngg, ngg);
    metric.SetSize(1 + mx + 2 * ngg, 1 + my + 2 * ngg, 1 + mz + 2 * ngg, ngg, 3, 3);
    inverse_metric.SetSize(1 + mx + 2 * ngg, 1 + my + 2 * ngg, 1 + mz + 2 * ngg, ngg, 3, 3);

    /*初始化*/
    jacobi = 1;
    metric = 0;
    inverse_metric = 0;

    // d(x,y,z)/d(1，2，3)
    int ghostdir[3] = {0};
    // bool isghostmesh = false;
    bool isghost_i, isghost_j, isghost_k;

    for (int ii = rmin[0]; ii < rmax[0]; ii++)
    {
        int i = ii;
        if (ii < 2)
            i = -ngg - ii + 1;
        // 将-ngg到0的顺序调整一下从而能够先算0再算-1、-2、...、-ngg，这样就能用已算的修正虚网格负Jacob的点

        if (i >= mx)
            ghostdir[0] = 1;
        // 大号面虚网格方向为正
        else if (i <= 0)
            ghostdir[0] = -1;
        // 小号面虚网格方向为负
        else
            ghostdir[0] = 0;

        isghost_i = bool(ghostdir[0]);

        for (int jj = rmin[1]; jj < rmax[1]; jj++)
        {
            int j = jj;
            if (jj < 2)
                j = -ngg - jj + 1;
            // 将-ngg到0的顺序调整一下从而能够先算0再算-1、-2、...、-ngg，这样就能用已算的修正虚网格负Jacob的点

            if (j >= my)
                ghostdir[1] = 1;
            // 大号面虚网格方向为正
            else if (j <= 0)
                ghostdir[1] = -1;
            // 小号面虚网格方向为负
            else
                ghostdir[1] = 0;

            isghost_j = bool(ghostdir[1]);

            for (int kk = rmin[2]; kk < rmax[2]; kk++)
            {
                int k = kk;
                if (kk < 2)
                    k = -ngg - kk + 1;
                // 将-ngg到0的顺序调整一下从而能够先算0再算-1、-2、...、-ngg，这样就能用已算的修正虚网格负Jacob的点

                if (k >= mz)
                    ghostdir[2] = 1;
                // 大号面虚网格方向为正
                else if (k <= 0)
                    ghostdir[2] = -1;
                // 小号面虚网格方向为负
                else
                    ghostdir[2] = 0;

                isghost_k = bool(ghostdir[2]);

                int iu1 = i + 1, ju1 = j + 1, ku1 = k + 1, il1 = i - 1, jl1 = j - 1, kl1 = k - 1; // u:upper;l:lower
                // double h[3] = { 2,2,2 };

                double dxd1 = (x(iu1, j, k) - x(il1, j, k)) / 2.0;
                double dxd2 = (x(i, ju1, k) - x(i, jl1, k)) / 2.0;
                double dxd3 = (x(i, j, ku1) - x(i, j, kl1)) / 2.0;
                double dyd1 = (y(iu1, j, k) - y(il1, j, k)) / 2.0;
                double dyd2 = (y(i, ju1, k) - y(i, jl1, k)) / 2.0;
                double dyd3 = (y(i, j, ku1) - y(i, j, kl1)) / 2.0;
                double dzd1 = (z(iu1, j, k) - z(il1, j, k)) / 2.0;
                double dzd2 = (z(i, ju1, k) - z(i, jl1, k)) / 2.0;
                double dzd3 = (z(i, j, ku1) - z(i, j, kl1)) / 2.0;

                inverse_metric(i, j, k, 0, 0) = dxd1;
                inverse_metric(i, j, k, 0, 1) = dxd2;
                inverse_metric(i, j, k, 0, 2) = dxd3;
                inverse_metric(i, j, k, 1, 0) = dyd1;
                inverse_metric(i, j, k, 1, 1) = dyd2;
                inverse_metric(i, j, k, 1, 2) = dyd3;
                inverse_metric(i, j, k, 2, 0) = dzd1;
                inverse_metric(i, j, k, 2, 1) = dzd2;
                inverse_metric(i, j, k, 2, 2) = dzd3;

                double j1 = dxd1 * dyd2 * dzd3;
                double j2 = dxd2 * dyd3 * dzd1;
                double j3 = dxd3 * dyd1 * dzd2;
                double j4 = dxd3 * dyd2 * dzd1;
                double j5 = dxd1 * dyd3 * dzd2;
                double j6 = dxd2 * dyd1 * dzd3;
                double ja = j1 + j2 + j3 - j4 - j5 - j6;

                double M1x = (dyd2 * dzd3 - dyd3 * dzd2) / ja;
                double M2x = -(dyd1 * dzd3 - dyd3 * dzd1) / ja;
                double M3x = (dyd1 * dzd2 - dyd2 * dzd1) / ja;
                double M1y = -(dxd2 * dzd3 - dxd3 * dzd2) / ja;
                double M2y = (dxd1 * dzd3 - dxd3 * dzd1) / ja;
                double M3y = -(dxd1 * dzd2 - dxd2 * dzd1) / ja;
                double M1z = (dxd2 * dyd3 - dxd3 * dyd2) / ja;
                double M2z = -(dxd1 * dyd3 - dxd3 * dyd1) / ja;
                double M3z = (dxd1 * dyd2 - dxd2 * dyd1) / ja;

                // dξ/dx dξ/dy dξ/dz
                metric(i, j, k, 0, 0) = M1x;
                metric(i, j, k, 0, 1) = M1y;
                metric(i, j, k, 0, 2) = M1z;

                // d_eta/dx d_eta/dy d_eta/dz
                metric(i, j, k, 1, 0) = M2x;
                metric(i, j, k, 1, 1) = M2y;
                metric(i, j, k, 1, 2) = M2z;

                // d_zeta/dx d_zeta/dy d_zeta/dz
                metric(i, j, k, 2, 0) = M3x;
                metric(i, j, k, 2, 1) = M3y;
                metric(i, j, k, 2, 2) = M3z;

                if (ja < 1e-20)
                {
                    if (isghost_i || isghost_j || isghost_k)
                    {
                        ja = jacobi(i - ghostdir[0], j - ghostdir[1], k - ghostdir[2]);
                        for (int32_t index_i = 0; index_i < 3; index_i++)
                            for (int32_t index_j = 0; index_j < 3; index_j++)
                                metric(i, j, k, index_i, index_j) = metric(i - ghostdir[0], j - ghostdir[1], k - ghostdir[2], index_i, index_j);
                    }
                    else
                    {
                        // std::cout << std::endl
                        //           << "*****************************************" << std::endl
                        //           << std::endl;
                        // std::cout << "The determination of Jacob is negative!Check The Grids!" << std::endl;
                        // std::cout << "size=\t" << mx << my << mz << std::endl;
                        // std::cout << i << "\t" << j << "\t" << k << "\t" << ja << std::endl;
                        // std::cout << x(i, j, k) << "\t" << y(i, j, k) << "\t" << x(-1, j, k) << "\t" << x(1, j, k) << y(i, -1, k) << "\t" << y(i, 1, k) << std::endl;
                        // std::cout << dzd3 << "\t" << dzd2 << "\t" << dzd1 << "\t" << dxd1 << "\t" << dyd2 << std::endl
                        //           << std::endl;
                        // std::cout << "*****************************************" << std::endl;
                        // exit(0);
                    }
                }

                jacobi(i, j, k) = ja;
            }
        }
    }
}

void Grid::MeshTrans_Parallel()
{
    GRID_TRANS::Para_data_transfer DATA;
    DATA.Initial_Allocate_Grid(this, 1);
    //===============================================================================
    int32_t index, num_parallel;
    //===============================================================================
    //----------------------------------------------------------------------
    // 收集所有数据
    index = 0;
    for (int i = 0; i < nblock; i++)
    {
        num_parallel = grids(i).parallel_bc.size();
        for (int j = 0; j < num_parallel; j++)
        {
            // 周期边界条件send recv flag均用奇数表示,故而只传偶数
            if (fmod(grids(i).parallel_bc[j].send_flag, 2) != 0)
                continue;
            // 修改2026/01/14 16：04 耦合面也需要传值
            // if (grids(i).parallel_bc[j].target_block_name != grids(i).parallel_bc[j].this_block_name)
            // {
            //     index++;
            //     continue;
            // }
            GRID_TRANS::Parallel_send_scalar(&grids(i).parallel_bc[j], ngg + 1, grids(i).x, DATA.buf_send[index]);
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
            // 周期边界条件send recv flag均用奇数表示,故而只传偶数
            if (fmod(grids(i).parallel_bc[j].send_flag, 2) != 0)
                continue;
            // 修改2026/01/14 16：04 耦合面也需要传值
            // if (grids(i).parallel_bc[j].target_block_name != grids(i).parallel_bc[j].this_block_name)
            // {
            //     DATA.request_s[index] = MPI_REQUEST_NULL;
            //     DATA.request_r[index] = MPI_REQUEST_NULL;
            //     index++;
            //     continue;
            // }
            PARALLEL::mpi_data_send(grids(i).parallel_bc[j].tar_myid, grids(i).parallel_bc[j].send_flag, DATA.buf_send[index], DATA.length[index], &DATA.request_s[index]);
            PARALLEL::mpi_data_recv(grids(i).parallel_bc[j].tar_myid, grids(i).parallel_bc[j].rece_flag, DATA.buf_recv[index], DATA.length[index], &DATA.request_r[index]);
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
            // 周期边界条件send recv flag均用奇数表示,故而只传偶数
            if (fmod(grids(i).parallel_bc[j].send_flag, 2) != 0)
                continue;
            // 修改2026/01/14 16：04 耦合面也需要传值
            // if (grids(i).parallel_bc[j].target_block_name != grids(i).parallel_bc[j].this_block_name)
            // {
            //     index++;
            //     continue;
            // }
            GRID_TRANS::Parallel_recv_scalar(&grids(i).parallel_bc[j], ngg + 1, grids(i).x, DATA.buf_recv[index]);
            index++;
        }
    }
    PARALLEL::mpi_barrier();
    //===============================================================================
    //===============================================================================
    //----------------------------------------------------------------------
    // 收集所有数据
    index = 0;
    for (int i = 0; i < nblock; i++)
    {
        num_parallel = grids(i).parallel_bc.size();
        for (int j = 0; j < num_parallel; j++)
        {
            // 周期边界条件send recv flag均用奇数表示,故而只传偶数
            if (fmod(grids(i).parallel_bc[j].send_flag, 2) != 0)
                continue;
            // 修改2026/01/14 16：04 耦合面也需要传值
            // if (grids(i).parallel_bc[j].target_block_name != grids(i).parallel_bc[j].this_block_name)
            // {
            //     index++;
            //     continue;
            // }
            GRID_TRANS::Parallel_send_scalar(&grids(i).parallel_bc[j], ngg + 1, grids(i).y, DATA.buf_send[index]);
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
            // 周期边界条件send recv flag均用奇数表示,故而只传偶数
            if (fmod(grids(i).parallel_bc[j].send_flag, 2) != 0)
                continue;
            // 修改2026/01/14 16：04 耦合面也需要传值
            // if (grids(i).parallel_bc[j].target_block_name != grids(i).parallel_bc[j].this_block_name)
            // {
            //     DATA.request_s[index] = MPI_REQUEST_NULL;
            //     DATA.request_r[index] = MPI_REQUEST_NULL;
            //     index++;
            //     continue;
            // }
            PARALLEL::mpi_data_send(grids(i).parallel_bc[j].tar_myid, grids(i).parallel_bc[j].send_flag, DATA.buf_send[index], DATA.length[index], &DATA.request_s[index]);
            PARALLEL::mpi_data_recv(grids(i).parallel_bc[j].tar_myid, grids(i).parallel_bc[j].rece_flag, DATA.buf_recv[index], DATA.length[index], &DATA.request_r[index]);
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
            // 周期边界条件send recv flag均用奇数表示,故而只传偶数
            if (fmod(grids(i).parallel_bc[j].send_flag, 2) != 0)
                continue;
            // 修改2026/01/14 16：04 耦合面也需要传值
            // if (grids(i).parallel_bc[j].target_block_name != grids(i).parallel_bc[j].this_block_name)
            // {
            //     index++;
            //     continue;
            // }
            GRID_TRANS::Parallel_recv_scalar(&grids(i).parallel_bc[j], ngg + 1, grids(i).y, DATA.buf_recv[index]);
            index++;
        }
    }
    PARALLEL::mpi_barrier();
    //===============================================================================
    //===============================================================================
    //----------------------------------------------------------------------
    // 收集所有数据
    index = 0;
    for (int i = 0; i < nblock; i++)
    {
        num_parallel = grids(i).parallel_bc.size();
        for (int j = 0; j < num_parallel; j++)
        {
            // 周期边界条件send recv flag均用奇数表示,故而只传偶数
            if (fmod(grids(i).parallel_bc[j].send_flag, 2) != 0)
                continue;
            // 修改2026/01/14 16：04 耦合面也需要传值
            // if (grids(i).parallel_bc[j].target_block_name != grids(i).parallel_bc[j].this_block_name)
            // {
            //     index++;
            //     continue;
            // }
            GRID_TRANS::Parallel_send_scalar(&grids(i).parallel_bc[j], ngg + 1, grids(i).z, DATA.buf_send[index]);
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
            // 周期边界条件send recv flag均用奇数表示,故而只传偶数
            if (fmod(grids(i).parallel_bc[j].send_flag, 2) != 0)
                continue;
            // 修改2026/01/14 16：04 耦合面也需要传值
            // if (grids(i).parallel_bc[j].target_block_name != grids(i).parallel_bc[j].this_block_name)
            // {
            //     DATA.request_s[index] = MPI_REQUEST_NULL;
            //     DATA.request_r[index] = MPI_REQUEST_NULL;
            //     index++;
            //     continue;
            // }
            PARALLEL::mpi_data_send(grids(i).parallel_bc[j].tar_myid, grids(i).parallel_bc[j].send_flag, DATA.buf_send[index], DATA.length[index], &DATA.request_s[index]);
            PARALLEL::mpi_data_recv(grids(i).parallel_bc[j].tar_myid, grids(i).parallel_bc[j].rece_flag, DATA.buf_recv[index], DATA.length[index], &DATA.request_r[index]);
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
            // 周期边界条件send recv flag均用奇数表示,故而只传偶数
            if (fmod(grids(i).parallel_bc[j].send_flag, 2) != 0)
                continue;
            // 修改2026/01/14 16：04 耦合面也需要传值
            // if (grids(i).parallel_bc[j].target_block_name != grids(i).parallel_bc[j].this_block_name)
            // {
            //     index++;
            //     continue;
            // }
            GRID_TRANS::Parallel_recv_scalar(&grids(i).parallel_bc[j], ngg + 1, grids(i).z, DATA.buf_recv[index]);
            index++;
        }
    }
    PARALLEL::mpi_barrier();
    //===============================================================================
    // 释放空间
    DATA.Deallocate();
    PARALLEL::mpi_barrier();
};

void Grid::MeshTrans_Inner()
{
    int32_t num_inner;
    int32_t tar_block_num;
    // 由于inner传值每次是处理一对，因此需要用inner_bc中的index标记是否已经处理，避免重复传值
    // 不仅仅是计算量的问题，重复传值会导致边界面的第0层或第max层数据错误
    // 将所有的index初始化为0
    for (int i = 0; i < nblock; i++)
    {
        num_inner = grids(i).inner_bc.size();
        for (int j = 0; j < num_inner; j++)
            grids(i).inner_bc[j].index = 0;
    }
    // 内部传值
    for (int i = 0; i < nblock; i++)
    {
        num_inner = grids(i).inner_bc.size();
        for (int j = 0; j < num_inner; j++)
        {
            // 块号小于零，表示周期边界条件，网格不传
            // if (grids(i).inner_bc[j].this_block_num < 0)
            // Modify(2024/04/26/17:13)该数为正，从0开始，是否为周期信息存在is_period
            if (grids(i).inner_bc[j].is_period)
                continue;
            // 若此inner面的index不为0，表明已经传值处理过了
            if (grids(i).inner_bc[j].index != 0)
                continue;
            // 耦合壁面不需要传 // 修改2026/01/14 16：04 耦合面也需要传值
            // if (grids(i).inner_bc[j].target_block_name != grids(i).inner_bc[j].this_block_name)
            //     continue;
            //---------------------------------------------------------------------------------
            // Modify(2024/04/26/17:13)该数为正，从0开始，是否为周期信息存在is_period
            tar_block_num = grids(i).inner_bc[j].tar_block_num;
            // tar_block_num = grids(i).inner_bc[j].tar_block_num - 1;

            GRID_TRANS::Inner_trans_scalar(&grids(i).inner_bc[j], ngg + 1, grids(i).x, grids(tar_block_num).x);
            GRID_TRANS::Inner_trans_scalar(&grids(i).inner_bc[j], ngg + 1, grids(i).y, grids(tar_block_num).y);
            GRID_TRANS::Inner_trans_scalar(&grids(i).inner_bc[j], ngg + 1, grids(i).z, grids(tar_block_num).z);
            // 将此inner面的index标记为1，表明已经传值处理过了
            grids(i).inner_bc[j].index = 1;
            // 将目标块的inner面的index标也记为1，表明已经传值处理过了
            grids(tar_block_num).inner_bc[grids(i).inner_bc[j].tar_block_index].index = 1;
        }
    }
}

/**
 * @brief 假虚网格（Faker）区域的几何信息部分处理，分别对Physical, Inner, Parallel的Faker区域开辟空间
 *        物理边界条件直接复制本物理快的虚网格信息；
 *        然后通过调用Faker_Process_MeshTrans_Inner Faker_Process_MeshTrans_Parallel对内部、并行的
 *        物理块耦合界面虚网格进行处理
 */
void Grid::Faker_Process()
{
    for (int32_t izone = 0; izone < nblock; ++izone)
    {
        // 对每一个块的边界条件中Faker区域网格部分开辟空间
        Block &iblock = grids(izone);

        //-------------------------------------------------------------------------------
        // Physical_Boundary
        int32_t bc_number = iblock.physical_bc.size();
        int32_t length[3], sub[3], sup[3];
        for (int32_t index = 0; index < bc_number; index++)
        {
            Physical_Boundary &phybnd = iblock.physical_bc[index];
            for (int i = 0; i < 3; i++)
            {
                length[i] = phybnd.sup[i] - phybnd.sub[i] + 1;
                sub[i] = phybnd.sub[i];
                sup[i] = phybnd.sup[i];
            }
            length[abs(phybnd.direction) - 1] += ngg;
            if (phybnd.direction < 0)
                sub[abs(phybnd.direction) - 1] -= ngg;
            else
                sup[abs(phybnd.direction) - 1] += ngg;

            phybnd.Faker.Fakerx.SetSize(length[0], length[1], length[2], 0);
            phybnd.Faker.Fakery.SetSize(length[0], length[1], length[2], 0);
            phybnd.Faker.Fakerz.SetSize(length[0], length[1], length[2], 0);
            phybnd.Faker.Fakerjac.SetSize(length[0], length[1], length[2], 0);
            phybnd.Faker.Fakermetric.SetSize(length[0], length[1], length[2], 0, 3, 3);
            phybnd.Faker.is_multi_phys = false;

            for (int32_t i = sub[0]; i <= sup[0]; i++)
                for (int32_t j = sub[1]; j <= sup[1]; j++)
                    for (int32_t k = sub[2]; k <= sup[2]; k++)
                    {
                        phybnd.Faker.Fakerx(i - sub[0], j - sub[1], k - sub[2]) = iblock.x(i, j, k);
                        phybnd.Faker.Fakery(i - sub[0], j - sub[1], k - sub[2]) = iblock.y(i, j, k);
                        phybnd.Faker.Fakerz(i - sub[0], j - sub[1], k - sub[2]) = iblock.z(i, j, k);
                        phybnd.Faker.Fakerjac(i - sub[0], j - sub[1], k - sub[2]) = iblock.jacobi(i, j, k);
                        for (int32_t ll = 0; ll < 3; ll++)
                            for (int32_t kk = 0; kk < 3; kk++)
                                phybnd.Faker.Fakermetric(i - sub[0], j - sub[1], k - sub[2], ll, kk) = iblock.metric(i, j, k, ll, kk);
                    }
        }
        //-------------------------------------------------------------------------------

        //-------------------------------------------------------------------------------
        // Inner_Boundary
        bc_number = iblock.inner_bc.size();
        for (int32_t index = 0; index < bc_number; index++)
        {
            Inner_Boundary &Innerbnd = iblock.inner_bc[index];
            if (Innerbnd.this_block_name != Innerbnd.target_block_name)
                Innerbnd.Faker.is_multi_phys = true;
            else
            {
                Innerbnd.Faker.is_multi_phys = false;
                continue;
            }

            for (int i = 0; i < 3; i++)
            {
                sub[i] = fmin(abs(Innerbnd.sub[i]), abs(Innerbnd.sup[i]));
                sup[i] = fmax(abs(Innerbnd.sub[i]), abs(Innerbnd.sup[i]));
                length[i] = sup[i] - sub[i] + 1;
            }
            length[abs(Innerbnd.direction) - 1] += ngg;
            if (Innerbnd.direction < 0)
                sub[abs(Innerbnd.direction) - 1] -= ngg;
            else
                sup[abs(Innerbnd.direction) - 1] += ngg;

            Innerbnd.Faker.Fakerx.SetSize(length[0], length[1], length[2], 0);
            Innerbnd.Faker.Fakery.SetSize(length[0], length[1], length[2], 0);
            Innerbnd.Faker.Fakerz.SetSize(length[0], length[1], length[2], 0);
            Innerbnd.Faker.Fakerjac.SetSize(length[0], length[1], length[2], 0);
            Innerbnd.Faker.Fakermetric.SetSize(length[0], length[1], length[2], 0, 3, 3);
        }
        //-------------------------------------------------------------------------------

        //-------------------------------------------------------------------------------
        // Para_Boundary
        bc_number = iblock.parallel_bc.size();
        for (int32_t index = 0; index < bc_number; index++)
        {
            Parallel_Boundary &Parabnd = iblock.parallel_bc[index];

            if (Parabnd.this_block_name != Parabnd.target_block_name)
                Parabnd.Faker.is_multi_phys = true;
            else
            {
                Parabnd.Faker.is_multi_phys = false;
                continue;
            }

            for (int i = 0; i < 3; i++)
            {
                sub[i] = fmin(abs(Parabnd.sub[i]), abs(Parabnd.sup[i]));
                sup[i] = fmax(abs(Parabnd.sub[i]), abs(Parabnd.sup[i]));
                length[i] = sup[i] - sub[i] + 1;
            }
            length[abs(Parabnd.direction) - 1] += ngg;
            if (Parabnd.direction < 0)
                sub[abs(Parabnd.direction) - 1] -= ngg;
            else
                sup[abs(Parabnd.direction) - 1] += ngg;

            Parabnd.Faker.Fakerx.SetSize(length[0], length[1], length[2], 0);
            Parabnd.Faker.Fakery.SetSize(length[0], length[1], length[2], 0);
            Parabnd.Faker.Fakerz.SetSize(length[0], length[1], length[2], 0);
            Parabnd.Faker.Fakerjac.SetSize(length[0], length[1], length[2], 0);
            Parabnd.Faker.Fakermetric.SetSize(length[0], length[1], length[2], 0, 3, 3);
        }
        //-------------------------------------------------------------------------------
    }

    // 对Faker区域的几何参数进行处理
    Faker_Process_MeshTrans_Inner();

    Faker_Process_MeshTrans_Parallel();
}

/**
 * @brief 假虚网格区域的内部网格传值处理，直接使用DATATRANS::Inner_flush_scalar _tensor传递几何信息
 */
void Grid::Faker_Process_MeshTrans_Inner()
{
    int32_t bc_number;
    for (int32_t izone = 0; izone < nblock; ++izone)
    {
        // 对每一个块的边界条件中Faker区域处理
        Block &iblock = grids(izone);
        //-------------------------------------------------------------------------------
        // Inner_Boundary
        bc_number = iblock.inner_bc.size();
        for (int32_t index = 0; index < bc_number; index++)
        {
            Inner_Boundary &Innerbnd = iblock.inner_bc[index];
            // Block &tar_block = grids(Innerbnd.tar_block_num - 1);
            Block &tar_block = grids(Innerbnd.tar_block_num);
            if (!Innerbnd.Faker.is_multi_phys)
            {
                continue;
            }

            GRID_TRANS::Inner_flush_scalar(&Innerbnd, Innerbnd.Faker.Fakerx, ngg, tar_block.x);
            GRID_TRANS::Inner_flush_scalar(&Innerbnd, Innerbnd.Faker.Fakery, ngg, tar_block.y);
            GRID_TRANS::Inner_flush_scalar(&Innerbnd, Innerbnd.Faker.Fakerz, ngg, tar_block.z);
            GRID_TRANS::Inner_flush_scalar(&Innerbnd, Innerbnd.Faker.Fakerjac, ngg, tar_block.jacobi);

            GRID_TRANS::Inner_flush_tensor(&Innerbnd, Innerbnd.Faker.Fakermetric, ngg, tar_block.metric);
        }
        //-------------------------------------------------------------------------------
    }
}

/**
 * @brief 假虚网格区域的并行网格传值处理，利用Initial_Allocate_Couple_Grid开辟并行缓冲区，变量数设置为9
 *        分别传递虚网格的几何信息，存入Faker区域。这里9是为了能够一次性传递度量系数
 */
void Grid::Faker_Process_MeshTrans_Parallel()
{
    GRID_TRANS::Para_data_transfer DATA;
    DATA.Initial_Allocate_Couple_Grid(this, 9);
    //===============================================================================
    int32_t index, num_parallel;
    //===============================================================================
    //----------------------------------------------------------------------
    // 收集所有数据
    index = 0;
    for (int i = 0; i < nblock; i++)
    {
        num_parallel = grids(i).parallel_bc.size();
        for (int j = 0; j < num_parallel; j++)
        {
            // 这里只考虑耦合边界条件
            if (!grids(i).parallel_bc[j].Faker.is_multi_phys)
                continue;
            GRID_TRANS::Parallel_send_scalar(&grids(i).parallel_bc[j], ngg, grids(i).x, DATA.buf_send[index]);
            index++;
        }
    }
    //----------------------------------------------------------------------
    // 发送接受所有数据
    index = 0;
    for (int i = 0; i < nblock; i++)
    {
        num_parallel = grids(i).parallel_bc.size();
        for (int j = 0; j < num_parallel; j++)
        {
            // 这里只考虑耦合边界条件
            if (!grids(i).parallel_bc[j].Faker.is_multi_phys)
                continue;
            PARALLEL::mpi_data_send(grids(i).parallel_bc[j].tar_myid, grids(i).parallel_bc[j].send_flag, DATA.buf_send[index], DATA.length[index], &DATA.request_s[index]);
            PARALLEL::mpi_data_recv(grids(i).parallel_bc[j].tar_myid, grids(i).parallel_bc[j].rece_flag, DATA.buf_recv[index], DATA.length[index], &DATA.request_r[index]);
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
            // 这里只考虑耦合边界条件
            if (!grids(i).parallel_bc[j].Faker.is_multi_phys)
                continue;
            GRID_TRANS::Parallel_flush_recv_scalar(&grids(i).parallel_bc[j], grids(i).parallel_bc[j].Faker.Fakerx, ngg, DATA.buf_recv[index]);
            index++;
        }
    }
    //===============================================================================

    //===============================================================================
    //----------------------------------------------------------------------
    // 收集所有数据
    index = 0;
    for (int i = 0; i < nblock; i++)
    {
        num_parallel = grids(i).parallel_bc.size();
        for (int j = 0; j < num_parallel; j++)
        {
            // 这里只考虑耦合边界条件
            if (!grids(i).parallel_bc[j].Faker.is_multi_phys)
                continue;
            GRID_TRANS::Parallel_send_scalar(&grids(i).parallel_bc[j], ngg, grids(i).y, DATA.buf_send[index]);
            index++;
        }
    }
    //----------------------------------------------------------------------
    // 发送接受所有数据
    index = 0;
    for (int i = 0; i < nblock; i++)
    {
        num_parallel = grids(i).parallel_bc.size();
        for (int j = 0; j < num_parallel; j++)
        {
            // 这里只考虑耦合边界条件
            if (!grids(i).parallel_bc[j].Faker.is_multi_phys)
                continue;
            PARALLEL::mpi_data_send(grids(i).parallel_bc[j].tar_myid, grids(i).parallel_bc[j].send_flag, DATA.buf_send[index], DATA.length[index], &DATA.request_s[index]);
            PARALLEL::mpi_data_recv(grids(i).parallel_bc[j].tar_myid, grids(i).parallel_bc[j].rece_flag, DATA.buf_recv[index], DATA.length[index], &DATA.request_r[index]);
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
            // 这里只考虑耦合边界条件
            if (!grids(i).parallel_bc[j].Faker.is_multi_phys)
                continue;
            GRID_TRANS::Parallel_flush_recv_scalar(&grids(i).parallel_bc[j], grids(i).parallel_bc[j].Faker.Fakery, ngg, DATA.buf_recv[index]);
            index++;
        }
    }
    //===============================================================================

    //===============================================================================
    //----------------------------------------------------------------------
    // 收集所有数据
    index = 0;
    for (int i = 0; i < nblock; i++)
    {
        num_parallel = grids(i).parallel_bc.size();
        for (int j = 0; j < num_parallel; j++)
        {
            // 这里只考虑耦合边界条件
            if (!grids(i).parallel_bc[j].Faker.is_multi_phys)
                continue;
            GRID_TRANS::Parallel_send_scalar(&grids(i).parallel_bc[j], ngg, grids(i).z, DATA.buf_send[index]);
            index++;
        }
    }
    //----------------------------------------------------------------------
    // 发送接受所有数据
    index = 0;
    for (int i = 0; i < nblock; i++)
    {
        num_parallel = grids(i).parallel_bc.size();
        for (int j = 0; j < num_parallel; j++)
        {
            // 这里只考虑耦合边界条件
            if (!grids(i).parallel_bc[j].Faker.is_multi_phys)
                continue;
            PARALLEL::mpi_data_send(grids(i).parallel_bc[j].tar_myid, grids(i).parallel_bc[j].send_flag, DATA.buf_send[index], DATA.length[index], &DATA.request_s[index]);
            PARALLEL::mpi_data_recv(grids(i).parallel_bc[j].tar_myid, grids(i).parallel_bc[j].rece_flag, DATA.buf_recv[index], DATA.length[index], &DATA.request_r[index]);
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
            // 这里只考虑耦合边界条件
            if (!grids(i).parallel_bc[j].Faker.is_multi_phys)
                continue;
            GRID_TRANS::Parallel_flush_recv_scalar(&grids(i).parallel_bc[j], grids(i).parallel_bc[j].Faker.Fakerz, ngg, DATA.buf_recv[index]);
            index++;
        }
    }
    //===============================================================================

    //===============================================================================
    //----------------------------------------------------------------------
    // 收集所有数据
    index = 0;
    for (int i = 0; i < nblock; i++)
    {
        num_parallel = grids(i).parallel_bc.size();
        for (int j = 0; j < num_parallel; j++)
        {
            // 这里只考虑耦合边界条件
            if (!grids(i).parallel_bc[j].Faker.is_multi_phys)
                continue;
            GRID_TRANS::Parallel_send_scalar(&grids(i).parallel_bc[j], ngg, grids(i).jacobi, DATA.buf_send[index]);
            index++;
        }
    }
    //----------------------------------------------------------------------
    // 发送接受所有数据
    index = 0;
    for (int i = 0; i < nblock; i++)
    {
        num_parallel = grids(i).parallel_bc.size();
        for (int j = 0; j < num_parallel; j++)
        {
            // 这里只考虑耦合边界条件
            if (!grids(i).parallel_bc[j].Faker.is_multi_phys)
                continue;
            PARALLEL::mpi_data_send(grids(i).parallel_bc[j].tar_myid, grids(i).parallel_bc[j].send_flag, DATA.buf_send[index], DATA.length[index], &DATA.request_s[index]);
            PARALLEL::mpi_data_recv(grids(i).parallel_bc[j].tar_myid, grids(i).parallel_bc[j].rece_flag, DATA.buf_recv[index], DATA.length[index], &DATA.request_r[index]);
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
            // 这里只考虑耦合边界条件
            if (!grids(i).parallel_bc[j].Faker.is_multi_phys)
                continue;
            GRID_TRANS::Parallel_flush_recv_scalar(&grids(i).parallel_bc[j], grids(i).parallel_bc[j].Faker.Fakerjac, ngg, DATA.buf_recv[index]);
            index++;
        }
    }
    //===============================================================================

    //===============================================================================
    //----------------------------------------------------------------------
    // 收集所有数据
    index = 0;
    for (int i = 0; i < nblock; i++)
    {
        num_parallel = grids(i).parallel_bc.size();
        for (int j = 0; j < num_parallel; j++)
        {
            // 这里只考虑耦合边界条件
            if (!grids(i).parallel_bc[j].Faker.is_multi_phys)
                continue;
            GRID_TRANS::Parallel_send_tensor(&grids(i).parallel_bc[j], ngg, grids(i).metric, DATA.buf_send[index]);
            index++;
        }
    }
    //----------------------------------------------------------------------
    // 发送接受所有数据
    index = 0;
    for (int i = 0; i < nblock; i++)
    {
        num_parallel = grids(i).parallel_bc.size();
        for (int j = 0; j < num_parallel; j++)
        {
            // 这里只考虑耦合边界条件
            if (!grids(i).parallel_bc[j].Faker.is_multi_phys)
                continue;
            PARALLEL::mpi_data_send(grids(i).parallel_bc[j].tar_myid, grids(i).parallel_bc[j].send_flag, DATA.buf_send[index], DATA.length[index], &DATA.request_s[index]);
            PARALLEL::mpi_data_recv(grids(i).parallel_bc[j].tar_myid, grids(i).parallel_bc[j].rece_flag, DATA.buf_recv[index], DATA.length[index], &DATA.request_r[index]);
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
            // 这里只考虑耦合边界条件
            if (!grids(i).parallel_bc[j].Faker.is_multi_phys)
                continue;
            GRID_TRANS::Parallel_flush_recv_tensor(&grids(i).parallel_bc[j], grids(i).parallel_bc[j].Faker.Fakermetric, ngg, DATA.buf_recv[index]);
            index++;
        }
    }
    //===============================================================================

    //===============================================================================
    // 释放空间
    DATA.Deallocate();
    PARALLEL::mpi_barrier();
}
