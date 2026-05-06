/*****************************************************************************
@Copyright: NLCFD
@File name: 1_MPCNS_Grid.h
@Author: Descartes
@Version: 1.0
@Date: 2022年09月11日
@Description:	利用自定义的可变数组ARRAY_DEFIEN形成的用于计算的结构网格数据结构，所有的
                几何信息全部存储在Grid这一类中，其中包含了多个Block(即网格块)
@History:		（修改历史记录列表， 每条修改记录应包括修改日期、修改者及修改内容简述。）
                1、修改了文件名称，基本的功能头文件名用数字1+全部首字母大写的方式
                2、2022/10/04 16：59修改Inner Parallel 网格传值，inner id<0不传
                parallel的flag<0也不传
*****************************************************************************/
#pragma once
#include "0_basic/ARRAY_DEFINE.h"
#include "0_basic/1_MPCNS_Parameter.h"
#include "1_grid/Grid_Boundary.h"
#include <math.h>

/**
 * @brief			        Block类存储着所有的与网格相关的信息
 * @param dimension	        网格维度
 * @param imax jmax kmax	三个方向网格点数（未加虚网格）1~imax
 * @param mx   my   mz      三个方向的网格范围0~mx，故mx=imax-1
 * @param block_name        网格块名称，用于判断使用什么物理长求解器
 * @param x y z		        三个三维数组，分别存储网格点的x y z坐标
 * @param jacobi	        网格的jacobi行列式值
 * @param metric	        度量系数，1代表ξ，2代表η，3代表ζ，分别有dkdx dkdy dkdz三个分量
 * @param inverse_metric	逆度量系数，1代表ξ，2代表η，3代表ζ，分别有dxdk dydk dzdk三个分量
 * @param parallel_bc		并行边界信息
 * @param inner_bc	        内边界信息
 * @param physical_bc	    物理边界信息
 */
class Block
{
public:
    int32_t dimension;
    int32_t imax, jmax, kmax;
    int32_t mx, my, mz;
    std::string block_name;
    double3D x, y, z, jacobi;
    Phy_Tensor metric, inverse_metric;
    std::vector<Parallel_Boundary> parallel_bc;
    std::vector<Inner_Boundary> inner_bc;
    std::vector<Physical_Boundary> physical_bc;
    std::vector<Edge> edge;
    std::vector<std::vector<Edge>> corner, corner_send; // 在DATATrans Initial中处理
    Param *par;
    Phy_Vector GCL_metric_xi, GCL_metric_eta, GCL_metric_zeta; // 半点界面的Geometric Conservation Law的面度量系数

    Phy_Vector GCL_covar_xi_xi;   // (i+1/2,j,k)，ξ面的ξ方向切向协变向量
    Phy_Vector GCL_covar_xi_eta;  // (i+1/2,j,k)，ξ面的η方向切向协变向量
    Phy_Vector GCL_covar_xi_zeta; // (i+1/2,j,k)，ξ面的ζ方向切向协变向量
    Phy_Tensor GCL_covar_xi;      // (i+1/2,j,k)，ξ面的协变度规

    Phy_Vector GCL_covar_eta_xi;   // (i,j+1/2,k)，η面的ξ方向切向协变向量
    Phy_Vector GCL_covar_eta_eta;  // (i,j+1/2,k)，η面的η方向切向协变向量
    Phy_Vector GCL_covar_eta_zeta; // (i,j+1/2,k)，η面的ζ方向切向协变向量
    Phy_Tensor GCL_covar_eta;      // (i,j+1/2,k)，η面的协变度规

    Phy_Vector GCL_covar_zeta_xi;   // (i,j,k+1/2)，ζ面的ξ方向切向协变向量
    Phy_Vector GCL_covar_zeta_eta;  // (i,j,k+1/2)，ζ面的η方向切向协变向量
    Phy_Vector GCL_covar_zeta_zeta; // (i,j,k+1/2)，ζ面的ζ方向切向协变向量
    Phy_Tensor GCL_covar_zeta;      // (i,j,k+1/2)，ζ面的协变度规

    double3D dual_x, dual_y, dual_z; // dual grid的i，j，k都向正向偏移半个单位即为正常网格坐标

public:
    Block() {};
    ~Block() {};

    /**
     * @brief   为该Block添加虚网格
     */
    void Add_ghostmesh(int32_t &ngg);

    /**
     * @brief   为该Block计算Jacobi和度量系数
     *          计算出雅可比行列式为负时需要进行判断，若是虚网格雅可比为负，则该点雅可比等于内层点的雅可比。
     *	        若为网格点雅可比为负，说明网格方向有误，需要修改网格ξ η ζ的方向，满足右手系。
     */
    void calc_Jacobi_and_Metrics(int32_t &ngg);

    /**
     * @brief   为该Block计算GCL的半点度量系数
     */
    void calc_Metrics_GCL(int32_t &ngg);

    /**
     * @brief   为该Block计算GCL的Jacobi
     */
    void calc_modify_Jacobi(int32_t &ngg);

    /**
     * @brief   为该Block计算GCL的半点切向度量系数
     */
    void calc_Face_Tangent_Vectors(int32_t &ngg);

    /**
     * @brief   为该Block计算对偶网格dual grids
     */
    void calc_Dual_Grids();

    /**
     * @brief 针对该Block的角区进行预处理，记录哪些边界涉及角区处理
     */
    void calc_Corner_Preprocess();

    void create_edges_from_boundary(const int32_t sub[3], const int32_t sup[3], int32_t direction,
                                    Parallel_Boundary *para_bound, Inner_Boundary *inner_bound, Physical_Boundary *phy_bound);
    void create_edges_from_boundary(const int32_t sub[3], const int32_t sup[3], int32_t direction,
                                    Parallel_Boundary *para_bound, Inner_Boundary *inner_bound, Physical_Boundary *phy_bound, int _index);

    // 用于比较两条边是否重叠或包含（方向无关）
    bool is_overlapping(const Edge &e1, const Edge &e2, int &dir);
    bool is_same(const Edge &e1, const Edge &e2);

    void add_edge(Edge &e1);
    void add_edge(Edge &e1, int _index);

    void integration_equal_edge(Edge &e2, Edge &e1);

    void judge_edges();
    void judge_edges(int _index);
};

/**
 * @brief			        Grid类存储着所有的与网格相关的信息
 * @param dimension	        网格维度
 * @param nblock	        本进程中网格块的个数
 * @param ngg	            虚网格个数
 * @param grids	            Block组成的一维数组，存储所有数据
 * @param boundary_name_index 根据名称（string）找序号
 * @param boundary_name     根据序号找名称（string）
 */
class Grid
{
public:
    Param *par;
    int32_t dimension;
    int32_t nblock; // 线程内block数
    int32_t ngg;    // 虚网格个数
    Array1D<Block> grids;
    std::map<std::string, int32_t> boundary_name_index;
    std::vector<std::string> boundary_name;

public:
    Grid() {};
    ~Grid() {};

public:
    void Grid_Preprocess(Param *_par)
    {
        par = _par;
        int32_t _myid = par->GetInt("myid");
        int32_t _ngg = par->GetInt("ngg");
        double _scale = par->GetDou("scale");
        int32_t Read_Method = par->GetInt("Read_Method");

        //-----------------------------------------------------------------------------------------
        if (_myid == 0)
            std::cout << "---->Reading the Grids...\n";
        ReadGeometry(_myid, _ngg, Read_Method);

        if (_myid == 0)
            std::cout << "*****************Finish the Grid reading Process! !*******************\n\n";
        //-----------------------------------------------------------------------------------------
        if (_myid == 0)
            std::cout << "---->Processing the Grids...\n";
        if (_myid == 0)
            std::cout << "\t-->scaling...\n";
        scale_mesh(_scale);

        if (_myid == 0)
            std::cout << "\t-->Adding ghost mesh...\n";
        Get_ghostmesh();

        if (_myid == 0)
            std::cout << "\t-->Synchronizing Inner interface mesh...\n";
        MeshTrans_Inner();

        if (_myid == 0)
            std::cout << "\t-->Synchronizing Paral interface mesh...\n";
        MeshTrans_Parallel();

        if (_myid == 0)
            std::cout << "\t-->Adding Pole ghost mesh...\n";
        MeshTrans_Inner_Pole();

        if (par->GetInt("corner_treat") == 2)
        {
            if (_myid == 0)
                std::cout << "---->Processing Corner...\n";
            for (int i = 0; i < nblock; i++)
                grids(i).par = par;
            Corner_Preprocess();

            if (_myid == 0)
                std::cout << "\t-->Synchronizing Inner Corner mesh...\n";
            MeshTrans_Corner_Inner();

            if (_myid == 0)
                std::cout << "\t-->Synchronizing Paral Corner mesh...\n";
            MeshTrans_Corner_Parallel();

            if (_myid == 0)
                std::cout << "\t-->Synchronizing Inner Corner3D mesh...\n";
            MeshTrans_Corner3D_Inner();

            if (_myid == 0)
                std::cout << "\t-->Synchronizing Paral Corner3D mesh...\n";
            MeshTrans_Corner3D_Parallel();
        }

        if (_myid == 0)
            std::cout << "\t-->Calculating Jacobi and Metrics...\n";
        Calc_Jacobi_and_Metrics();

        if (_myid == 0)
            std::cout << "\t-->Processing the Faker_Boundary Part...\n";
        Faker_Process();

        if (_myid == 0)
            std::cout << "***************Finish the Grid Manipulating Process! !****************\n\n";
        //-----------------------------------------------------------------------------------------
    }

private:
    //=============================================================================================
    /**
     * @brief 读取plot3D格式文件，grid n.grd、boundary n.txt、inner n.txt、parall n.txt、boundary.fvbnd
     *        ------->针对grids
     *        使用了##ngg##参数       设置网格数组尺度;
     *        设置Grid的dimesion nblock grids的长度;
     *        读入fvbnd文件，获得数组boundary_name_index boundary_name;
     *        ------->针对每一个网格块block
     *        设置网格名称block_name;
     *        设置网格的尺寸i，j，kmax mx my mz;
     *        设置网格的维度dimension;
     *        设置网格坐标xyz尺寸 jacobi尺寸 metric inver_metric尺寸;
     *        读入坐标xyz;
     *       【虚网格xyz jacobi 度规 在这计算】;
     *        -------->针对每一个网格块block中的physical_bc
     *        利用fvbnd获取的boundary_name_index boundary_name;
     *        读入parallel_bc中的起始 终止坐标sub sup;
     *        获取direction boundary_num this_block_num boundary_name;
     *        -------->针对每一个网格块block中的inner_bc
     *        每读入一个内部边界条件就push_back进inner_bc中
     *        其中包含了本块的边界起始位置，目标块的起始位置（利用Pre_process保证从0开始）
     *        本块、目标块在本进程中的block号码（绝对值从1开始，负表示周期边界条件）
     *        获取边界面的法向方向
     *        计算转化矩阵
     *        还处理了【tar_block_index this_block_index】，主要在内部传值时使用
     *        -------->针对每一个网格块block中的parallel_bc
     *        通过直接读入获取 sub sup this_myid tar_myid this_block_num
     *        以及send_flag rece_flag
     *        【sub sup以1为起点，direction， Transform未计算】均在Pre_process中处理
     */
    void ReadGeometry(int32_t myid, int32_t ghostmesh, int32_t Method);

    /**
     * @brief 读取plot3D格式文件，grid n.grd
     *        ------->针对每一个网格块block
     *        设置网格名称block_name;
     *        设置网格的尺寸i，j，kmax mx my mz;
     *        设置网格的维度dimension;
     *        设置网格坐标xyz尺寸 jacobi尺寸 metric inver_metric尺寸;
     *        读入坐标xyz;
     */
    void Read_Grid(std::string _my_id_s, int32_t my_id);

    /**
     * @brief 读取plot3D格式二进制文件，grid n.grd
     *        ------->针对每一个网格块block
     *        设置网格名称block_name;
     *        设置网格的尺寸i，j，kmax mx my mz;
     *        设置网格的维度dimension;
     *        设置网格坐标xyz尺寸 jacobi尺寸 metric inver_metric尺寸;
     *        读入坐标xyz;
     */
    void Read_Grid_Binary(std::string _my_id_s, int32_t my_id);

    /**
     * @brief -------->针对每一个网格块block中的physical_bc读入
     *        读入parallel_bc中的起始 终止坐标sub sup;
     *        获取direction boundary_num this_block_num boundary_name;
     */
    void Read_Phy_Boundary(std::string _my_id_s);

    /**
     * @brief -------->针对每一个网格块block中的physical_bc
     *        利用fvbnd获取的boundary_name_index boundary_name;
     */
    void Read_Phy_Name(std::string _my_id_s, int32_t my_id);

    /**
     * @brief -------->针对每一个网格块block中的physical_bc
     *        利用 boundary_name, 根据intList:Boundary_Priority调整优先级
     *        即各个边界条件在grids(izone).physical_bc顺序;
     */
    void Adjust_Phy_Boundary();

    /**
     * @brief -------->针对每一个网格块block中的inner_bc
     *        每读入一个内部边界条件就push_back进inner_bc中
     *        其中包含了本块的边界起始位置，目标块的起始位置（利用Pre_process保证从0开始）
     *        本块、目标块在本进程中的block号码（#绝对值从1开始，负表示周期边界条件#）
     *        获取边界面的法向方向
     *        计算转化矩阵
     *        还处理了【tar_block_index this_block_index】，主要在内部传值时使用
     */
    void Read_Inner_Boundary(std::string _my_id_s, int32_t my_id);

    /**
     * @brief -------->针对每一个网格块block中的parallel_bc
     *        通过直接读入获取 sub sup this_myid tar_myid this_block_num
     *        以及send_flag rece_flag
     *        【sub sup以1为起点，direction， Transform未计算】均在Pre_process中处理
     */
    void Read_Parallel_Boundary(std::string _my_id_s, int32_t my_id);
    //=============================================================================================
    /**
     * @brief 循环所有block，调用grid[i].ghostmesh();
     */
    void Get_ghostmesh()
    {
        for (int32_t i = 0; i < nblock; i++)
            grids(i).Add_ghostmesh(ngg);
    };

    /**
     * @brief 循环所有block，调用grids[i].calc_Jacobi_and_Metrics();
     */
    void Calc_Jacobi_and_Metrics()
    {
        for (int i = 0; i < nblock; i++)
            grids(i).calc_Jacobi_and_Metrics(ngg);

        for (int i = 0; i < nblock; i++)
            grids(i).calc_Dual_Grids();

        for (int i = 0; i < nblock; i++)
            grids(i).calc_Face_Tangent_Vectors(ngg);

        for (int i = 0; i < nblock; i++)
            grids(i).calc_Metrics_GCL(ngg);

        for (int i = 0; i < nblock; i++)
            grids(i).calc_modify_Jacobi(ngg);
    };

public:
    /**
     * @brief 循环所有block，face调用Datatrans_para();
     */
    void MeshTrans_Parallel();

    /**
     * @brief 循环所有block，face调用Datatrans_inner();
     */
    void MeshTrans_Inner();

    /**
     * @brief 循环所有block，给所有的Polar面/轴添加虚网格坐标;
     */
    void MeshTrans_Inner_Pole();
    /**
     * @brief 网格整体放缩
     */
    void scale_mesh(double scale)
    {
        for (int i = 0; i < nblock; i++)
        {
            grids(i).x *= scale;
            grids(i).y *= scale;
            grids(i).z *= scale;
        }
    };

    /**
     * @brief 假虚网格（Faker）区域的几何信息部分处理，分别对Physical, Inner, Parallel的Faker区域开辟空间
     *        物理边界条件直接复制本物理快的虚网格信息；
     *        然后通过调用Faker_Process_MeshTrans_Inner Faker_Process_MeshTrans_Parallel对内部、并行的
     *        物理块耦合界面虚网格进行处理
     */
    void Faker_Process();

    /**
     * @brief 假虚网格区域的内部网格传值处理，直接使用DATATRANS::Inner_flush_scalar _tensor传递几何信息
     */
    void Faker_Process_MeshTrans_Inner();

    /**
     * @brief 假虚网格区域的并行网格传值处理，利用Initial_Allocate_Couple_Grid开辟并行缓冲区，变量数设置为9
     *        分别传递虚网格的几何信息，存入Faker区域。这里9是为了能够一次性传递度量系数
     */
    void Faker_Process_MeshTrans_Parallel();

    /**
     * @brief 针对角区进行预处理，记录哪些边界涉及角区处理
     */
    void Corner_Preprocess()
    {
        for (int i = 0; i < nblock; i++)
            grids(i).calc_Corner_Preprocess();
    }

    /**
     * @brief 针对包含内部边界条件的角区(棱边)进行传值处理
     */
    void MeshTrans_Corner_Inner();

    /**
     * @brief 针对包含内部边界条件的角区(顶点)进行传值处理
     */
    void MeshTrans_Corner3D_Inner();

    /**
     * @brief 针对包含并行边界条件的角区(棱边)进行传值处理
     */
    void MeshTrans_Corner_Parallel();

    /**
     * @brief 针对包含并行边界条件的角区(顶点)进行传值处理
     */
    void MeshTrans_Corner3D_Parallel();
};