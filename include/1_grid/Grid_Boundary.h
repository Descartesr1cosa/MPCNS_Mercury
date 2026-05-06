/*****************************************************************************
@Copyright: NLCFD
@File name: 1_MPCNS_Boundary.h
@Author: Descartes
@Version: 1.0
@Date: 2023年09月05日
@Description:	基于自定义的可变长度Array数组定义了实现边界条件的数据结构，用于MPCNS计算
@History:		（修改历史记录列表， 每条修改记录应包括修改日期、修改者及修改内容简述。）
                1、修改了文件名称，基本的功能头文件名用数字1+全部首字母大写的方式
                2、2022/10/04/16：52修改了Inner_Boundary的this_block_num,tar_block_num
                他们均从1开始，若为负则表示周期边界条件
                3、2023/09/05 21：28添加了Faker_Boundary边界条件，用于多物理场面耦合假
                虚网格，也可添加来流速度型、吹吸边界条件等。
                4、2023/04/26 将cycle加入便于调用，表示外法线方向；修改this_block_num，tar_block_num均为正从0开始
                5、2024/04/27 将is_corner加入便于判断角区，Para中还有is_corner_send表示是否需要发送消息
*****************************************************************************/
#pragma once
#include "0_basic/ARRAY_DEFINE.h"
#include <map>

/**
 * @brief 定义了用于存储多物理场面耦合缓冲区（假虚网格）的类，用于存储和传递数据，也可用于物理边界条件添加，注意内存空间的控制
 * @param    FakerField      假虚网格的物理场，vec_num长度由该耦合面所对应物理场的物理量数量决定，坐标与sub sup范围相同
 * @param    FakerPV         假虚网格的临时原始变量场，作为原始变量用于添加耦合/物理面的边界条件
 * @param    Fakerx          假虚网格的坐标x
 * @param    Fakery          假虚网格的坐标y
 * @param    Fakerz          假虚网格的坐标z
 * @param    Fakerjac        假虚网格的Jacobian
 * @param    Fakermetric     假虚网格的度规（度量系数）
 * @param    is_multi_phys   是否为多物理场面耦合边界条件
 */
class Faker_Boundary
{
public:
    std::map<std::string, Phy_Vector> FakerField;
    Phy_Vector FakerPV;
    double3D Fakerx, Fakery, Fakerz, Fakerjac;
    Phy_Tensor Fakermetric;
    bool is_multi_phys;

public:
    Faker_Boundary() = default;
    ~Faker_Boundary() = default;
    Faker_Boundary &operator=(const Faker_Boundary &fake)
    {
        FakerPV = fake.FakerPV;
        Fakerx = fake.Fakerx;
        Fakery = fake.Fakery;
        Fakerz = fake.Fakerz;
        Fakerjac = fake.Fakerjac;
        Fakermetric = fake.Fakermetric;
        is_multi_phys = fake.is_multi_phys;

        if (fake.FakerField.size() != 0)
            for (const auto &pair : fake.FakerField)
            {
                Phy_Vector tmp = pair.second;
                std::string name = pair.first;
                FakerField[name] = tmp;
            }
        return *this;
    }
};

/**
 * @brief 定义了并行边界条件的类，用于存储和MPI并行传递数据// 注意：以下所有说明，均指调用 Pre_process(dimension) 之后的内部表示。
 * @param    sub             边界块的i,j,k坐标范围一: 0-based 正负皆可
 * @param    sup             边界块的i,j,k坐标范围二: 0-based 正负皆可
 * @param    cycle           该边界网格i,j,k的法向,I,j,k+cycle[]指向网格外法向：e.g. 0 ±1 0
 * @param    direction       边界法向的方向，正负表示面的大小号如1表示i的大号面，-1为小号面: 1-based
 * @param    tar_myid        目标边界块进程号: 0-based
 * @param    this_myid       本块的进程号: 0-based
 * @param    this_block_num  本块的块序号: 0-based
 * @param    send_flag       发送消息编号
 * @param    rece_flag       接收消息编号
 * @param    Transform       发送时本块按照k j i时数组排布顺序信息以1表示坐标相同的方向（法），以2表示为负的方向， 以0表示剩下的方向
 * @param    target_block_name 目标块名称
 * @param    this_block_name 本块名称
 * @param    Faker           假虚网格
 */
class Parallel_Boundary
{
public:
    int32_t sub[3], sup[3], cycle[3];
    int32_t direction;
    int32_t this_myid, tar_myid;
    int32_t this_block_num;
    int32_t send_flag, rece_flag;
    int32_t Transform[3];
    std::string target_block_name, this_block_name;
    Faker_Boundary Faker;

public:
    Parallel_Boundary() {};
    ~Parallel_Boundary() {};

    void Pre_process(int32_t &dimension)
    {
        // 绝对数值减1,保证编号从零开始
        minus1(sub);
        minus1(sup);

        // 设置边界的法向，123表示方向，+-表示大小号面
        for (int32_t i = 0; i < 3; i++)
        {
            cycle[i] = 0;
            if (i == 2 && dimension == 2)
                continue;
            if (sub[i] == sup[i])
            {
                direction = (sub[i] == 0) ? -(i + 1) : (i + 1);
                cycle[i] = (sub[i] == 0) ? -1 : 1;
            }
        }

        // 计算Transform，以1表示坐标相同的方向（法），以2表示为负的方向， 以0表示剩下的方向
        // 均从sub到sup,因此不存在正负问题
        //----------------------------------------
        // 对于三维问题
        if (dimension == 3)
        {
            if (sup[2] == sub[2])
            {
                Transform[2] = 1;
                if (sup[1] < 0 || sub[1] < 0)
                {
                    Transform[1] = 2;
                    Transform[0] = 0;
                }
                else
                {
                    Transform[1] = 0;
                    Transform[0] = 2;
                }
            }
            else if (sup[2] < 0 || sub[2] < 0)
            {
                Transform[2] = 2;
                if (sup[1] == sub[1])
                {
                    Transform[1] = 1;
                    Transform[0] = 0;
                }
                else
                {
                    Transform[1] = 0;
                    Transform[0] = 1;
                }
            }
            else
            {
                Transform[2] = 0;
                if (sup[1] < 0 || sub[1] < 0)
                {
                    Transform[1] = 2;
                    Transform[0] = 1;
                }
                else
                {
                    Transform[1] = 1;
                    Transform[0] = 2;
                }
            }
        }
        // 对于二维问题 找为相等的方向为1,剩下为2
        else if (dimension == 2)
        {
            Transform[2] = 0;
            // 二维的情况下不用负数表示，用是否相等即可判断
            if (sup[1] == sub[1])
            {
                Transform[0] = 2;
                Transform[1] = 1;
            }
            else
            {
                Transform[1] = 2;
                Transform[0] = 1;
            }
        }
    }

private:
    void minus1(int32_t *p)
    {
        // p[0]p[1]p[2] 绝对数值减1
        for (int i = 0; i < 3; i++)
        {
            if (p[i] > 0)
            {
                p[i] = p[i] - 1;
            }
            else
            {
                p[i] = p[i] + 1;
            }
        }
    }
};

/**
 * @brief 定义了内部多块边界条件的类，用于存储和同进程传递数据（赋值）// 注意：以下所有说明，均指调用 Pre_process(dimension) 之后的内部表示。
 * @param    sub             边界块的i,j,k坐标范围一: 0-based 正负皆可
 * @param    sup             边界块的i,j,k坐标范围二: 0-based 正负皆可
 * @param    cycle           该边界网格i,j,k的法向,I,j,k+cycle[]指向网格外法向：e.g. 0 ±1 0
 * @param    tar_cycle       目标边界网格i,j,k的法向,I,j,k+cycle[]指向网格外法向：e.g. 0 ±1 0
 * @param    tar_sub         所对应的目标边界块的i,j,k坐标范围一: 0-based 正负皆可
 * @param    tar_sup         所对应的目标边界块的i,j,k坐标范围二: 0-based 正负皆可
 * @param    Transform       发送时本块按照k j i时数组排布顺序信息
 * @param    tar_Transform   目标块按照k j i时数组排布顺序信息
 * @param    direction       边界法向的方向，正负表示面的大小号如1表示i的大号面，-1为小号面: 1-based
 * @param    tar_direction   目标边界块边界法向的方向，正负表示面的大小号如1表示i的大号面，-1为小号面: 1-based
 * @param    tar_block_num   目标边界块在其进程中块的序号：0-based，一定为正，周期边界由is_period判断
 * @param    this_block_num  本块的块序号：0-based，一定为正，周期边界由is_period判断
 * @param    tar_block_index   目标边界在inner_bound数组中的序号//需要在grids中处理
 * @param    this_block_index  本边界的在inner_bound数组中的序号//需要在grids中处理
 * @param    index           用来标记本内部边界条件是否已经传值，未传值初始化为0（e.g.于grids中MeshTrans_Inner()初始化）
 * @param    target_block_name 目标块名称
 * @param    this_block_name 本块名称
 * @param    is_period       是否为周期边界
 * @param    Faker           假虚网格
 */
class Inner_Boundary
{
public:
    int32_t sub[3], sup[3], cycle[3], tar_cycle[3];
    int32_t tar_sub[3], tar_sup[3], Transform[3], tar_Transform[3];
    int32_t direction, tar_direction;
    int32_t this_block_num, tar_block_num;
    int32_t this_block_index, tar_block_index;
    int32_t index;
    std::string target_block_name, this_block_name;
    Faker_Boundary Faker;
    bool is_period;

public:
    Inner_Boundary()
    {
        tar_block_index = -1;
        this_block_index = -1;
    };
    ~Inner_Boundary() {};

    void Pre_process(int32_t &dimension)
    {
        // 绝对数值减1,保证编号从零开始
        minus1(sub);
        minus1(sup);
        minus1(tar_sub);
        minus1(tar_sup);

        // 设置边界的法向，123表示方向，+-表示大小号面
        for (int32_t i = 0; i < 3; i++)
        {
            cycle[i] = 0;
            tar_cycle[i] = 0;
            if (i == 2 && dimension == 2)
                continue;
            if (sub[i] == sup[i])
            {
                direction = (sub[i] == 0) ? -(i + 1) : (i + 1);
                cycle[i] = (sub[i] == 0) ? -1 : 1;
            }
        }
        for (int32_t i = 0; i < 3; i++)
        {
            if (i == 2 && dimension == 2)
                continue;
            if (tar_sub[i] == tar_sup[i])
            {
                tar_direction = (tar_sub[i] == 0) ? -(i + 1) : (i + 1);
                tar_cycle[i] = (tar_sub[i] == 0) ? -1 : 1;
            }
        }
        //===========================================================
        // 计算转化矩阵，transform适用于求转化矩阵的，transform[]表示本块ijk对应目标块+-i、+-j、+-k
        // 额外需要一个Transorm
        // 计算Transform，以1表示坐标相同的方向（法），以2表示为负的方向， 以0表示剩下的方向
        // 初始化为-1
        for (int32_t i = 0; i < dimension; i++)
            Transform[i] = -1;

        for (int32_t i = 0; i < dimension; i++)
        {
            if (sub[i] == sup[i])
            {
                Transform[i] = 1;
            }
            else if (sub[i] < 0 || sup[i] < 0)
            {
                Transform[i] = 2;
            }
        }
        if (dimension == 3)
        {
            for (int32_t i = 0; i < 3; i++)
            {
                if (Transform[i] == -1)
                {
                    Transform[i] = 0;
                }
            }
        }
        else
        {
            Transform[2] = 0;
            for (int32_t i = 0; i < 3; i++)
            {
                if (Transform[i] == -1)
                {
                    Transform[i] = 2;
                }
            }
        }
        //-------------------------------------------------------------
        // 处理目标块的转化关系
        for (int32_t i = 0; i < dimension; i++)
            tar_Transform[i] = -1;

        for (int32_t i = 0; i < dimension; i++)
        {
            if (tar_sub[i] == tar_sup[i])
            {
                tar_Transform[i] = 1;
            }
            else if (tar_sub[i] < 0 || tar_sup[i] < 0)
            {
                tar_Transform[i] = 2;
            }
        }
        if (dimension == 3)
        {
            for (int32_t i = 0; i < 3; i++)
            {
                if (tar_Transform[i] == -1)
                {
                    tar_Transform[i] = 0;
                }
            }
        }
        else
        {
            tar_Transform[2] = 0;
            for (int32_t i = 0; i < 3; i++)
            {
                if (tar_Transform[i] == -1)
                {
                    tar_Transform[i] = 2;
                }
            }
        }
    }

private:
    void minus1(int32_t *p)
    {
        // p[0]p[1]p[2] 绝对数值减1
        for (int i = 0; i < 3; i++)
        {
            if (p[i] > 0)
            {
                p[i] = p[i] - 1;
            }
            else
            {
                p[i] = p[i] + 1;
            }
        }
    }
};

/**
 * @brief 定义了物理边界条件的类，用于存储
 * @param    sub             边界块的i,j,k坐标下界: 0-based 正值 sub<=sup
 * @param    sup             边界块的i,j,k坐标上界: 0-based 正值 sup>=sub
 * @param    cycle           该边界网格i,j,k的法向,I,j,k+cycle[]指向网格外法向：e.g. 0 ±1 0
 * @param    direction       边界法向的方向，正负表示面的大小号如1表示i的大号面，-1为小号面: 1-based
 * @param    this_block_num  本块的块序号：0-based，一定为正
 * @param    boundary_num    物理边界条件的序号：0-based，一定为正
 * @param    boundary_name   物理边界条件的名称
 * @param    Faker           假虚网格
 */
class Physical_Boundary
{
public:
    int32_t sub[3], sup[3], cycle[3];
    int32_t direction;
    int32_t this_block_num;
    int32_t boundary_num;
    std::string boundary_name;
    Faker_Boundary Faker;
    std::string target_block_name, this_block_name; // 仅仅用于统一三种边界条件类型

public:
    Physical_Boundary() {};
    ~Physical_Boundary() {};

    Physical_Boundary(const Physical_Boundary &phy)
    {
        for (int i = 0; i < 3; i++)
        {
            sub[i] = phy.sub[i];
            sup[i] = phy.sup[i];
            cycle[i] = phy.cycle[i];
        }
        direction = phy.direction;
        this_block_num = phy.this_block_num;
        boundary_num = phy.boundary_num;
        boundary_name = phy.boundary_name;
        target_block_name = phy.target_block_name;
        this_block_name = phy.this_block_name;
        Faker = phy.Faker;
    }
    Physical_Boundary &operator=(const Physical_Boundary &phy)
    {
        for (int i = 0; i < 3; i++)
        {
            sub[i] = phy.sub[i];
            sup[i] = phy.sup[i];
            cycle[i] = phy.cycle[i];
        }
        direction = phy.direction;
        this_block_num = phy.this_block_num;
        boundary_num = phy.boundary_num;
        boundary_name = phy.boundary_name;
        target_block_name = phy.target_block_name;
        this_block_name = phy.this_block_name;
        Faker = phy.Faker;

        return *this;
    };
};

class Edge
{
public:
    int32_t sub[3], sup[3], cycle1[3], cycle2[3];
    int32_t direction1, direction2, this_block_num;
    int32_t index; // 0-inner 1-para 2-phy
    bool is_singular, if_large_singular;
    int32_t sub_mid[3], sup_mid[3]; // 仅用于para传值

    Inner_Boundary *inner_bound;
    Parallel_Boundary *para_bound;
    Physical_Boundary *phy_bound;

    Edge()
    {
        for (int i = 0; i < 3; i++)
        {
            sub[i] = 0;
            sup[i] = 0;
            sub_mid[i] = 0;
            sup_mid[i] = 0;
            cycle1[i] = 0;
            cycle2[i] = 0;
        }
        direction1 = 0;
        direction2 = 0;
        this_block_num = 0;
        index = 0;
        is_singular = 0;
        inner_bound = nullptr;
        para_bound = nullptr;
        phy_bound = nullptr;
    };
    Edge(const Edge &e)
    {
        for (int i = 0; i < 3; i++)
        {
            sub[i] = e.sub[i];
            sup[i] = e.sup[i];
            sub_mid[i] = e.sub_mid[i];
            sup_mid[i] = e.sup_mid[i];
            cycle1[i] = e.cycle1[i];
            cycle2[i] = e.cycle2[i];
        }
        direction1 = e.direction1;
        direction2 = e.direction2;
        this_block_num = e.this_block_num;
        index = e.index;
        is_singular = e.is_singular;
        inner_bound = e.inner_bound;
        para_bound = e.para_bound;
        phy_bound = e.phy_bound;
    }
    Edge &operator=(const Edge &e)
    {
        for (int i = 0; i < 3; i++)
        {
            sub[i] = e.sub[i];
            sup[i] = e.sup[i];
            sub_mid[i] = e.sub_mid[i];
            sup_mid[i] = e.sup_mid[i];

            cycle1[i] = e.cycle1[i];
            cycle2[i] = e.cycle2[i];
        }
        direction1 = e.direction1;
        direction2 = e.direction2;
        this_block_num = e.this_block_num;
        index = e.index;
        is_singular = e.is_singular;
        inner_bound = e.inner_bound;
        para_bound = e.para_bound;
        phy_bound = e.phy_bound;

        return *this;
    };
    ~Edge() {};
};