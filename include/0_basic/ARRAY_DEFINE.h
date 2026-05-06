/*****************************************************************************
@Copyright: NLCFD
@File name: 0_ARRAY_DEFINE.h
@Author: WangYH Descartes
@Version: 1.0
@Date: 2022年09月11日
@Description:	基于WangYH的程序使用了容器vector定义了可动态分配的多维数组
                重载了多维数组*= /= += -=与其他数的作用
                重载了多维数组与自身相同尺寸的多维数组的+= -=
                Array1D<type> type类型的一维数组
                Array2D<type> type类型的二维数组
                Array3D<type> type类型的三维数组
                （type为任意输入的类型，包括但不限于int、double、string等）
                对括号进行了重载，索引数组某个元素时，采用Value(i,j,k)的方式。
                这里没有添加Array4D 5D 而是根据实际运用选择构造3D空间上的标量（double1D）
                矢量（VectorField）和1-1（0-2 2-0）型张量场（TensorField）
@History:		（修改历史记录列表， 每条修改记录应包括修改日期、修改者及修改内容简述。）
                1、修改了文件名称，最基本的功能头文件名用数字0+全部大写的方式
                2. 添加Array1D operater=重载，防止Array1D = 0.0出错 2024/04/26
                3. 添加Array1D GetA1()引用形式 2024/05/06
*****************************************************************************/

#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <stdint.h>

constexpr double SMALL = 1e-14; // 小量，用于判断=0，或分母中防止NAN

/**
* @brief		自定义一维数组
                起始位置为0，终点位置为dim1-1

* @param dim1	数组大小
* @param a1		存储数据的vector数组
*/
template <typename T>
class Array1D
{
public:
    Array1D() = default;

    Array1D(int32_t n1)
    {
        dim1 = n1;
        a1.resize(dim1);
    };

    ~Array1D() {};

    T &operator()(int32_t n1)
    {
        return a1[n1];
    };

    Array1D<T> &operator=(T A1)
    {
        int32_t num = dim1;
        for (int32_t i = 0; i < num; i++)
            this->a1[i] = A1;
        return *this;
    }
    //-------------------------------------------------------------------------
    template <typename T2>
    Array1D<T> &operator*=(T2 A1)
    {
        for (int32_t i = 0; i < dim1; i++)
            this->operator()(i) = A1 * this->operator()(i);
        return *this;
    };

    template <typename T2>
    Array1D<T> &operator/=(T2 A1)
    {
        if (A1 < SMALL && A1 > -SMALL)
        {
            std::cout << "Warning！除以小量！" << std::endl;
            exit(0);
        }
        for (int32_t i = 0; i < dim1; i++)
            this->operator()(i) = this->operator()(i) / A1;
        return *this;
    };

    template <typename T2>
    Array1D<T> &operator+=(T2 A1)
    {
        for (int32_t i = 0; i < dim1; i++)
            this->operator()(i) = A1 + this->operator()(i);
        return *this;
    };

    template <typename T2>
    Array1D<T> &operator-=(T2 A1)
    {
        for (int32_t i = 0; i < dim1; i++)
            this->operator()(i) = this->operator()(i) - A1;
        return *this;
    };
    //-------------------------------------------------------------------------
    Array1D<T> &operator+=(Array1D<T> &A2)
    {
        if (A2.dim1 != dim1)
        {
            std::cout << "相加的两个数组大小不同" << std::endl;
            exit(0);
        }
        for (int32_t i = 0; i < dim1; i++)
            this->operator()(i) = this->operator()(i) + A2(i);
        return *this;
    };

    Array1D<T> &operator-=(Array1D<T> &A2)
    {
        if (A2.dim1 != dim1)
        {
            std::cout << "相减的两个数组大小不同" << std::endl;
            exit(0);
        }
        for (int32_t i = 0; i < dim1; i++)
            this->operator()(i) = this->operator()(i) - A2(i);
        return *this;
    };
    //-------------------------------------------------------------------------
    std::vector<T> &GetA1() { return a1; };
    void SetA1(std::vector<T> &a1_in, int32_t n1)
    {
        a1 = a1_in;
        dim1 = n1;
    };

    void SetA1(Array1D<T> &a1_in)
    {
        a1 = a1_in.GetA1();
        dim1 = a1_in.Getsize1();
    }

    int32_t Getsize1() { return dim1; };
    void SetSize(int32_t setdim1)
    {
        dim1 = setdim1;
        a1.clear();
        a1.resize(setdim1);
    };

private:
    int32_t dim1;
    std::vector<T> a1;
};

/**
* @brief		自定义二维数组
                索引为0～dim1-1 0~dim2-1
* @remark		数组定义的时候第二维度j为连续存储，循环时最好先第二维度再第一维度提高效率
* @param dim1	数组第一维度大小
* @param dim2	数组第二维度大小
* @param a2		存储数据的vector数组
*/
template <typename T>
class Array2D
{
public:
    Array2D() = default;

    Array2D(int32_t n1, int32_t n2)
    {
        dim1 = n1;
        dim2 = n2;
        a2.resize(dim1 * dim2);
    };

    ~Array2D() {};

    T &operator()(int32_t n1, int32_t n2)
    {
        return a2[n1 * dim2 + n2];
    };
    //-------------------------------------------------------------------------
    template <typename T2>
    Array2D<T> &operator*=(T2 A1)
    {
        for (int32_t i = 0; i < dim1; i++)
            for (int32_t j = 0; j < dim2; j++)
                this->operator()(i, j) = A1 * this->operator()(i, j);
        return *this;
    };

    template <typename T2>
    Array2D<T> &operator/=(T2 A1)
    {
        if (A1 < SMALL && A1 > -SMALL)
        {
            std::cout << "Warning！除以小量！" << std::endl;
            exit(0);
        }
        for (int32_t i = 0; i < dim1; i++)
            for (int32_t j = 0; j < dim2; j++)
                this->operator()(i, j) = this->operator()(i, j) / A1;
        return *this;
    };

    template <typename T2>
    Array2D<T> &operator+=(T2 A1)
    {
        for (int32_t i = 0; i < dim1; i++)
            for (int32_t j = 0; j < dim2; j++)
                this->operator()(i, j) = A1 + this->operator()(i, j);
        return *this;
    };

    template <typename T2>
    Array2D<T> &operator-=(T2 A1)
    {
        for (int32_t i = 0; i < dim1; i++)
            for (int32_t j = 0; j < dim2; j++)
                this->operator()(i, j) = this->operator()(i, j) - A1;
        return *this;
    };

    Array2D<T> &operator+=(Array2D<T> &A2)
    {
        if (A2.dim1 != dim1 || A2.dim2 != dim2)
        {
            std::cout << "相加的两个数组大小不同" << std::endl;
            exit(0);
        }
        for (int32_t i = 0; i < dim1; i++)
            for (int32_t j = 0; j < dim2; j++)
                this->operator()(i, j) = this->operator()(i, j) + A2(i, j);
        return *this;
    };

    Array2D<T> &operator-=(Array2D<T> &A2)
    {
        if (A2.dim1 != dim1 || A2.dim2 != dim2)
        {
            std::cout << "相减的两个数组大小不同" << std::endl;
            exit(0);
        }
        for (int32_t i = 0; i < dim1; i++)
            for (int32_t j = 0; j < dim2; j++)
                this->operator()(i, j) = this->operator()(i, j) - A2(i, j);
        return *this;
    };

    std::vector<T> GetA2() { return a2; };
    void SetA2(std::vector<T> &a2_in, int32_t n1, int32_t n2)
    {
        a2 = a2_in;
        dim1 = n1;
        dim2 = n2;
    };
    void SetA2(Array2D<T> &a2_in)
    {
        a2 = a2_in.GetA2();
        dim1 = a2_in.Getsize1();
        dim2 = a2_in.Getsize2();
    }

    int32_t Getsize1() { return dim1; };
    int32_t Getsize2() { return dim2; };
    void SetSize(int32_t setdim1, int32_t setdim2)
    {
        dim1 = setdim1;
        dim2 = setdim2;
        a2.clear();
        a2.resize(setdim1 * setdim2);
    };

private:
    int32_t dim1, dim2;
    std::vector<T> a2;
};

/**
* @brief		自定义三维数组，三个维度分别索引网格点 i、j、k

* @param dim1	数组第一维度大小
* @param dim2	数组第二维度大小
* @param dim3	数组第三维度大小
* @param ghostmesh	索引的起始位置从-ghostmesh开始，不加改变的时候默认为0
* @param a3		存储数据的vector数组
* @brief 		重载了将一个数用=赋给整体a3 以及+= -= *= /=
                重载了与自身相同尺寸的 += -=
* @remarks 		注意这里的ghostmesh不一定是计算中的ngg，两者应该区分开
* 				e.g. ghostmesh对于网格数组应为ngg+1而流场数组应为ngg
*/
template <typename T>
class Array3D
{
private:
    int32_t dim1 = 0, dim2 = 0, dim3 = 0;
    int32_t ghostmesh = 0;
    int32_t multip1 = 0, addconst = 0;
    std::vector<T> a3;

public:
    int32_t Getsize1() { return dim1; };
    int32_t Getsize2() { return dim2; };
    int32_t Getsize3() { return dim3; };
    int32_t Getghostmesh() { return ghostmesh; };

    std::vector<T> &GetA3() { return a3; };

public:
    void SetGhostmesh(int32_t gmesh) { ghostmesh = gmesh; }

    void SetSize(int32_t setdim1, int32_t setdim2, int32_t setdim3, int32_t ghost)
    {
        dim1 = setdim1;
        dim2 = setdim2;
        dim3 = setdim3;
        a3.resize(setdim1 * setdim2 * setdim3);
        ghostmesh = ghost;
        multip1 = dim2 * dim3;
        addconst = ghostmesh * multip1 + ghostmesh * dim3 + ghostmesh;
    }

    // 从一个已有的1D数组创建，告知各个方向尺度和虚网格层数
    void SetA3(std::vector<T> &a3_in, int32_t n1, int32_t n2, int32_t n3, int ghost)
    {
        a3 = a3_in;
        dim1 = n1;
        dim2 = n2;
        dim3 = n3;
        ghostmesh = ghost;
        multip1 = dim2 * dim3;
        addconst = ghostmesh * multip1 + ghostmesh * dim3 + ghostmesh;
    };

    // 从一个已有的3D数组创建
    void SetA3(Array3D<T> &a3_in)
    {
        a3 = a3_in.GetA3();
        dim1 = a3_in.Getsize1();
        dim2 = a3_in.Getsize2();
        dim3 = a3_in.Getsize3();
        ghostmesh = a3_in.Getghostmesh();
        multip1 = dim2 * dim3;
        addconst = ghostmesh * multip1 + ghostmesh * dim3 + ghostmesh;
    }

public:
    // Array3D() {};
    Array3D() = default;
    // 构造函数，输入各方向尺度和索引开始位置进行构造
    Array3D(int32_t n1, int32_t n2, int32_t n3, int ghost)
    {
        dim1 = n1;
        dim2 = n2;
        dim3 = n3;
        ghostmesh = ghost;
        a3.resize(dim1 * dim2 * dim3);
        multip1 = dim2 * dim3;
        addconst = ghostmesh * multip1 + ghostmesh * dim3 + ghostmesh;
    };

    ~Array3D() = default;

public:
    T &operator()(int32_t n1, int32_t n2, int32_t n3)
    {
        // return a3[(n1 + ghostmesh) * multip1 + (n2 + ghostmesh) * dim3 + n3 + ghostmesh];
        return a3[n1 * multip1 + n2 * dim3 + n3 + addconst];
    };

    Array3D<T> &operator=(T A1)
    {
        int32_t num = dim1 * dim2 * dim3;
        for (int32_t i = 0; i < num; i++)
            this->a3[i] = A1;
        return *this;
    }

    Array3D<T> &operator=(const Array3D<T> &A1)
    {
        SetSize(A1.dim1, A1.dim2, A1.dim3, A1.ghostmesh);
        int32_t num = dim1 * dim2 * dim3;
        for (int32_t i = 0; i < num; i++)
            this->a3[i] = A1.a3[i];
        return *this;
    }

    template <typename T2>
    Array3D<T> &operator*=(T2 A1)
    {
        int32_t num = dim1 * dim2 * dim3;
        for (int32_t i = 0; i < num; i++)
            this->a3[i] = A1 * this->a3[i];
        return *this;
    };

    template <typename T2>
    Array3D<T> &operator/=(T2 A1)
    {
        if (A1 < SMALL && A1 > -SMALL)
        {
            std::cout << "Warning！除以小量！" << std::endl;
            exit(0);
        }
        int32_t num = dim1 * dim2 * dim3;
        for (int32_t i = 0; i < num; i++)
            this->a3[i] = this->a3[i] / A1;
        return *this;
    };

    template <typename T2>
    Array3D<T> &operator+=(T2 A1)
    {
        int32_t num = dim1 * dim2 * dim3;
        for (int32_t i = 0; i < num; i++)
            this->a3[i] = this->a3[i] + A1;
        return *this;
    };

    template <typename T2>
    Array3D<T> &operator-=(T2 A1)
    {
        int32_t num = dim1 * dim2 * dim3;
        for (int32_t i = 0; i < num; i++)
            this->a3[i] = this->a3[i] - A1;
        return *this;
    };
    //-------------------------------------------------------------------------
    Array3D<T> &operator+=(Array3D<T> &A2)
    {
        if (A2.dim1 != dim1 || A2.dim2 != dim2 || A2.dim3 != dim3)
        {
            std::cout << "相加的两个数组大小不同" << std::endl;
            exit(0);
        }
        std::vector<T> &a32 = A2.GetA3();
        int32_t num = dim1 * dim2 * dim3;
        for (int32_t i = 0; i < num; i++)
            this->a3[i] = this->a3[i] + a32[i];
        return *this;
    };

    Array3D<T> &operator-=(Array3D<T> &A2)
    {
        if (A2.dim1 != dim1 || A2.dim2 != dim2 || A2.dim3 != dim3)
        {
            std::cout << "相减的两个数组大小不同" << std::endl;
            exit(0);
        }
        std::vector<T> &a32 = A2.GetA3();
        int32_t num = dim1 * dim2 * dim3;
        for (int32_t i = 0; i < num; i++)
            this->a3[i] = this->a3[i] - a32[i];
        return *this;
    };
};

/**
* @brief			自定义矢量场，三个维度分别索引网格点 i、j、k，存在虚网格，剩下维度为矢量的分量

* @param dim1		数组第一维度大小
* @param dim2		数组第二维度大小
* @param dim3		数组第三维度大小
* @param dim_vec	矢量维度的大小
* @param ghostmesh	索引的起始位置从-ghostmesh开始，不加改变的时候默认为0
* @param v3			存储数据的vector数组
* @brief 			重载了将一个数用=赋给整体v3 以及+= -= *= /=
                    重载了与自身相同尺寸的 += -=
* @remarks 			注意这里的ghostmesh不一定是计算中的ngg，两者应该区分开
* 					e.g. ghostmesh对于网格数组应为ngg+1而流场数组应为ngg
*					矢量维度不用虚网格
*/
template <typename T>
class VectorField
{
private:
    int32_t dim1 = 0, dim2 = 0, dim3 = 0;
    int32_t dim_vec = 0;
    int32_t ghostmesh = 0;
    int32_t multip1 = 0, multip2 = 0, addconst = 0;
    std::vector<T> v3;

public:
    int32_t Getsize1() { return dim1; };
    int32_t Getsize2() { return dim2; };
    int32_t Getsize3() { return dim3; };
    int32_t Getsizevec() { return dim_vec; };
    int32_t Getghostmesh() { return ghostmesh; };

    std::vector<T> &GetV3() { return v3; };

public:
    void SetGhostmesh(int32_t gmesh) { ghostmesh = gmesh; }

    void SetSize(int32_t setdim1, int32_t setdim2, int32_t setdim3, int32_t ghost, int32_t vec_length)
    {
        dim1 = setdim1;
        dim2 = setdim2;
        dim3 = setdim3;
        dim_vec = vec_length;
        v3.resize(setdim1 * setdim2 * setdim3 * vec_length);
        ghostmesh = ghost;
        multip1 = dim2 * dim3 * dim_vec;
        multip2 = dim3 * dim_vec;
        addconst = ghostmesh * (multip1 + multip2 + dim_vec);
    }

    // 从一个已有的1D数组创建，告知各个方向尺度和虚网格层数
    void SetV3(std::vector<T> &v3_in, int32_t n1, int32_t n2, int32_t n3, int32_t ghost, int32_t vec_length)
    {
        v3 = v3_in;
        dim1 = n1;
        dim2 = n2;
        dim3 = n3;
        dim_vec = vec_length;
        ghostmesh = ghost;
        multip1 = dim2 * dim3 * dim_vec;
        multip2 = dim3 * dim_vec;
        addconst = ghostmesh * (multip1 + multip2 + dim_vec);
    };

    // 从一个已有的VectorField创建
    void SetV3(VectorField<T> &v3_in)
    {
        v3 = v3_in.GetV3();
        dim1 = v3_in.Getsize1();
        dim2 = v3_in.Getsize2();
        dim3 = v3_in.Getsize3();
        dim_vec = v3_in.Getsizevec();
        ghostmesh = v3_in.Getghostmesh();
        multip1 = dim2 * dim3 * dim_vec;
        multip2 = dim3 * dim_vec;
        addconst = ghostmesh * (multip1 + multip2 + dim_vec);
    }

public:
    VectorField() = default;
    // 构造函数，输入各方向尺度和索引开始位置进行构造
    VectorField(int32_t n1, int32_t n2, int32_t n3, int ghost, int32_t vec)
    {
        dim1 = n1;
        dim2 = n2;
        dim3 = n3;
        dim_vec = vec;
        ghostmesh = ghost;
        v3.resize(dim1 * dim2 * dim3 * dim_vec);
        multip1 = dim2 * dim3 * dim_vec;
        multip2 = dim3 * dim_vec;
        addconst = ghostmesh * (multip1 + multip2 + dim_vec);
    };

    ~VectorField() = default;

public:
    T &operator()(int32_t n1, int32_t n2, int32_t n3, int32_t n4)
    {
        // return v3[(n1 + ghostmesh) * multip1 + (n2 + ghostmesh) * multip2 + (n3 + ghostmesh) * dim_vec + n4];
        return v3[n1 * multip1 + n2 * multip2 + n3 * dim_vec + n4 + addconst];
    };

    VectorField<T> &operator=(const VectorField<T> &A1)
    {
        SetSize(A1.dim1, A1.dim2, A1.dim3, A1.ghostmesh, A1.dim_vec);
        int32_t num = dim1 * dim2 * dim3 * dim_vec;
        for (int32_t i = 0; i < num; i++)
            v3[i] = A1.v3[i];
        return *this;
    }

    VectorField<T> &operator=(T A1)
    {
        int32_t num = dim1 * dim2 * dim3 * dim_vec;
        for (int32_t i = 0; i < num; i++)
            v3[i] = A1;
        return *this;
    }

    template <typename T2>
    VectorField<T> &operator*=(T2 A1)
    {
        int32_t num = dim1 * dim2 * dim3 * dim_vec;
        for (int32_t i = 0; i < num; i++)
            this->v3[i] = A1 * this->v3[i];
        return *this;
    };

    template <typename T2>
    VectorField<T> &operator/=(T2 A1)
    {
        if (A1 < SMALL && A1 > -SMALL)
        {
            std::cout << "Warning！除以小量！" << std::endl;
            exit(0);
        }
        int32_t num = dim1 * dim2 * dim3 * dim_vec;
        for (int32_t i = 0; i < num; i++)
            this->v3[i] = this->v3[i] / A1;
        return *this;
    };

    template <typename T2>
    VectorField<T> &operator+=(T2 A1)
    {
        int32_t num = dim1 * dim2 * dim3 * dim_vec;
        for (int32_t i = 0; i < num; i++)
            this->v3[i] = this->v3[i] + A1;
        return *this;
    };

    template <typename T2>
    VectorField<T> &operator-=(T2 A1)
    {
        int32_t num = dim1 * dim2 * dim3 * dim_vec;
        for (int32_t i = 0; i < num; i++)
            this->v3[i] = this->v3[i] - A1;
        return *this;
    };
    //-------------------------------------------------------------------------
    VectorField<T> &operator+=(VectorField<T> &A2)
    {
        if (A2.dim1 != dim1 || A2.dim2 != dim2 || A2.dim3 != dim3 || A2.dim_vec != dim_vec)
        {
            std::cout << "相加的两个数组大小不同" << std::endl;
            exit(0);
        }
        std::vector<T> &v32 = A2.GetV3();
        int32_t num = dim1 * dim2 * dim3 * dim_vec;
        for (int32_t i = 0; i < num; i++)
            this->v3[i] = this->v3[i] + v32[i];
        return *this;
    };

    VectorField<T> &operator-=(VectorField<T> &A2)
    {
        if (A2.dim1 != dim1 || A2.dim2 != dim2 || A2.dim3 != dim3 || A2.dim_vec != dim_vec)
        {
            std::cout << "相减的两个数组大小不同" << std::endl;
            exit(0);
        }
        std::vector<T> &v32 = A2.GetV3();
        int32_t num = dim1 * dim2 * dim3 * dim_vec;
        for (int32_t i = 0; i < num; i++)
            this->v3[i] = this->v3[i] - v32[i];
        return *this;
    };
};

/**
* @brief			自定义张量场，三个维度分别索引网格点 i、j、k，存在虚网格，剩下维度为张量的分量

* @param dim1		数组第一维度大小
* @param dim2		数组第二维度大小
* @param dim3		数组第三维度大小
* @param dim_ten1	张量第一维度的大小
* @param dim_ten2	张量第二维度的大小
* @param ghostmesh	索引的起始位置从-ghostmesh开始，不加改变的时候默认为0
* @param t3			存储数据的vector数组
* @brief 			重载了将一个数用=赋给整体t3 以及+= -= *= /=
                    重载了与自身相同尺寸的 += -=
* @remarks 			注意这里的ghostmesh不一定是计算中的ngg，两者应该区分开
* 					e.g. ghostmesh对于网格数组应为ngg+1而流场数组应为ngg
*					张量维度不用虚网格
*/
template <typename T>
class TensorField
{
private:
    int32_t dim1 = 0, dim2 = 0, dim3 = 0;
    int32_t dim_ten1 = 0, dim_ten2 = 0;
    int32_t ghostmesh = 0;
    int32_t multip1 = 0, multip2 = 0, multip3 = 0, addconst = 0;
    std::vector<T> t3;

public:
    int32_t Getsize1() { return dim1; };
    int32_t Getsize2() { return dim2; };
    int32_t Getsize3() { return dim3; };
    int32_t Getsizeten1() { return dim_ten1; };
    int32_t Getsizeten2() { return dim_ten2; };
    int32_t Getghostmesh() { return ghostmesh; };

    std::vector<T> &GetT3() { return t3; };

public:
    void SetGhostmesh(int32_t gmesh) { ghostmesh = gmesh; }

    void SetSize(int32_t setdim1, int32_t setdim2, int32_t setdim3, int32_t ghost, int32_t ten_length1, int32_t ten_length2)
    {
        dim1 = setdim1;
        dim2 = setdim2;
        dim3 = setdim3;
        dim_ten1 = ten_length1;
        dim_ten2 = ten_length2;
        t3.resize(setdim1 * setdim2 * setdim3 * ten_length1 * ten_length2);
        ghostmesh = ghost;
        multip1 = dim2 * dim3 * dim_ten1 * dim_ten2;
        multip2 = dim3 * dim_ten1 * dim_ten2;
        multip3 = dim_ten1 * dim_ten2;
        addconst = ghostmesh * (multip1 + multip2 + multip3);
    }

    // 从一个已有的1D数组创建，告知各个方向尺度和虚网格层数
    void SetT3(std::vector<T> &t3_in, int32_t n1, int32_t n2, int32_t n3, int32_t ghost, int32_t n4, int32_t n5)
    {
        t3 = t3_in;
        dim1 = n1;
        dim2 = n2;
        dim3 = n3;
        dim_ten1 = n4;
        dim_ten2 = n5;
        ghostmesh = ghost;
        multip1 = dim2 * dim3 * dim_ten1 * dim_ten2;
        multip2 = dim3 * dim_ten1 * dim_ten2;
        multip3 = dim_ten1 * dim_ten2;
        addconst = ghostmesh * (multip1 + multip2 + multip3);
    };

    // 从一个已有的TensorField创建
    void SetT3(TensorField<T> &t3_in)
    {
        t3 = t3_in.GetT3();
        dim1 = t3_in.Getsize1();
        dim2 = t3_in.Getsize2();
        dim3 = t3_in.Getsize3();
        dim_ten1 = t3_in.Getsizeten1();
        dim_ten2 = t3_in.Getsizeten2();
        ghostmesh = t3_in.Getghostmesh();
        multip1 = dim2 * dim3 * dim_ten1 * dim_ten2;
        multip2 = dim3 * dim_ten1 * dim_ten2;
        multip3 = dim_ten1 * dim_ten2;
        addconst = ghostmesh * (multip1 + multip2 + multip3);
    }

public:
    TensorField() = default;
    // 构造函数，输入各方向尺度和索引开始位置进行构造
    TensorField(int32_t n1, int32_t n2, int32_t n3, int ghost, int32_t ten1, int32_t ten2)
    {
        dim1 = n1;
        dim2 = n2;
        dim3 = n3;
        dim_ten1 = ten1;
        dim_ten2 = ten2;
        ghostmesh = ghost;
        t3.resize(dim1 * dim2 * dim3 * dim_ten1 * dim_ten2);
        multip1 = dim2 * dim3 * dim_ten1 * dim_ten2;
        multip2 = dim3 * dim_ten1 * dim_ten2;
        multip3 = dim_ten1 * dim_ten2;
        addconst = ghostmesh * (multip1 + multip2 + multip3);
    };

    ~TensorField() = default;

public:
    T &operator()(int32_t n1, int32_t n2, int32_t n3, int32_t n4, int32_t n5)
    {
        // return t3[(n1 + ghostmesh) * dim2 * dim3 * dim_ten1 * dim_ten2 + (n2 + ghostmesh) * dim3 * dim_ten1 * dim_ten2 + (n3 + ghostmesh) * dim_ten1 * dim_ten2 + n4 * dim_ten2 + n5];
        return t3[n1 * multip1 + n2 * multip2 + n3 * multip3 + n4 * dim_ten2 + n5 + addconst];
    };

    TensorField<T> &operator=(T T1)
    {
        int32_t num = dim1 * dim2 * dim3 * dim_ten1 * dim_ten2;
        for (int32_t i = 0; i < num; i++)
            this->t3[i] = T1;
        return *this;
    }

    TensorField<T> &operator=(const TensorField<T> &T1)
    {
        SetSize(T1.dim1, T1.dim2, T1.dim3, T1.ghostmesh, T1.dim_ten1, T1.dim_ten2);
        int32_t num = dim1 * dim2 * dim3 * dim_ten1 * dim_ten2;
        for (int32_t i = 0; i < num; i++)
            t3[i] = T1.t3[i];
        return *this;
    }

    template <typename T2>
    TensorField<T> &operator*=(T2 T1)
    {
        int32_t num = dim1 * dim2 * dim3 * dim_ten1 * dim_ten2;
        for (int32_t i = 0; i < num; i++)
            this->t3[i] = T1 * this->t3[i];
        return *this;
    };

    template <typename T2>
    TensorField<T> &operator/=(T2 T1)
    {
        if (T1 < SMALL && T1 > -SMALL)
        {
            std::cout << "Warning！除以小量！" << std::endl;
            exit(0);
        }
        int32_t num = dim1 * dim2 * dim3 * dim_ten1 * dim_ten2;
        for (int32_t i = 0; i < num; i++)
            this->t3[i] = this->t3[i] / T1;
        return *this;
    };

    template <typename T2>
    TensorField<T> &operator+=(T2 T1)
    {
        int32_t num = dim1 * dim2 * dim3 * dim_ten1 * dim_ten2;
        for (int32_t i = 0; i < num; i++)
            this->t3[i] = this->t3[i] + T1;
        return *this;
    };

    template <typename T2>
    TensorField<T> &operator-=(T2 T1)
    {
        int32_t num = dim1 * dim2 * dim3 * dim_ten1 * dim_ten2;
        for (int32_t i = 0; i < num; i++)
            this->t3[i] = this->t3[i] - T1;
        return *this;
    };
    //-------------------------------------------------------------------------
    TensorField<T> &operator+=(TensorField<T> &T2)
    {
        if (T2.dim1 != dim1 || T2.dim2 != dim2 || T2.dim3 != dim3 || T2.dim_ten1 != dim_ten1 || T2.dim_ten2 != dim_ten2)
        {
            std::cout << "相加的两个数组大小不同" << std::endl;
            exit(0);
        }
        std::vector<T> &v32 = T2.GetT3();
        int32_t num = dim1 * dim2 * dim3 * dim_ten1 * dim_ten2;
        for (int32_t i = 0; i < num; i++)
            this->t3[i] = this->t3[i] + v32[i];
        return *this;
    };

    TensorField<T> &operator-=(TensorField<T> &T2)
    {
        if (T2.dim1 != dim1 || T2.dim2 != dim2 || T2.dim3 != dim3 || T2.dim_ten1 != dim_ten1 || T2.dim_ten2 != dim_ten2)
        {
            std::cout << "相减的两个数组大小不同" << std::endl;
            exit(0);
        }
        std::vector<T> &v32 = T2.GetT3();
        int32_t num = dim1 * dim2 * dim3 * dim_ten1 * dim_ten2;
        for (int32_t i = 0; i < num; i++)
            this->t3[i] = this->t3[i] - v32[i];
        return *this;
    };
};

// typedef::
typedef Array1D<int> int1D;
typedef Array1D<double> double1D;
typedef Array1D<std::string> string1D;
typedef Array1D<char> char1D;

typedef Array2D<int> int2D;
typedef Array2D<double> double2D;
typedef Array2D<std::string> string2D;
typedef Array2D<char> char2D;

typedef Array3D<int> int3D;
typedef Array3D<double> double3D;
typedef Array3D<std::string> string3D;
typedef Array3D<char> char3D;

typedef VectorField<double> Phy_Vector;
typedef TensorField<double> Phy_Tensor;
