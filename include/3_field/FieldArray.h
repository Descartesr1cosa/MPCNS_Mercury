#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <stdint.h>
#include "0_basic/DEFINE.h"

template <typename T>
class FieldScalar
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

    // 从一个已有的3D数组创建
    void SetA3(FieldScalar<T> &a3_in)
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
    FieldScalar() = default;
    // 构造函数，输入各方向尺度和索引开始位置进行构造
    FieldScalar(int32_t n1, int32_t n2, int32_t n3, int ghost)
    {
        dim1 = n1;
        dim2 = n2;
        dim3 = n3;
        ghostmesh = ghost;
        a3.resize(dim1 * dim2 * dim3);
        multip1 = dim2 * dim3;
        addconst = ghostmesh * multip1 + ghostmesh * dim3 + ghostmesh;
    };

    ~FieldScalar() = default;

public:
    T &operator()(int32_t n1, int32_t n2, int32_t n3)
    {
#if if_Debug_Field_Array == 1
        if ((n1 < -ghostmesh || n1 > dim1 - ghostmesh - 1) ||
            (n2 < -ghostmesh || n2 > dim2 - ghostmesh - 1) ||
            (n3 < -ghostmesh || n3 > dim3 - ghostmesh - 1))
        {
            std::cout << "Error ! FieldScalar out of range !!\n";
            exit(-1);
        }
#endif
        return a3[n1 * multip1 + n2 * dim3 + n3 + addconst];
    };

    FieldScalar<T> &operator=(T A1)
    {
        int32_t num = dim1 * dim2 * dim3;
        for (int32_t i = 0; i < num; i++)
            this->a3[i] = A1;
        return *this;
    }

    FieldScalar<T> &operator=(const FieldScalar<T> &A1)
    {
        SetSize(A1.dim1, A1.dim2, A1.dim3, A1.ghostmesh);
        int32_t num = dim1 * dim2 * dim3;
        for (int32_t i = 0; i < num; i++)
            this->a3[i] = A1.a3[i];
        return *this;
    }
};

template <typename T>
class FieldVector
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

    // 从一个已有的FieldVector创建
    void SetV3(FieldVector<T> &v3_in)
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
    FieldVector() = default;
    // 构造函数，输入各方向尺度和索引开始位置进行构造
    FieldVector(int32_t n1, int32_t n2, int32_t n3, int ghost, int32_t vec)
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

    ~FieldVector() = default;

public:
    T &operator()(int32_t n1, int32_t n2, int32_t n3, int32_t n4)
    {
#if if_Debug_Field_Array == 1
        if ((n1 < -ghostmesh || n1 > dim1 - ghostmesh - 1) ||
            (n2 < -ghostmesh || n2 > dim2 - ghostmesh - 1) ||
            (n3 < -ghostmesh || n3 > dim3 - ghostmesh - 1) ||
            (n4 > dim_vec - 1))
        {
            std::cout << "Error ! FieldVector out of range !!\n";
            exit(-1);
        }
#endif
        // return v3[(n1 + ghostmesh) * multip1 + (n2 + ghostmesh) * multip2 + (n3 + ghostmesh) * dim_vec + n4];
        return v3[n1 * multip1 + n2 * multip2 + n3 * dim_vec + n4 + addconst];
    };

    FieldVector<T> &operator=(const FieldVector<T> &A1)
    {
        SetSize(A1.dim1, A1.dim2, A1.dim3, A1.ghostmesh, A1.dim_vec);
        int32_t num = dim1 * dim2 * dim3 * dim_vec;
        for (int32_t i = 0; i < num; i++)
            v3[i] = A1.v3[i];
        return *this;
    }

    FieldVector<T> &operator=(T A1)
    {
        int32_t num = dim1 * dim2 * dim3 * dim_vec;
        for (int32_t i = 0; i < num; i++)
            v3[i] = A1;
        return *this;
    }
};

typedef FieldScalar<double> Scalar;
typedef FieldVector<double> Vector;