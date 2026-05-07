#pragma once
#include <cstddef>
#include <cmath>
#include <stdexcept>

class Vec3
{
public:
    double vec[3]; // 你要求的内部数据结构：vec[3]

    Vec3()
    {
        vec[0] = 0.0;
        vec[1] = 0.0;
        vec[2] = 0.0;
    };
    Vec3(double x, double y, double z)
    {
        vec[0] = x;
        vec[1] = y;
        vec[2] = z;
    };

    // 直接访问
    double &operator[](std::size_t i)
    {
        if (i > 2)
            throw std::out_of_range("Vec3 index out of range");
        return vec[i];
    };
    const double &operator[](std::size_t i) const
    {
        if (i > 2)
            throw std::out_of_range("Vec3 index out of range");
        return vec[i];
    };

    // 像 double 一样的算术
    Vec3 operator+() const { return *this; };
    Vec3 operator-() const { return Vec3(-vec[0], -vec[1], -vec[2]); }

    Vec3 &operator+=(const Vec3 &rhs)
    {
        vec[0] += rhs.vec[0];
        vec[1] += rhs.vec[1];
        vec[2] += rhs.vec[2];
        return *this;
    };
    Vec3 &operator-=(const Vec3 &rhs)
    {
        vec[0] -= rhs.vec[0];
        vec[1] -= rhs.vec[1];
        vec[2] -= rhs.vec[2];
        return *this;
    };
    Vec3 &operator*=(double s)
    {
        vec[0] *= s;
        vec[1] *= s;
        vec[2] *= s;
        return *this;
    };
    Vec3 &operator/=(double s)
    {
        vec[0] /= s;
        vec[1] /= s;
        vec[2] /= s;
        return *this;
    };

    // 点积：v * w -> double（按你想要的 * 作为内积）
    double operator*(const Vec3 &rhs) const
    {
        return vec[0] * rhs.vec[0] + vec[1] * rhs.vec[1] + vec[2] * rhs.vec[2];
    };

    // 叉积：v ^ w -> Vec3（按你想要的 ^ 作为叉积）
    Vec3 operator^(const Vec3 &rhs) const
    {
        return Vec3(
            vec[1] * rhs.vec[2] - vec[2] * rhs.vec[1],
            vec[2] * rhs.vec[0] - vec[0] * rhs.vec[2],
            vec[0] * rhs.vec[1] - vec[1] * rhs.vec[0]);
    };

    // 范数
    double norm2() const { return (*this) * (*this); };
    double norm() const { return std::sqrt(norm2()); };
};

// 非成员二元运算（用于支持 Vec3 + Vec3 等）
inline Vec3 operator+(Vec3 lhs, const Vec3 &rhs) { return lhs += rhs; };
inline Vec3 operator-(Vec3 lhs, const Vec3 &rhs) { return lhs -= rhs; };

// 标量乘除（用于支持 Vec3 * double, double * Vec3）
inline Vec3 operator*(Vec3 v, double s) { return v *= s; };
inline Vec3 operator*(double s, Vec3 v) { return v *= s; };
inline Vec3 operator/(Vec3 v, double s) { return v /= s; };
