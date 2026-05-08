#pragma once
#include <cstdio>
#include <cstdlib>

// 用于缓存“字段 id 的三分量”，这里存 int
struct IdTriplet
{
    int xi = -1;
    int eta = -1;
    int zeta = -1;

    bool all_valid() const { return (xi >= 0) && (eta >= 0) && (zeta >= 0); }

    void require_all(const char *what) const
    {
        if (!all_valid())
        {
            std::fprintf(stderr, "[IdTriplet] %s not fully bound (xi=%d eta=%d zeta=%d)\n",
                         what ? what : "(null)", xi, eta, zeta);
            std::abort();
        }
    }

    const int &at(int dir) const
    {
        switch (dir)
        {
        case 1:
            return xi;
        case 2:
            return eta;
        case 3:
            return zeta;
        default:
            std::fprintf(stderr, "[IdTriplet] invalid dir=%d (expect 0/1/2)\n", dir);
            std::abort();
        }
    }
};