// 0_basic/types.h
#pragma once
struct Int3
{
    int i, j, k;
};

struct Box3
{
    Int3 lo, hi; //[lo,hi)半开区间
};
