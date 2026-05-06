#include "1_grid/1_MPCNS_Grid.h"
#include <algorithm>

void Block::calc_Corner_Preprocess()
{

    edge.clear(); // 清空当前棱边列表

    if (dimension == 2)
    {
        // 遍历 Parallel_Boundary，生成棱边
        for (auto &parallel : parallel_bc)
            create_edges_from_boundary(parallel.sub, parallel.sup, parallel.direction, &parallel, nullptr, nullptr, dimension);

        // 遍历 Inner_Boundary，生成棱边
        for (auto &inner : inner_bc)
            create_edges_from_boundary(inner.sub, inner.sup, inner.direction, nullptr, &inner, nullptr, dimension);

        // 遍历 Physical_Boundary，生成棱边
        for (auto &physical : physical_bc)
            create_edges_from_boundary(physical.sub, physical.sup, physical.direction, nullptr, nullptr, &physical, dimension);

        judge_edges(dimension);
    }
    else
    {
        // 遍历 Parallel_Boundary，生成棱边
        for (auto &parallel : parallel_bc)
            create_edges_from_boundary(parallel.sub, parallel.sup, parallel.direction, &parallel, nullptr, nullptr);

        // 遍历 Inner_Boundary，生成棱边
        for (auto &inner : inner_bc)
            create_edges_from_boundary(inner.sub, inner.sup, inner.direction, nullptr, &inner, nullptr);

        // 遍历 Physical_Boundary，生成棱边
        for (auto &physical : physical_bc)
            create_edges_from_boundary(physical.sub, physical.sup, physical.direction, nullptr, nullptr, &physical);

        judge_edges();
    }
}

void Block::create_edges_from_boundary(const int32_t sub[3], const int32_t sup[3], int32_t direction,
                                       Parallel_Boundary *para_bound, Inner_Boundary *inner_bound, Physical_Boundary *phy_bound)
{

    int a[3], b[3], c[3], d[3], temp = 0;
    for (int i = 0; i < 3; ++i)
    {
        if (sub[i] != sup[i])
        {
            temp++;
            if (temp == 1)
            {
                a[i] = sub[i];
                b[i] = sub[i];
                c[i] = sup[i];
                d[i] = sup[i];
            }
            else if (temp == 2)
            {
                a[i] = sub[i];
                b[i] = sup[i];
                c[i] = sub[i];
                d[i] = sup[i];
            }
        }
        else
        {
            a[i] = sub[i];
            b[i] = sub[i];
            c[i] = sub[i];
            d[i] = sub[i];
        }
    }

    Edge new_edge[4];
    for (int i = 0; i < 3; ++i)
    {
        // ac
        new_edge[0].sub[i] = a[i];
        new_edge[0].sup[i] = c[i];

        // ab
        new_edge[1].sub[i] = a[i];
        new_edge[1].sup[i] = b[i];

        // bd
        new_edge[2].sub[i] = b[i];
        new_edge[2].sup[i] = d[i];

        // cd
        new_edge[3].sub[i] = c[i];
        new_edge[3].sup[i] = d[i];
    }

    int face_dir;
    if (para_bound)
    {
        for (int j = 0; j < 4; ++j)
        {
            new_edge[j].para_bound = para_bound;
            new_edge[j].index = 1; // 并行边界
            new_edge[j].this_block_num = para_bound->this_block_num;
            face_dir = para_bound->direction;
        }
    }
    else if (inner_bound)
    {
        for (int j = 0; j < 4; ++j)
        {
            new_edge[j].inner_bound = inner_bound;
            new_edge[j].index = 0; // 内部边界
            new_edge[j].this_block_num = inner_bound->this_block_num;
            face_dir = inner_bound->direction;
        }
    }
    else if (phy_bound)
    {
        for (int j = 0; j < 4; ++j)
        {
            new_edge[j].phy_bound = phy_bound;
            new_edge[j].index = 2; // 物理边界
            new_edge[j].this_block_num = phy_bound->this_block_num;
            face_dir = phy_bound->direction;
        }
    }

    int tempdir;
    for (int j = 0; j < 4; ++j)
    {
        for (int i = 0; i < 3; ++i)
        {
            if (new_edge[j].sub[i] != new_edge[j].sup[i])
            {
                tempdir = i + 1;
                new_edge[j].cycle1[i] = 0;
            }
            else
            {
                if (abs(new_edge[j].sub[i]) == std::max(abs(sub[i]), abs(sup[i])))
                    new_edge[j].cycle1[i] = 1;
                else
                    new_edge[j].cycle1[i] = -1;
            }
        }
        new_edge[j].cycle1[abs(face_dir) - 1] = 0;
        new_edge[j].direction1 = 6 - abs(face_dir) - tempdir;

        new_edge[j].direction1 = (new_edge[j].cycle1[new_edge[j].direction1 - 1] > 0) ? new_edge[j].direction1 : -new_edge[j].direction1;
        add_edge(new_edge[j]);
    }
}

void Block::create_edges_from_boundary(const int32_t sub[3], const int32_t sup[3], int32_t direction,
                                       Parallel_Boundary *para_bound, Inner_Boundary *inner_bound, Physical_Boundary *phy_bound, int _index)
{
    int a[3], b[3];
    for (int i = 0; i < 3; ++i)
    {
        if (sub[i] != sup[i])
        {
            a[i] = sub[i];
            b[i] = sup[i];
        }
        else
        {
            a[i] = sub[i];
            b[i] = sub[i];
        }
    }

    Edge new_edge[2];
    for (int i = 0; i < 3; ++i)
    {
        // a
        new_edge[0].sub[i] = a[i];
        new_edge[0].sup[i] = a[i];

        // b
        new_edge[1].sub[i] = b[i];
        new_edge[1].sup[i] = b[i];
    }

    int face_dir;
    if (para_bound)
    {
        for (int j = 0; j < 2; ++j)
        {
            new_edge[j].para_bound = para_bound;
            new_edge[j].index = 1; // 并行边界
            new_edge[j].this_block_num = para_bound->this_block_num;
            face_dir = para_bound->direction;
        }
    }
    else if (inner_bound)
    {
        for (int j = 0; j < 2; ++j)
        {
            new_edge[j].inner_bound = inner_bound;
            new_edge[j].index = 0; // 内部边界
            new_edge[j].this_block_num = inner_bound->this_block_num;
            face_dir = inner_bound->direction;
        }
    }
    else if (phy_bound)
    {
        for (int j = 0; j < 2; ++j)
        {
            new_edge[j].phy_bound = phy_bound;
            new_edge[j].index = 2; // 物理边界
            new_edge[j].this_block_num = phy_bound->this_block_num;
            face_dir = phy_bound->direction;
        }
    }

    for (int j = 0; j < 2; ++j)
    {
        new_edge[j].cycle1[2] = 0;
        new_edge[j].cycle1[abs(face_dir) - 1] = 0;
        new_edge[j].direction1 = 3 - abs(face_dir);

        int direction1 = new_edge[j].direction1;
        if (new_edge[j].sub[direction1 - 1] == std::min(abs(a[direction1 - 1]), abs(b[direction1 - 1])))
            new_edge[j].direction1 = -direction1;

        new_edge[j].cycle1[direction1 - 1] = (new_edge[j].direction1 > 0) ? 1 : -1;
        add_edge(new_edge[j], _index);
    }
}

// 用于比较两条边是否重叠或包含
bool Block::is_overlapping(const Edge &e1, const Edge &e2, int &dir)
{
    int dir1 = -1, dir2 = -1;
    int sub1[3] = {e1.sub[0], e1.sub[1], e1.sub[2]};
    int sup1[3] = {e1.sup[0], e1.sup[1], e1.sup[2]};
    int sub2[3] = {e2.sub[0], e2.sub[1], e2.sub[2]};
    int sup2[3] = {e2.sup[0], e2.sup[1], e2.sup[2]};

    // 查找dir1
    for (size_t i = 0; i < dimension; i++)
    {
        if (e1.sub[i] != e1.sup[i])
        {
            dir1 = i;
            sub1[dir1] = 0; // 压缩 sub1 和 sup1 成一个点
            break;
        }
    }

    // 查找dir2
    for (size_t i = 0; i < dimension; i++)
    {
        if (e2.sub[i] != e2.sup[i])
        {
            dir2 = i;
            sub2[dir2] = 0; // 压缩 sub2 和 sup2 成一个点
            break;
        }
    }

    // 如果方向不同，则不重叠
    if (dir1 != dir2)
        return false;

    // 判断在dir1方向上的重叠
    int e1_min = std::min(abs(e1.sub[dir1]), abs(e1.sup[dir1]));
    int e1_max = std::max(abs(e1.sub[dir1]), abs(e1.sup[dir1]));
    int e2_min = std::min(abs(e2.sub[dir1]), abs(e2.sup[dir1]));
    int e2_max = std::max(abs(e2.sub[dir1]), abs(e2.sup[dir1]));

    bool overlap1 = (e1_min < e2_max) && (e1_max > e2_min);

    // 判断在非dir1方向上是否重合
    bool overlap2 = (abs(sub1[0]) == abs(sub2[0])) &&
                    (abs(sub1[1]) == abs(sub2[1])) &&
                    (abs(sub1[2]) == abs(sub2[2]));

    if (overlap1 && overlap2)
    {
        dir = dir1;
        return true;
    }

    return false;
}

// 用于比较两点是否重合， For 2D
bool Block::is_same(const Edge &e1, const Edge &e2)
{
    int sub1[3] = {e1.sub[0], e1.sub[1], e1.sub[2]};
    // int sup1[3] = {e1.sup[0], e1.sup[1], e1.sup[2]};
    int sub2[3] = {e2.sub[0], e2.sub[1], e2.sub[2]};
    // int sup2[3] = {e2.sup[0], e2.sup[1], e2.sup[2]};

    bool overlap2 = (abs(sub1[0]) == abs(sub2[0])) &&
                    (abs(sub1[1]) == abs(sub2[1])) &&
                    (abs(sub1[2]) == abs(sub2[2]));

    return overlap2;
}

void Block::add_edge(Edge &e1)
{
    bool is_find = false;
    int dir, index_find;

    for (int i = 0; i < edge.size(); i++)
    {
        Edge &temp = edge[i];
        if (is_overlapping(e1, temp, dir))
        {
            is_find = true;
            index_find = i;
            break;
        }
    }
    if (!is_find) // 没有找到edge，直接添加入edge即可
        edge.push_back(e1);
    else
    {
        // 已知棱的方向为dir;
        Edge &e2 = edge[index_find];
        // 获取线1和线2在 dir 方向的范围
        int range1_start = std::min(abs(e1.sub[dir]), abs(e1.sup[dir]));
        int range1_end = std::max(abs(e1.sub[dir]), abs(e1.sup[dir]));
        int range2_start = std::min(abs(e2.sub[dir]), abs(e2.sup[dir]));
        int range2_end = std::max(abs(e2.sub[dir]), abs(e2.sup[dir]));

        // 计算重叠的起始和结束
        int overlap_start = std::max(range1_start, range2_start);
        int overlap_end = std::min(range1_end, range2_end);

        auto copy_edge_with_split = [&](Edge &edge, int split_start, int split_end)
        {
            if (abs(edge.sub[dir]) < abs(edge.sup[dir]))
            {
                edge.sub[dir] = (edge.sup[dir] < 0) ? -split_start : split_start;
                edge.sup[dir] = (edge.sup[dir] < 0) ? -split_end : split_end;
            }
            else
            {
                edge.sup[dir] = (edge.sub[dir] < 0) ? -split_start : split_start;
                edge.sub[dir] = (edge.sub[dir] < 0) ? -split_end : split_end;
            }
        };

        // 非重合部分生成新线段
        if (overlap_start > range1_start)
        {
            Edge e1_copy(e1);
            copy_edge_with_split(e1_copy, range1_start, overlap_start);
            add_edge(e1_copy);
        }
        if (overlap_end < range1_end)
        {
            Edge e1_copy(e1);
            copy_edge_with_split(e1_copy, overlap_end, range1_end);
            add_edge(e1_copy);
        }
        if (overlap_start > range2_start)
        {
            // Edge e2_copy = e2;
            Edge e2_copy(edge[index_find]);
            copy_edge_with_split(e2_copy, range2_start, overlap_start);
            edge.push_back(e2_copy);
        }
        if (overlap_end < range2_end)
        {
            // Edge e2_copy = e2;
            Edge e2_copy(edge[index_find]);
            copy_edge_with_split(e2_copy, overlap_end, range2_end);
            edge.push_back(e2_copy);
        }

        // 重叠部分融合
        copy_edge_with_split(edge[index_find], overlap_start, overlap_end);
        copy_edge_with_split(e1, overlap_start, overlap_end);
        integration_equal_edge(edge[index_find], e1);
    }
}

void Block::add_edge(Edge &e1, int _index)
{
    // For 2D
    bool is_find = false;
    Edge *temp_edge = nullptr;

    for (Edge &temp : edge)
    {
        if (is_same(e1, temp))
        {
            is_find = true;
            temp_edge = &temp;
        }
    }
    if (!is_find) // 没有找到edge，直接添加入edge即可
        edge.push_back(e1);
    else
    {
        Edge &e2 = *temp_edge;
        // 融合
        integration_equal_edge(e2, e1);
    }
}

void Block::integration_equal_edge(Edge &e2, Edge &e1)
{

    auto e2_higher_priority = [](Edge &e2, Edge &e1)
    {
        for (int i = 0; i < 3; i++)
        {
            e2.cycle2[i] = e1.cycle1[i];
        }
        e2.direction2 = e1.direction1;
    };
    auto e1_higher_priority = [](Edge &e2, Edge &e1)
    {
        for (int i = 0; i < 3; i++)
        {
            e2.cycle2[i] = e2.cycle1[i];
            e2.cycle1[i] = e1.cycle1[i];
            e2.sub[i] = e1.sub[i];
            e2.sup[i] = e1.sup[i];
        }
        e2.direction2 = e2.direction1;
        e2.direction1 = e1.direction1;

        e2.index = e1.index;
        e2.para_bound = e1.para_bound;
        e2.inner_bound = e1.inner_bound;
        e2.phy_bound = e1.phy_bound;
    };

    int index1 = e1.index, index2 = e2.index;
    if (e1.index == 0)
        index1 += (e1.inner_bound->this_block_name != e1.inner_bound->target_block_name) ? 3 : 0;
    else if (e1.index == 1)
        index1 += (e1.para_bound->this_block_name != e1.para_bound->target_block_name) ? 3 : 0;

    if (e2.index == 0)
        index2 += (e2.inner_bound->this_block_name != e2.inner_bound->target_block_name) ? 3 : 0;
    else if (e2.index == 1)
        index2 += (e2.para_bound->this_block_name != e2.para_bound->target_block_name) ? 3 : 0;

    if (index2 <= index1)
    {
        if (index2 < index1)
            e2_higher_priority(e2, e1);
        else if (index2 != 2)
            e2_higher_priority(e2, e1);
        else
        {
            List<int> Priority = par->GetInt_List("Boundary_Priority");

            std::string bdname1, bdname2;
            bdname1 = e1.phy_bound->boundary_name;
            bdname2 = e2.phy_bound->boundary_name;

            if (Priority.data[bdname2] <= Priority.data[bdname1])
                e2_higher_priority(e2, e1);
            else
                e1_higher_priority(e2, e1);
        }
    }
    else
        e1_higher_priority(e2, e1);
}

void Block::judge_edges()
{
    auto judge_singular = [&](Edge &e) -> bool
    {
        int cycle11[3] = {0, 0, 0}, cycle22[3] = {0, 0, 0};
        cycle11[abs(e.direction1) - 1] = e.cycle1[abs(e.direction1) - 1];

        if (e.para_bound)
            cycle22[abs(e.para_bound->direction) - 1] = (e.para_bound->direction > 0) ? 1 : -1;
        else if (e.inner_bound)
            cycle22[abs(e.inner_bound->direction) - 1] = (e.inner_bound->direction > 0) ? 1 : -1;
        else if (e.phy_bound)
            cycle22[abs(e.phy_bound->direction) - 1] = (e.phy_bound->direction > 0) ? 1 : -1;

        double coor11[3], coor22[3];
        coor11[0] = x(abs(e.sub[0]) + cycle11[0], abs(e.sub[1]) + cycle11[1], abs(e.sub[2]) + cycle11[2]);
        coor11[1] = y(abs(e.sub[0]) + cycle11[0], abs(e.sub[1]) + cycle11[1], abs(e.sub[2]) + cycle11[2]);
        coor11[2] = z(abs(e.sub[0]) + cycle11[0], abs(e.sub[1]) + cycle11[1], abs(e.sub[2]) + cycle11[2]);
        coor22[0] = x(abs(e.sub[0]) + cycle22[0], abs(e.sub[1]) + cycle22[1], abs(e.sub[2]) + cycle22[2]);
        coor22[1] = y(abs(e.sub[0]) + cycle22[0], abs(e.sub[1]) + cycle22[1], abs(e.sub[2]) + cycle22[2]);
        coor22[2] = z(abs(e.sub[0]) + cycle22[0], abs(e.sub[1]) + cycle22[1], abs(e.sub[2]) + cycle22[2]);

        if (sqrt((coor22[0] - coor11[0]) * (coor22[0] - coor11[0]) + (coor22[1] - coor11[1]) * (coor22[1] - coor11[1]) + (coor22[2] - coor11[2]) * (coor22[2] - coor11[2])) < 1e-10)
        {
            e.if_large_singular = false;
            return true;
        }

        for (int i = 0; i < 3; i++)
        {
            cycle11[i] = 0;
            cycle22[i] = 0;
        }

        cycle11[abs(e.direction1) - 1] = e.cycle1[abs(e.direction1) - 1];
        cycle22[abs(e.direction1) - 1] = -e.cycle1[abs(e.direction1) - 1];
        if (e.para_bound)
        {
            cycle11[abs(e.para_bound->direction) - 1] = (e.para_bound->direction > 0) ? 1 : -1;
            cycle22[abs(e.para_bound->direction) - 1] = (e.para_bound->direction > 0) ? 1 : -1;
        }
        else if (e.inner_bound)
        {
            cycle11[abs(e.inner_bound->direction) - 1] = (e.inner_bound->direction > 0) ? 1 : -1;
            cycle22[abs(e.inner_bound->direction) - 1] = (e.inner_bound->direction > 0) ? 1 : -1;
        }
        else if (e.phy_bound)
        {
            cycle11[abs(e.phy_bound->direction) - 1] = (e.phy_bound->direction > 0) ? 1 : -1;
            cycle22[abs(e.phy_bound->direction) - 1] = (e.phy_bound->direction > 0) ? 1 : -1;
        }

        coor11[0] = x(abs(e.sub[0]) + cycle11[0], abs(e.sub[1]) + cycle11[1], abs(e.sub[2]) + cycle11[2]);
        coor11[1] = y(abs(e.sub[0]) + cycle11[0], abs(e.sub[1]) + cycle11[1], abs(e.sub[2]) + cycle11[2]);
        coor11[2] = z(abs(e.sub[0]) + cycle11[0], abs(e.sub[1]) + cycle11[1], abs(e.sub[2]) + cycle11[2]);
        coor22[0] = x(abs(e.sub[0]) + cycle22[0], abs(e.sub[1]) + cycle22[1], abs(e.sub[2]) + cycle22[2]);
        coor22[1] = y(abs(e.sub[0]) + cycle22[0], abs(e.sub[1]) + cycle22[1], abs(e.sub[2]) + cycle22[2]);
        coor22[2] = z(abs(e.sub[0]) + cycle22[0], abs(e.sub[1]) + cycle22[1], abs(e.sub[2]) + cycle22[2]);

        // if (par->GetInt("myid") == 3)
        // {
        //     std::cout << cycle11[0] << "   " << cycle11[1] << "   " << cycle11[2] << std::endl;
        //     std::cout << cycle22[0] << "   " << cycle22[1] << "   " << cycle22[2] << "   " << std::endl;
        //     std::cout << "judge_singular:  " << sqrt((coor22[0] - coor11[0]) * (coor22[0] - coor11[0]) + (coor22[1] - coor11[1]) * (coor22[1] - coor11[1]) + (coor22[2] - coor11[2]) * (coor22[2] - coor11[2])) << std::endl;
        //     if (sqrt((coor22[0] - coor11[0]) * (coor22[0] - coor11[0]) + (coor22[1] - coor11[1]) * (coor22[1] - coor11[1]) + (coor22[2] - coor11[2]) * (coor22[2] - coor11[2])) < 1e-10)
        //         std::cout << e.sub[0] << "  " << e.sup[0] << "  " << e.sub[1] << "  " << e.sup[1] << "  " << e.sub[2] << "  " << e.sup[2] << std::endl;
        // }

        if (sqrt((coor22[0] - coor11[0]) * (coor22[0] - coor11[0]) + (coor22[1] - coor11[1]) * (coor22[1] - coor11[1]) + (coor22[2] - coor11[2]) * (coor22[2] - coor11[2])) < 1e-10)
        {
            e.if_large_singular = true;
            return true;
        }

        else
            return false;
    };

    std::vector<Edge> temp_edge;
    temp_edge.resize(edge.size());
    for (int i = 0; i < edge.size(); i++)
        temp_edge[i] = edge[i];

    edge.clear();
    int edge_max[3] = {mx, my, mz};
    bool judge[3];
    for (Edge &e : temp_edge)
    {
        e.is_singular = judge_singular(e);
        if (e.is_singular)
            edge.push_back(e);
        else
        {
            for (int i = 0; i < 3; i++)
            {
                judge[i] = false;
                if (e.sub[i] == e.sup[i])
                    judge[i] = abs(e.sub[i]) == 0 || abs(e.sub[i]) == edge_max[i];
                else
                    judge[i] = true;
            }
            if (judge[0] && judge[1] && judge[2])
                edge.push_back(e);
        }
    }
}

void Block::judge_edges(int _index)
{
    // For 2D 但其实和三维一模一样
    auto judge_singular = [&](Edge e) -> bool
    {
        int cycle11[3] = {0}, cycle22[3] = {0};
        cycle11[abs(e.direction1) - 1] = e.cycle1[abs(e.direction1) - 1];

        if (e.para_bound)
            cycle22[abs(e.para_bound->direction) - 1] = (e.para_bound->direction > 0) ? 1 : -1;
        else if (e.inner_bound)
            cycle22[abs(e.inner_bound->direction) - 1] = (e.inner_bound->direction > 0) ? 1 : -1;
        else if (e.phy_bound)
            cycle22[abs(e.phy_bound->direction) - 1] = (e.phy_bound->direction > 0) ? 1 : -1;

        double coor11[3], coor22[3];
        coor11[0] = x(abs(e.sub[0]) + cycle11[0], abs(e.sub[1]) + cycle11[1], abs(e.sub[2]) + cycle11[2]);
        coor11[1] = y(abs(e.sub[0]) + cycle11[0], abs(e.sub[1]) + cycle11[1], abs(e.sub[2]) + cycle11[2]);
        coor11[2] = z(abs(e.sub[0]) + cycle11[0], abs(e.sub[1]) + cycle11[1], abs(e.sub[2]) + cycle11[2]);
        coor22[0] = x(abs(e.sub[0]) + cycle22[0], abs(e.sub[1]) + cycle22[1], abs(e.sub[2]) + cycle22[2]);
        coor22[1] = y(abs(e.sub[0]) + cycle22[0], abs(e.sub[1]) + cycle22[1], abs(e.sub[2]) + cycle22[2]);
        coor22[2] = z(abs(e.sub[0]) + cycle22[0], abs(e.sub[1]) + cycle22[1], abs(e.sub[2]) + cycle22[2]);

        if (sqrt((coor22[0] - coor11[0]) * (coor22[0] - coor11[0]) + (coor22[1] - coor11[1]) * (coor22[1] - coor11[1]) + (coor22[2] - coor11[2]) * (coor22[2] - coor11[2])) < 1e-10)
            return true;

        cycle11[3] = {0}, cycle22[3] = {0};

        cycle11[abs(e.direction1) - 1] = e.cycle1[abs(e.direction1) - 1];
        cycle22[abs(e.direction1) - 1] = -e.cycle1[abs(e.direction1) - 1];
        if (e.para_bound)
        {
            cycle11[abs(e.para_bound->direction) - 1] = (e.para_bound->direction > 0) ? 1 : -1;
            cycle22[abs(e.para_bound->direction) - 1] = (e.para_bound->direction > 0) ? 1 : -1;
        }
        else if (e.inner_bound)
        {
            cycle11[abs(e.inner_bound->direction) - 1] = (e.inner_bound->direction > 0) ? 1 : -1;
            cycle22[abs(e.inner_bound->direction) - 1] = (e.inner_bound->direction > 0) ? 1 : -1;
        }
        else if (e.phy_bound)
        {
            cycle11[abs(e.phy_bound->direction) - 1] = (e.phy_bound->direction > 0) ? 1 : -1;
            cycle22[abs(e.phy_bound->direction) - 1] = (e.phy_bound->direction > 0) ? 1 : -1;
        }

        coor11[0] = x(abs(e.sub[0]) + cycle11[0], abs(e.sub[1]) + cycle11[1], abs(e.sub[2]) + cycle11[2]);
        coor11[1] = y(abs(e.sub[0]) + cycle11[0], abs(e.sub[1]) + cycle11[1], abs(e.sub[2]) + cycle11[2]);
        coor11[2] = z(abs(e.sub[0]) + cycle11[0], abs(e.sub[1]) + cycle11[1], abs(e.sub[2]) + cycle11[2]);
        coor22[0] = x(abs(e.sub[0]) + cycle22[0], abs(e.sub[1]) + cycle22[1], abs(e.sub[2]) + cycle22[2]);
        coor22[1] = y(abs(e.sub[0]) + cycle22[0], abs(e.sub[1]) + cycle22[1], abs(e.sub[2]) + cycle22[2]);
        coor22[2] = z(abs(e.sub[0]) + cycle22[0], abs(e.sub[1]) + cycle22[1], abs(e.sub[2]) + cycle22[2]);

        if (sqrt((coor22[0] - coor11[0]) * (coor22[0] - coor11[0]) + (coor22[1] - coor11[1]) * (coor22[1] - coor11[1]) + (coor22[2] - coor11[2]) * (coor22[2] - coor11[2])) < 1e-10)
            return true;
        else
            return false;
    };

    std::vector<Edge> temp_edge = edge;
    edge.clear();
    int edge_max[3] = {mx, my, 0};
    bool judge[2];
    for (Edge &e : temp_edge)
    {
        e.is_singular = judge_singular(e);
        if (e.is_singular)
            edge.push_back(e);
        else
        {
            for (int i = 0; i < 2; i++)
            {
                judge[i] = (abs(e.sub[i]) == 0 || abs(e.sub[i]) == edge_max[i]);
            }
            if (judge[0] && judge[1])
                edge.push_back(e);
        }
    }
}
