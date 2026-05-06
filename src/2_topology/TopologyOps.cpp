#include "2_topology/TopologyOps.h"

#include "0_basic/BoxOps.h"
#include "0_basic/Direction.h"
#include "0_basic/Error.h"

#include <algorithm>
#include <cstdlib>
#include <sstream>

namespace TOPO
{
    int get_comp(const Int3 &x, int axis0)
    {
        if (axis0 == 0)
            return x.i;
        if (axis0 == 1)
            return x.j;
        if (axis0 == 2)
            return x.k;

        ERROR::Abort("TOPO::get_comp: invalid axis");
        return 0;
    }

    void set_comp(Int3 &x, int axis0, int value)
    {
        if (axis0 == 0)
            x.i = value;
        else if (axis0 == 1)
            x.j = value;
        else if (axis0 == 2)
            x.k = value;
        else
            ERROR::Abort("TOPO::set_comp: invalid axis");
    }

    Box3 make_node_box_from_subsup(const int sub[3], const int sup[3])
    {
        Box3 box{};

        for (int d = 0; d < 3; ++d)
        {
            const int s = std::abs(sub[d]);
            const int t = std::abs(sup[d]);

            const int a = std::min(s, t);
            const int b = std::max(s, t);

            set_comp(box.lo, d, a);

            // Grid_Boundary 里的 sub/sup 是闭区间端点；
            // Topology 统一使用 half-open node-space box。
            set_comp(box.hi, d, b + 1);
        }

        return box;
    }

    void fill_node_box_from_subsup(const int sub[3], const int sup[3], Box3 &box)
    {
        box = make_node_box_from_subsup(sub, sup);
    }

    IndexTransform identity_transform()
    {
        IndexTransform tr{};

        tr.perm[0] = 0;
        tr.perm[1] = 1;
        tr.perm[2] = 2;

        tr.sign[0] = 1;
        tr.sign[1] = 1;
        tr.sign[2] = 1;

        tr.offset = {0, 0, 0};

        return tr;
    }

    Int3 map_node_point(const Int3 &p_local, const IndexTransform &tr)
    {
        Int3 p_nb{0, 0, 0};

        for (int a = 0; a < 3; ++a)
        {
            const int b = tr.perm[a];
            const int v = tr.sign[a] * get_comp(p_local, a) + get_comp(tr.offset, a);
            set_comp(p_nb, b, v);
        }

        return p_nb;
    }

    Box3 transform_node_box(const Box3 &box, const IndexTransform &tr)
    {
        if (BOX::empty(box))
            return box;

        const Int3 lo_closed = box.lo;
        const Int3 hi_closed{
            box.hi.i - 1,
            box.hi.j - 1,
            box.hi.k - 1};

        const Int3 p0 = map_node_point(lo_closed, tr);
        const Int3 p1 = map_node_point(hi_closed, tr);

        Box3 out{};

        out.lo.i = std::min(p0.i, p1.i);
        out.lo.j = std::min(p0.j, p1.j);
        out.lo.k = std::min(p0.k, p1.k);

        out.hi.i = std::max(p0.i, p1.i) + 1;
        out.hi.j = std::max(p0.j, p1.j) + 1;
        out.hi.k = std::max(p0.k, p1.k) + 1;

        return out;
    }

    IndexTransform make_index_transform_from_boundary_arrays(
        const int this_sub[3],
        const int this_sup[3],
        const int this_transform[3],
        const int nb_sub[3],
        const int nb_sup[3],
        const int nb_transform[3],
        int this_direction,
        int nb_direction,
        int dimension,
        const char *context)
    {
        IndexTransform tr = identity_transform();

        int inverse_nb_transform[3] = {0, 0, 0};

        for (int d = 0; d < 3; ++d)
        {
            const int td = nb_transform[d];

            if (td < 0 || td > 2)
            {
                std::ostringstream oss;
                oss << context << ": invalid nb_transform[" << d << "]=" << td;
                ERROR::Abort(oss.str());
            }

            inverse_nb_transform[td] = d;
        }

        for (int d = 0; d < 3; ++d)
        {
            const int td = this_transform[d];

            if (td < 0 || td > 2)
            {
                std::ostringstream oss;
                oss << context << ": invalid this_transform[" << d << "]=" << td;
                ERROR::Abort(oss.str());
            }

            tr.perm[d] = inverse_nb_transform[td];
        }

        // sign
        for (int d = 0; d < 3; ++d)
        {
            const int my_sub = std::abs(this_sub[d]);
            const int my_sup = std::abs(this_sup[d]);

            const int nb_axis = tr.perm[d];

            const int tar_sub = std::abs(nb_sub[nb_axis]);
            const int tar_sup = std::abs(nb_sup[nb_axis]);

            if (dimension == 2 && d == 2)
            {
                tr.sign[d] = 0;
                continue;
            }

            // 法向维：本侧是一个面，sub == sup
            if (my_sub == my_sup && d < dimension)
            {
                if (tar_sub != tar_sup)
                {
                    std::ostringstream oss;
                    oss << context << ": normal direction mismatch at d=" << d;
                    ERROR::Abort(oss.str());
                }

                tr.sign[d] = (this_direction * nb_direction > 0) ? -1 : 1;
                continue;
            }

            if (std::abs(my_sup - my_sub) != std::abs(tar_sup - tar_sub))
            {
                std::ostringstream oss;
                oss << context << ": interval length mismatch at d=" << d
                    << ", this=[" << my_sub << "," << my_sup << "]"
                    << ", nb=[" << tar_sub << "," << tar_sup << "]";
                ERROR::Abort(oss.str());
            }

            tr.sign[d] = ((my_sup - my_sub) * (tar_sup - tar_sub) > 0) ? 1 : -1;
        }

        // offset
        int offset[3] = {0, 0, 0};

        for (int d = 0; d < 3; ++d)
        {
            if (dimension == 2 && d == 2)
            {
                offset[d] = 0;
                continue;
            }

            const int my_sub = std::abs(this_sub[d]);
            const int nb_axis = tr.perm[d];
            const int tar_sub = std::abs(nb_sub[nb_axis]);

            offset[d] = -tr.sign[d] * my_sub + tar_sub;
        }

        tr.offset = {offset[0], offset[1], offset[2]};

        return tr;
    }

    int patch_kind_priority(PatchKind k)
    {
        if (k == PatchKind::Inner)
            return 3;
        if (k == PatchKind::Parallel)
            return 2;
        if (k == PatchKind::Physical)
            return 1;
        return 0;
    }

    int patch_priority(PatchKind k, bool is_coupling)
    {
        const int non_coupling = is_coupling ? 0 : 1;
        return non_coupling * 10 + patch_kind_priority(k);
    }

    int edge_axis1_from_face_dirs(int dir1, int dir2)
    {
        if (!DIR::distinct_axes(dir1, dir2))
            ERROR::Abort("TOPO::edge_axis1_from_face_dirs: directions must have distinct axes");

        bool used[4] = {false, false, false, false};

        used[DIR::axis1(dir1)] = true;
        used[DIR::axis1(dir2)] = true;

        if (!used[1])
            return 1;
        if (!used[2])
            return 2;
        if (!used[3])
            return 3;

        ERROR::Abort("TOPO::edge_axis1_from_face_dirs: cannot infer edge axis");
        return 0;
    }

    void order_edge_dirs_by_owner(int owner_dir, int &dir1, int &dir2)
    {
        if (!DIR::is_valid(owner_dir))
            ERROR::Abort("TOPO::order_edge_dirs_by_owner: invalid owner_dir");

        if (dir1 == owner_dir)
            return;

        if (dir2 == owner_dir)
        {
            std::swap(dir1, dir2);
            return;
        }

        ERROR::Abort("TOPO::order_edge_dirs_by_owner: owner_dir is not in edge dirs");
    }

    void order_vertex_dirs_by_owner(int owner_dir,
                                    int &dir1,
                                    int &dir2,
                                    int &dir3)
    {
        if (!DIR::is_valid(owner_dir))
            ERROR::Abort("TOPO::order_vertex_dirs_by_owner: invalid owner_dir");

        if (dir1 == owner_dir)
            return;

        if (dir2 == owner_dir)
        {
            std::swap(dir1, dir2);
            return;
        }

        if (dir3 == owner_dir)
        {
            std::swap(dir1, dir3);
            return;
        }

        ERROR::Abort("TOPO::order_vertex_dirs_by_owner: owner_dir is not in vertex dirs");
    }
}
