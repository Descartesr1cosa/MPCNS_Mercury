#include "0_basic/LayoutTraits.h"

#include "0_basic/BoxOps.h"
#include "0_basic/Error.h"

namespace
{
    int abs_int(int v)
    {
        return v < 0 ? -v : v;
    }

    int get_axis(const Int3 &v, int axis)
    {
        switch (axis)
        {
        case 0:
            return v.i;
        case 1:
            return v.j;
        case 2:
            return v.k;
        default:
            ERROR::Abort("LAYOUT: invalid axis");
        }
    }

    void set_axis(Int3 &v, int axis, int value)
    {
        switch (axis)
        {
        case 0:
            v.i = value;
            return;
        case 1:
            v.j = value;
            return;
        case 2:
            v.k = value;
            return;
        default:
            ERROR::Abort("LAYOUT: invalid axis");
        }
    }

    void get_tangent_axes(int normal_axis, int &t1, int &t2)
    {
        if (normal_axis == 0)
        {
            t1 = 1;
            t2 = 2;
        }
        else if (normal_axis == 1)
        {
            t1 = 0;
            t2 = 2;
        }
        else if (normal_axis == 2)
        {
            t1 = 0;
            t2 = 1;
        }
        else
        {
            ERROR::Abort("LAYOUT: invalid normal axis");
        }
    }

    void convert_node_range_to_dof(int lo_node, int hi_node, int delta, int &lo, int &hi)
    {
        lo = lo_node;
        hi = (delta == 0) ? hi_node : (hi_node - 1);
    }

    void set_axis_range(Box3 &b, int axis, int lo, int hi)
    {
        set_axis(b.lo, axis, lo);
        set_axis(b.hi, axis, hi);
    }

    void set_tangent_from_node_box(Box3 &b, const Box3 &node_box, const Int3 &delta, int axis)
    {
        int lo = 0;
        int hi = 0;
        convert_node_range_to_dof(get_axis(node_box.lo, axis),
                                  get_axis(node_box.hi, axis),
                                  get_axis(delta, axis),
                                  lo,
                                  hi);
        set_axis_range(b, axis, lo, hi);
    }

    void apply_inner_strip(Box3 &b, StaggerLocation loc, int dir_code, int nghost)
    {
        const int axis = LAYOUT::axis_from_dir(dir_code);
        const int delta = get_axis(LAYOUT::dof_delta(loc), axis);
        const int exclude = (delta == 0) ? 1 : 0;

        if (LAYOUT::sign_from_dir(dir_code) < 0)
        {
            const int lo = get_axis(b.lo, axis) + exclude;
            set_axis_range(b, axis, lo, lo + nghost);
        }
        else
        {
            const int hi = get_axis(b.hi, axis) - exclude;
            set_axis_range(b, axis, hi - nghost, hi);
        }
    }

    void apply_ghost_strip(Box3 &b, int dir_code, int nghost)
    {
        const int axis = LAYOUT::axis_from_dir(dir_code);

        if (LAYOUT::sign_from_dir(dir_code) < 0)
        {
            set_axis_range(b, axis, -nghost, 0);
        }
        else
        {
            const int lo = get_axis(b.hi, axis);
            set_axis_range(b, axis, lo, lo + nghost);
        }
    }

    void set_ghost_range_from_inner_hi(Box3 &b, const Int3 &inner_hi, int dir_code, int nghost)
    {
        const int axis = LAYOUT::axis_from_dir(dir_code);
        const int hi_axis = get_axis(inner_hi, axis);

        if (LAYOUT::sign_from_dir(dir_code) < 0)
            set_axis_range(b, axis, -nghost, 0);
        else
            set_axis_range(b, axis, hi_axis, hi_axis + nghost);
    }

    void assert_distinct_axes(int a, int b, const char *where)
    {
        if (LAYOUT::axis_from_dir(a) == LAYOUT::axis_from_dir(b))
            ERROR::Abort(where);
    }

    void assert_distinct_axes(int a, int b, int c, const char *where)
    {
        const int ax_a = LAYOUT::axis_from_dir(a);
        const int ax_b = LAYOUT::axis_from_dir(b);
        const int ax_c = LAYOUT::axis_from_dir(c);
        if (ax_a == ax_b || ax_a == ax_c || ax_b == ax_c)
            ERROR::Abort(where);
    }
}

namespace LAYOUT
{
    Int3 node_size_from_cells(const Int3 &ncells)
    {
        return {ncells.i + 1, ncells.j + 1, ncells.k + 1};
    }

    Int3 dof_delta(StaggerLocation loc)
    {
        switch (loc)
        {
        case StaggerLocation::Cell:
            return {1, 1, 1};
        case StaggerLocation::Node:
            return {0, 0, 0};
        case StaggerLocation::FaceXi:
            return {0, 1, 1};
        case StaggerLocation::FaceEt:
            return {1, 0, 1};
        case StaggerLocation::FaceZe:
            return {1, 1, 0};
        case StaggerLocation::EdgeXi:
            return {1, 0, 0};
        case StaggerLocation::EdgeEt:
            return {0, 1, 0};
        case StaggerLocation::EdgeZe:
            return {0, 0, 1};
        default:
            ERROR::Abort("LAYOUT::dof_delta: unsupported StaggerLocation");
        }
    }

    const char *location_name(StaggerLocation loc)
    {
        switch (loc)
        {
        case StaggerLocation::Cell:
            return "Cell";
        case StaggerLocation::Node:
            return "Node";
        case StaggerLocation::FaceXi:
            return "FaceXi";
        case StaggerLocation::FaceEt:
            return "FaceEt";
        case StaggerLocation::FaceZe:
            return "FaceZe";
        case StaggerLocation::EdgeXi:
            return "EdgeXi";
        case StaggerLocation::EdgeEt:
            return "EdgeEt";
        case StaggerLocation::EdgeZe:
            return "EdgeZe";
        default:
            return "Unknown";
        }
    }

    bool is_cell_location(StaggerLocation loc)
    {
        return loc == StaggerLocation::Cell;
    }

    bool is_node_location(StaggerLocation loc)
    {
        return loc == StaggerLocation::Node;
    }

    bool is_face_location(StaggerLocation loc)
    {
        return loc == StaggerLocation::FaceXi ||
               loc == StaggerLocation::FaceEt ||
               loc == StaggerLocation::FaceZe;
    }

    bool is_edge_location(StaggerLocation loc)
    {
        return loc == StaggerLocation::EdgeXi ||
               loc == StaggerLocation::EdgeEt ||
               loc == StaggerLocation::EdgeZe;
    }

    int face_axis(StaggerLocation loc)
    {
        switch (loc)
        {
        case StaggerLocation::FaceXi:
            return 0;
        case StaggerLocation::FaceEt:
            return 1;
        case StaggerLocation::FaceZe:
            return 2;
        default:
            ERROR::Abort("LAYOUT::face_axis: location is not a face");
        }
    }

    int edge_axis(StaggerLocation loc)
    {
        switch (loc)
        {
        case StaggerLocation::EdgeXi:
            return 0;
        case StaggerLocation::EdgeEt:
            return 1;
        case StaggerLocation::EdgeZe:
            return 2;
        default:
            ERROR::Abort("LAYOUT::edge_axis: location is not an edge");
        }
    }

    int codim(StaggerLocation loc)
    {
        if (is_cell_location(loc))
            return 0;
        if (is_face_location(loc))
            return 1;
        if (is_edge_location(loc))
            return 2;
        if (is_node_location(loc))
            return 3;

        ERROR::Abort("LAYOUT::codim: unsupported StaggerLocation");
    }

    Box3 owned_box_from_cells(const Int3 &ncells, StaggerLocation loc)
    {
        const Int3 nodes = node_size_from_cells(ncells);
        const Int3 delta = dof_delta(loc);
        return Box3{{0, 0, 0},
                    {nodes.i - delta.i,
                     nodes.j - delta.j,
                     nodes.k - delta.k}};
    }

    Box3 allocated_box_from_cells(const Int3 &ncells, StaggerLocation loc, int nghost)
    {
        Box3 b = owned_box_from_cells(ncells, loc);
        b.lo = {-nghost, -nghost, -nghost};
        b.hi = {b.hi.i + nghost, b.hi.j + nghost, b.hi.k + nghost};
        return b;
    }

    Box3 node_box_to_dof_box(StaggerLocation loc, const Box3 &node_box)
    {
        const Int3 delta = dof_delta(loc);
        return Box3{node_box.lo,
                    {node_box.hi.i - delta.i,
                     node_box.hi.j - delta.j,
                     node_box.hi.k - delta.k}};
    }

    bool is_valid_dir(int dir_code)
    {
        const int a = abs_int(dir_code);
        return a >= 1 && a <= 3;
    }

    int axis_from_dir(int dir_code)
    {
        if (!is_valid_dir(dir_code))
            ERROR::Abort("LAYOUT::axis_from_dir: invalid direction code");
        return abs_int(dir_code) - 1;
    }

    int sign_from_dir(int dir_code)
    {
        if (!is_valid_dir(dir_code))
            ERROR::Abort("LAYOUT::sign_from_dir: invalid direction code");
        return dir_code < 0 ? -1 : 1;
    }

    Box3 neighbor_inner_strip_from_node_box(StaggerLocation loc,
                                            const Box3 &node_box,
                                            int dir_code,
                                            int nghost)
    {
        Box3 b = node_box_to_dof_box(loc, node_box);
        apply_inner_strip(b, loc, dir_code, nghost);
        return b;
    }

    Box3 ghost_strip_from_node_box(StaggerLocation loc,
                                   const Box3 &node_box,
                                   int dir_code,
                                   int nghost)
    {
        Box3 b = node_box_to_dof_box(loc, node_box);
        apply_ghost_strip(b, dir_code, nghost);
        return b;
    }

    Box3 corner2_ghost_from_node_box(StaggerLocation loc,
                                     const Box3 &node_box,
                                     int dir1,
                                     int dir2,
                                     int nghost)
    {
        assert_distinct_axes(dir1, dir2, "LAYOUT::corner2_ghost_from_node_box: directions share an axis");
        Box3 b = node_box_to_dof_box(loc, node_box);
        apply_ghost_strip(b, dir1, nghost);
        apply_ghost_strip(b, dir2, nghost);
        return b;
    }

    Box3 corner2_innerghost_from_node_box(StaggerLocation loc,
                                          const Box3 &node_box,
                                          int inner_dir,
                                          int ghost_dir,
                                          int nghost)
    {
        assert_distinct_axes(inner_dir, ghost_dir, "LAYOUT::corner2_innerghost_from_node_box: directions share an axis");
        Box3 b = node_box_to_dof_box(loc, node_box);
        apply_inner_strip(b, loc, inner_dir, nghost);
        apply_ghost_strip(b, ghost_dir, nghost);
        return b;
    }

    Box3 corner3_ghost_from_node_box(StaggerLocation loc,
                                     const Box3 &node_box,
                                     int dir1,
                                     int dir2,
                                     int dir3,
                                     int nghost)
    {
        assert_distinct_axes(dir1, dir2, dir3, "LAYOUT::corner3_ghost_from_node_box: directions share an axis");
        Box3 b = node_box_to_dof_box(loc, node_box);
        apply_ghost_strip(b, dir1, nghost);
        apply_ghost_strip(b, dir2, nghost);
        apply_ghost_strip(b, dir3, nghost);
        return b;
    }

    Box3 corner3_innerghost_from_node_box(StaggerLocation loc,
                                          const Box3 &node_box,
                                          int inner_dir,
                                          int ghost_dir1,
                                          int ghost_dir2,
                                          int nghost)
    {
        assert_distinct_axes(inner_dir,
                             ghost_dir1,
                             ghost_dir2,
                             "LAYOUT::corner3_innerghost_from_node_box: directions share an axis");
        Box3 b = node_box_to_dof_box(loc, node_box);
        apply_inner_strip(b, loc, inner_dir, nghost);
        apply_ghost_strip(b, ghost_dir1, nghost);
        apply_ghost_strip(b, ghost_dir2, nghost);
        return b;
    }

    Box3 boundary_inner_slab_one_layer_from_cells(const Int3 &ncells,
                                                  StaggerLocation loc,
                                                  const Box3 &face_node_box,
                                                  int dir_code)
    {
        const int axis = axis_from_dir(dir_code);
        const int sign = sign_from_dir(dir_code);
        const Int3 inner_hi = owned_box_from_cells(ncells, loc).hi;
        const Int3 delta = dof_delta(loc);

        int t1 = 0;
        int t2 = 0;
        get_tangent_axes(axis, t1, t2);

        Box3 b{};
        set_tangent_from_node_box(b, face_node_box, delta, t1);
        set_tangent_from_node_box(b, face_node_box, delta, t2);

        if (sign < 0)
        {
            set_axis_range(b, axis, 0, 1);
        }
        else
        {
            const int hi_axis = get_axis(inner_hi, axis);
            set_axis_range(b, axis, hi_axis - 1, hi_axis);
        }

        return b;
    }

    Box3 ghost_slab_from_inner(const Box3 &inner_slab, int dir_code, int nghost)
    {
        const int axis = axis_from_dir(dir_code);
        Box3 b = inner_slab;

        if (sign_from_dir(dir_code) < 0)
        {
            const int hi = get_axis(inner_slab.lo, axis);
            set_axis_range(b, axis, hi - nghost, hi);
        }
        else
        {
            const int lo = get_axis(inner_slab.hi, axis);
            set_axis_range(b, axis, lo, lo + nghost);
        }

        return b;
    }

    Box3 coupling_face_ghost_slab_from_cells(const Int3 &ncells,
                                             StaggerLocation loc,
                                             const Box3 &face_node_box,
                                             int dir_code,
                                             int nghost)
    {
        const Box3 inner = boundary_inner_slab_one_layer_from_cells(ncells, loc, face_node_box, dir_code);
        return ghost_slab_from_inner(inner, dir_code, nghost);
    }

    Box3 coupling_edge_ghost_slab_from_cells(const Int3 &ncells,
                                             StaggerLocation loc,
                                             const Box3 &edge_node_box,
                                             int dir1,
                                             int dir2,
                                             int nghost)
    {
        assert_distinct_axes(dir1, dir2, "LAYOUT::coupling_edge_ghost_slab_from_cells: directions share an axis");

        const int axis1 = axis_from_dir(dir1);
        const int axis2 = axis_from_dir(dir2);
        const int edge_axis = 3 - axis1 - axis2;
        const Int3 inner_hi = owned_box_from_cells(ncells, loc).hi;
        const Int3 delta = dof_delta(loc);

        Box3 b{};
        set_ghost_range_from_inner_hi(b, inner_hi, dir1, nghost);
        set_ghost_range_from_inner_hi(b, inner_hi, dir2, nghost);
        set_tangent_from_node_box(b, edge_node_box, delta, edge_axis);
        return b;
    }

    Box3 coupling_vertex_ghost_slab_from_cells(const Int3 &ncells,
                                               StaggerLocation loc,
                                               int dir1,
                                               int dir2,
                                               int dir3,
                                               int nghost)
    {
        assert_distinct_axes(dir1, dir2, dir3, "LAYOUT::coupling_vertex_ghost_slab_from_cells: directions share an axis");
        const Int3 inner_hi = owned_box_from_cells(ncells, loc).hi;

        Box3 b{};
        set_ghost_range_from_inner_hi(b, inner_hi, dir1, nghost);
        set_ghost_range_from_inner_hi(b, inner_hi, dir2, nghost);
        set_ghost_range_from_inner_hi(b, inner_hi, dir3, nghost);
        return b;
    }

    bool box_empty(const Box3 &b)
    {
        return (b.hi.i <= b.lo.i) ||
               (b.hi.j <= b.lo.j) ||
               (b.hi.k <= b.lo.k);
    }

    Int3 box_size(const Box3 &b)
    {
        return {b.hi.i - b.lo.i,
                b.hi.j - b.lo.j,
                b.hi.k - b.lo.k};
    }

    void assert_valid_box(const Box3 &b, const char *where)
    {
        if (b.hi.i < b.lo.i ||
            b.hi.j < b.lo.j ||
            b.hi.k < b.lo.k)
        {
            ERROR::Abort(where ? where : "LAYOUT::assert_valid_box");
        }
    }

    void assert_nonempty_box(const Box3 &b, const char *where)
    {
        BOX::assert_nonempty(b, where ? where : "LAYOUT::assert_nonempty_box");
    }

    void assert_box_inside(const Box3 &inner, const Box3 &outer, const char *where)
    {
        if (inner.lo.i < outer.lo.i ||
            inner.lo.j < outer.lo.j ||
            inner.lo.k < outer.lo.k ||
            inner.hi.i > outer.hi.i ||
            inner.hi.j > outer.hi.j ||
            inner.hi.k > outer.hi.k)
        {
            ERROR::Abort(where ? where : "LAYOUT::assert_box_inside");
        }
    }
}
