#include "2_topology/TopologyBuilder.h"
#include "2_topology/TopologyBuildDetail.h"
#include "2_topology/TopologyOps.h"

#include "0_basic/MPI_WRAPPER.h"
#include "0_basic/Error.h"
#include "1_grid/1_MPCNS_Grid.h"
#include "1_grid/Grid_Boundary.h"

#include <unordered_map>
#include <algorithm>
#include <limits>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>

namespace TOPO
{
    namespace
    {
        // ============================================================
        // Small utilities
        // ============================================================

        std::string parallel_rank_suffix(int r)
        {
            if (r < 10)
                return "   " + std::to_string(r);
            if (r < 100)
                return "  " + std::to_string(r);
            if (r < 1000)
                return " " + std::to_string(r);

            return std::to_string(r);
        }

        std::string parallel_boundary_filename(int rank)
        {
            return "./CASE/geometry/boundary_condition/parallel" +
                   parallel_rank_suffix(rank) + ".txt";
        }

        // ============================================================
        // Patch factories
        // ============================================================

        PhysicalPatch make_physical_patch(const Physical_Boundary &phy,
                                          int my_rank)
        {
            PhysicalPatch p;

            p.this_rank = my_rank;
            p.this_block = phy.this_block_num;
            p.this_block_name = phy.this_block_name;

            p.bc_id = phy.boundary_num;
            p.bc_name = phy.boundary_name;

            p.direction = phy.direction;
            p.raw = &phy;

            p.this_box_node = make_node_box_from_subsup(phy.sub, phy.sup);

            return p;
        }

        InterfacePatch make_inner_interface_patch(const Inner_Boundary &inner,
                                                  int my_rank,
                                                  int dimension)
        {
            InterfacePatch interface;

            interface.kind = PatchKind::Inner;

            interface.this_rank = my_rank;
            interface.nb_rank = my_rank;

            interface.this_block = inner.this_block_num;
            interface.nb_block = inner.tar_block_num;

            interface.this_block_name = inner.this_block_name;
            interface.nb_block_name = inner.target_block_name;

            interface.is_coupling =
                (interface.this_block_name != interface.nb_block_name);

            interface.this_box_node =
                make_node_box_from_subsup(inner.sub, inner.sup);

            interface.nb_box_node =
                make_node_box_from_subsup(inner.tar_sub, inner.tar_sup);

            interface.direction = inner.direction;
            interface.nb_direction = inner.tar_direction;

            interface.trans = make_index_transform_from_boundary_arrays(
                inner.sub,
                inner.sup,
                inner.Transform,
                inner.tar_sub,
                inner.tar_sup,
                inner.tar_Transform,
                inner.direction,
                inner.tar_direction,
                dimension,
                "build_topology: inner interface");

            interface.send_flag = 0;
            interface.recv_flag = 0;

            return interface;
        }

        InterfacePatch make_parallel_interface_patch(const Parallel_Boundary &para,
                                                     const Parallel_Boundary &tar_para,
                                                     int my_rank,
                                                     int dimension)
        {
            InterfacePatch interface;

            interface.kind = PatchKind::Parallel;

            interface.this_rank = my_rank;
            interface.nb_rank = para.tar_myid;

            interface.this_block = para.this_block_num;
            interface.nb_block = tar_para.this_block_num;

            interface.this_block_name = para.this_block_name;
            interface.nb_block_name = para.target_block_name;

            interface.is_coupling =
                (interface.this_block_name != interface.nb_block_name);

            interface.this_box_node =
                make_node_box_from_subsup(para.sub, para.sup);

            interface.nb_box_node =
                make_node_box_from_subsup(tar_para.sub, tar_para.sup);

            interface.direction = para.direction;
            interface.nb_direction = tar_para.direction;

            interface.trans = make_index_transform_from_boundary_arrays(
                para.sub,
                para.sup,
                para.Transform,
                tar_para.sub,
                tar_para.sup,
                tar_para.Transform,
                para.direction,
                tar_para.direction,
                dimension,
                "build_topology: parallel interface");

            interface.send_flag = para.send_flag;
            interface.recv_flag = para.rece_flag;

            return interface;
        }

        // ============================================================
        // Physical patch stage
        // ============================================================

        void append_physical_patches(Grid &grid,
                                     Topology &topo,
                                     int my_rank)
        {
            for (int ib = 0; ib < grid.nblock; ++ib)
            {
                const Block &blk = grid.grids(ib);

                for (const auto &phy : blk.physical_bc)
                {
                    topo.physical_patches.push_back(
                        make_physical_patch(phy, my_rank));
                }
            }
        }

        void sort_physical_patches_by_priority(Grid &grid,
                                               Topology &topo)
        {
            auto Priority = grid.par->GetInt_List("Boundary_Priority");

            std::unordered_map<std::string, int32_t> pri;
            pri.reserve(Priority.data.size());

            for (const auto &kv : Priority.data)
                pri.emplace(kv.first, kv.second);

            auto get_pri = [&](const std::string &name) -> int32_t
            {
                auto it = pri.find(name);
                if (it != pri.end())
                    return it->second;

                return std::numeric_limits<int32_t>::min();
            };

            // priority 越大越靠后，也就是后应用/覆盖。
            std::stable_sort(
                topo.physical_patches.begin(),
                topo.physical_patches.end(),
                [&](const PhysicalPatch &a, const PhysicalPatch &b)
                {
                    const int32_t pa = get_pri(a.bc_name);
                    const int32_t pb = get_pri(b.bc_name);

                    if (pa != pb)
                        return pa < pb;

                    if (a.this_block != b.this_block)
                        return a.this_block < b.this_block;

                    if (a.direction != b.direction)
                        return a.direction < b.direction;

                    return a.bc_id < b.bc_id;
                });
        }

        // ============================================================
        // Inner interface stage
        // ============================================================

        void append_inner_interface_patches(Grid &grid,
                                            Topology &topo,
                                            int my_rank,
                                            int dimension)
        {
            for (int ib = 0; ib < grid.nblock; ++ib)
            {
                const Block &blk = grid.grids(ib);

                for (const auto &inner : blk.inner_bc)
                {
                    topo.inner_patches.push_back(
                        make_inner_interface_patch(inner, my_rank, dimension));
                }
            }
        }

        // ============================================================
        // Parallel boundary reading stage
        // ============================================================

        std::vector<Parallel_Boundary>
        read_parallel_boundaries_for_rank(int rank,
                                          int dimension,
                                          int my_rank_for_error)
        {
            const std::string filename = parallel_boundary_filename(rank);

            std::ifstream grdfile(filename, std::ios_base::in);
            if (!grdfile.is_open())
            {
                std::cout << "ERROR: cannot open " << filename
                          << " when Read_All_Para(), myid = "
                          << my_rank_for_error << std::endl;
                exit(-1);
            }

            int blk_num_file = 0;
            grdfile >> blk_num_file;

            int num_parallel_face = 0;
            grdfile >> num_parallel_face;

            std::vector<int> nface(blk_num_file);
            for (int izone = 0; izone < blk_num_file; ++izone)
                grdfile >> nface[izone];

            std::string dummy;
            std::getline(grdfile, dummy);
            std::getline(grdfile, dummy);

            std::vector<Parallel_Boundary> result;
            result.reserve(num_parallel_face);

            for (int izone = 0; izone < blk_num_file; ++izone)
            {
                const int iface_num = nface[izone];

                for (int iiface = 0; iiface < iface_num; ++iiface)
                {
                    Parallel_Boundary pbc;

                    int pointst[3], pointed[3];
                    int srid, sflag, rflag;
                    std::string nameed;

                    grdfile >> pointst[0] >> pointed[0] >> pointst[1] >> pointed[1] >> pointst[2] >> pointed[2] >> srid >> sflag >> rflag;

                    grdfile >> nameed;

                    pbc.this_myid = rank;
                    pbc.tar_myid = srid;
                    pbc.this_block_num = izone;

                    for (int d = 0; d < 3; ++d)
                    {
                        pbc.sub[d] = pointst[d];
                        pbc.sup[d] = pointed[d];
                    }

                    pbc.send_flag = sflag;
                    pbc.rece_flag = rflag;

                    pbc.this_block_name = "";
                    pbc.target_block_name = nameed;

                    int32_t dim32 = static_cast<int32_t>(dimension);
                    pbc.Pre_process(dim32);

                    result.push_back(pbc);
                }
            }

            grdfile.close();
            return result;
        }

        std::vector<std::vector<Parallel_Boundary>>
        read_all_parallel_boundaries(int my_rank,
                                     int dimension)
        {
            int rank_num = 0;
            PARALLEL::mpi_size(&rank_num);

            std::vector<std::vector<Parallel_Boundary>> parallel_bc_all;
            parallel_bc_all.resize(rank_num);

            for (int r = 0; r < rank_num; ++r)
            {
                parallel_bc_all[r] =
                    read_parallel_boundaries_for_rank(r, dimension, my_rank);
            }

            if (my_rank == 0)
            {
                std::cout << "\t--> Building TOPO Para Interface: "
                          << "Read all parallel_*.txt successfully!\n";
            }

            return parallel_bc_all;
        }

        Parallel_Boundary *
        find_matching_parallel_boundary(std::vector<Parallel_Boundary> &candidates,
                                        const Parallel_Boundary &para)
        {
            // 保持原有行为：只用 send/recv flag 互相匹配。
            // 后续可以升级为 flag + block name + direction + box 的强匹配。
            for (auto &tar_para_temp : candidates)
            {
                if (tar_para_temp.rece_flag == para.send_flag &&
                    tar_para_temp.send_flag == para.rece_flag)
                {
                    return &tar_para_temp;
                }
            }

            return nullptr;
        }

        // ============================================================
        // Parallel interface stage
        // ============================================================

        void append_parallel_interface_patches(
            Grid &grid,
            Topology &topo,
            std::vector<std::vector<Parallel_Boundary>> &parallel_bc_all,
            int my_rank,
            int dimension)
        {
            for (int ib = 0; ib < grid.nblock; ++ib)
            {
                const Block &blk = grid.grids(ib);

                for (const auto &para : blk.parallel_bc)
                {
                    if (para.tar_myid < 0 ||
                        para.tar_myid >= static_cast<int>(parallel_bc_all.size()))
                    {
                        ERROR::Abort("build_topology: para.tar_myid out of range");
                    }

                    Parallel_Boundary *tar_para_pointer =
                        find_matching_parallel_boundary(
                            parallel_bc_all[para.tar_myid],
                            para);

                    if (!tar_para_pointer)
                    {
                        std::cout << "FATAL Error, Can not find corresponding PARA info!\n";
                        exit(-1);
                    }

                    Parallel_Boundary &tar_para = *tar_para_pointer;

                    topo.parallel_patches.push_back(
                        make_parallel_interface_patch(
                            para,
                            tar_para,
                            my_rank,
                            dimension));
                }
            }
        }

        // ============================================================
        // Derived patch stage
        // ============================================================

        void build_derived_topology_patches(Grid &grid,
                                            Topology &topo,
                                            int dimension)
        {
            detail::build_edge_patches(grid, topo, dimension);
            detail::build_vertex_patches(grid, topo, dimension);

            // 追加 Coupled-* 到 physical_patches。
            // 保持原有行为：不重新排序 physical_patches。
            detail::append_coupling_faces_as_physical_patches(
                grid,
                topo,
                dimension,
                "Coupled-");
        }
    }

    Topology build_topology(Grid &grid, int my_rank, int dimension)
    {
        Topology topo;

        append_physical_patches(grid, topo, my_rank);
        sort_physical_patches_by_priority(grid, topo);

        append_inner_interface_patches(grid, topo, my_rank, dimension);

        auto parallel_bc_all =
            read_all_parallel_boundaries(my_rank, dimension);

        append_parallel_interface_patches(
            grid,
            topo,
            parallel_bc_all,
            my_rank,
            dimension);

        build_derived_topology_patches(grid, topo, dimension);
        detail::build_equivalence(topo, grid, my_rank, dimension);

        if (my_rank == 0)
        {
            std::cout << "*************Finish the Topology Manipulating Process! !**************\n\n";
        }

        return topo;
    }
}
