#include "2_topology/TopologyBuildDetail.h"
#include "2_topology/TopologyEquivDetail.h"

#include "0_basic/MPI_WRAPPER.h"
#include "0_basic/Error.h"

#include <algorithm>
#include <sstream>
#include <stdexcept>

namespace TOPO::detail
{
    void pack_node(std::vector<int> &buf, const EntityKey &x)
    {
        buf.push_back(x.rank);
        buf.push_back(x.block);
        buf.push_back(x.i);
        buf.push_back(x.j);
        buf.push_back(x.k);
    }

    EntityKey unpack_node(const int *p)
    {
        return make_node(p[0], p[1], p[2], p[3], p[4]);
    }

    std::vector<int> allgather_packed_records(std::vector<int> send_buf,
                                              int record_size,
                                              const char *context)
    {
        int nrank = 1;
        PARALLEL::mpi_size(&nrank);

        int send_count = static_cast<int>(send_buf.size());
        std::vector<int> recv_counts(nrank, 0);
        std::vector<int> displs(nrank, 0);
        PARALLEL::mpi_gather(&send_count, 1, recv_counts.data(), 1, 0);
        PARALLEL::mpi_bcast(recv_counts.data(), nrank, 0);

        int total_recv = 0;
        for (int r = 0; r < nrank; ++r)
        {
            displs[r] = total_recv;
            total_recv += recv_counts[r];
        }

        std::vector<int> recv_buf(total_recv, 0);
        PARALLEL::mpi_allgatherv(
            send_count > 0 ? send_buf.data() : nullptr,
            send_count,
            total_recv > 0 ? recv_buf.data() : nullptr,
            recv_counts.data(),
            displs.data());

        if (total_recv % record_size != 0)
        {
            std::ostringstream oss;
            oss << context << ": gathered int count is not multiple of " << record_size << ".";
            throw std::runtime_error(oss.str());
        }
        return recv_buf;
    }

    namespace
    {
        EquivMember make_edge_member(const EntityKey &e, int orient_sign, bool is_owner)
        {
            EquivMember m;
            m.entity = e;
            m.orient_sign = orient_sign;
            m.is_owner = is_owner;
            return m;
        }

        EquivMember make_face_member(const EntityKey &f, int orient_sign, bool is_owner)
        {
            EquivMember m;
            m.entity = f;
            m.orient_sign = orient_sign;
            m.is_owner = is_owner;
            return m;
        }
    }

    void clear_equivalence(Topology &equiv)
    {
        equiv.nodes.local_to_rep.clear();
        equiv.nodes.rep_to_qid.clear();
        equiv.nodes.rep_count.clear();
        equiv.edges.local_to_qkey.clear();
        equiv.edges.local_to_qsign.clear();
        equiv.edges.qkey_to_members.clear();
        equiv.edges.qkey_to_owner.clear();
        equiv.edges.local_is_owner.clear();

        equiv.edges.owner_to_gid.clear();
        equiv.edges.gid_to_owner.clear();
        equiv.edges.qkey_to_qid.clear();
        equiv.edges.n_local_owner = 0;
        equiv.edges.n_global_owner = 0;
        equiv.edges.owner_gid_begin = 0;
        equiv.edges.owner_gid_end = 0;

        equiv.faces.local_to_qkey.clear();
        equiv.faces.local_to_qsign.clear();
        equiv.faces.qkey_to_members.clear();
        equiv.faces.qkey_to_owner.clear();
        equiv.faces.local_is_owner.clear();

        equiv.faces.owner_to_gid.clear();
        equiv.faces.gid_to_owner.clear();
        equiv.faces.qkey_to_qid.clear();
        equiv.faces.n_local_owner = 0;
        equiv.faces.n_global_owner = 0;
        equiv.faces.owner_gid_begin = 0;
        equiv.faces.owner_gid_end = 0;

        equiv.cells.local_to_qid.clear();

        equiv.nodes.classes.clear();
        equiv.edges.classes.clear();
        equiv.faces.classes.clear();
    }

    void rebuild_edge_classes(Topology &equiv)
    {
        equiv.edges.classes.clear();
        equiv.edges.classes.reserve(equiv.edges.qkey_to_members.size());

        for (const auto &[key, members] : equiv.edges.qkey_to_members)
        {
            EquivClass cls;
            cls.dim = EntityDim::Edge;

            auto owner_it = equiv.edges.qkey_to_owner.find(key);
            const bool has_owner = (owner_it != equiv.edges.qkey_to_owner.end());
            EntityKey owner{};

            if (has_owner)
            {
                owner = owner_it->second;
                auto gid_it = equiv.edges.owner_to_gid.find(owner);
                if (gid_it != equiv.edges.owner_to_gid.end())
                    cls.global_id = gid_it->second;

                int owner_sign = +1;
                auto sign_it = equiv.edges.local_to_qsign.find(owner);
                if (sign_it != equiv.edges.local_to_qsign.end())
                    owner_sign = static_cast<int>(sign_it->second);

                cls.owner = make_edge_member(owner, owner_sign, true);
            }

            cls.members.reserve(members.size());
            for (const auto &e : members)
            {
                int orient_sign = +1;
                auto sign_it = equiv.edges.local_to_qsign.find(e);
                if (sign_it != equiv.edges.local_to_qsign.end())
                    orient_sign = static_cast<int>(sign_it->second);

                bool is_owner = false;
                auto owner_flag_it = equiv.edges.local_is_owner.find(e);
                if (owner_flag_it != equiv.edges.local_is_owner.end())
                    is_owner = owner_flag_it->second;
                else if (has_owner)
                    is_owner = (e == owner);

                cls.members.push_back(make_edge_member(e, orient_sign, is_owner));
            }

            equiv.edges.classes.push_back(cls);
        }
    }

    void rebuild_face_classes(Topology &equiv)
    {
        equiv.faces.classes.clear();
        equiv.faces.classes.reserve(equiv.faces.qkey_to_members.size());

        for (const auto &[key, members] : equiv.faces.qkey_to_members)
        {
            EquivClass cls;
            cls.dim = EntityDim::Face;

            auto owner_it = equiv.faces.qkey_to_owner.find(key);
            const bool has_owner = (owner_it != equiv.faces.qkey_to_owner.end());
            EntityKey owner{};

            if (has_owner)
            {
                owner = owner_it->second;
                auto gid_it = equiv.faces.owner_to_gid.find(owner);
                if (gid_it != equiv.faces.owner_to_gid.end())
                    cls.global_id = gid_it->second;

                int owner_sign = +1;
                auto sign_it = equiv.faces.local_to_qsign.find(owner);
                if (sign_it != equiv.faces.local_to_qsign.end())
                    owner_sign = static_cast<int>(sign_it->second);

                cls.owner = make_face_member(owner, owner_sign, true);
            }

            cls.members.reserve(members.size());
            for (const auto &f : members)
            {
                int orient_sign = +1;
                auto sign_it = equiv.faces.local_to_qsign.find(f);
                if (sign_it != equiv.faces.local_to_qsign.end())
                    orient_sign = static_cast<int>(sign_it->second);

                bool is_owner = false;
                auto owner_flag_it = equiv.faces.local_is_owner.find(f);
                if (owner_flag_it != equiv.faces.local_is_owner.end())
                    is_owner = owner_flag_it->second;
                else if (has_owner)
                    is_owner = (f == owner);

                cls.members.push_back(make_face_member(f, orient_sign, is_owner));
            }

            equiv.faces.classes.push_back(cls);
        }
    }

    void build_equivalence(
        Topology &equiv,
        Grid &grid,
        int my_rank,
        int dimension)
    {
        clear_equivalence(equiv);
        equiv.nodes.rep_count.clear();

        auto all_local_nodes = collect_all_local_nodes(grid, my_rank);
        reconcile_node_equivalence_parallel(equiv, my_rank, all_local_nodes, equiv, equiv.nodes.rep_count);
        build_node_entity_ids(equiv);

        EdgeBuildScratch edge_scratch;
        auto all_local_edges = collect_all_local_edges(grid, my_rank, dimension);
        build_edge_equivalence_from_nodes(all_local_edges, equiv, edge_scratch);
        build_edge_entity_ids(equiv);

        select_edge_owner_parallel(equiv, edge_scratch, equiv.nodes.rep_count);
        build_edge_owner_gid(my_rank, equiv);
        rebuild_edge_classes(equiv);

        FaceBuildScratch face_scratch;
        auto all_local_faces = collect_all_local_faces(grid, my_rank, dimension);
        build_face_equivalence_from_nodes(all_local_faces, equiv, face_scratch);
        build_face_entity_ids(equiv);

        select_face_owner_parallel(equiv, face_scratch, equiv.nodes.rep_count);
        build_face_owner_gid(my_rank, equiv);
        rebuild_face_classes(equiv);

        auto all_local_cells = collect_all_local_cells(grid, my_rank, dimension);
        build_cell_entity_ids(all_local_cells, equiv);
    }
} // namespace TOPO::detail
