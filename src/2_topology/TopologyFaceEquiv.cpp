#include "2_topology/TopologyEquivDetail.h"
#include "2_topology/LocalIncidence.h"

#include "0_basic/MPI_WRAPPER.h"
#include "1_grid/1_MPCNS_Grid.h"

#include <algorithm>
#include <map>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace TOPO::detail
{
    using FaceEdgeStencil = std::map<EdgeKey, int>;
    constexpr int kFaceStencilMaxEntries = 4;
    constexpr int kPackedFaceStencilEntrySize = 11;
    constexpr int kPackedFaceOwnerCandidateSize =
        20 + 6 + 1 + kFaceStencilMaxEntries * kPackedFaceStencilEntrySize;

    static int node_equiv_count(const std::unordered_map<EntityKey, int, EntityKey::Hash> &rep_count,
                                const EntityKey &node)
    {
        auto it = rep_count.find(node);
        return (it != rep_count.end()) ? it->second : 1;
    }

    std::vector<EntityKey> collect_all_local_faces(
        Grid &grid,
        int my_rank,
        int dimension)
    {
        std::vector<EntityKey> faces;

        for (int ib = 0; ib < grid.nblock; ++ib)
        {
        const auto &blk = grid.grids(ib);

        if (dimension >= 3)
        {
            // FaceXi
            for (int i = 0; i <= blk.mx; ++i)
            for (int j = 0; j < blk.my; ++j)
                for (int k = 0; k < blk.mz; ++k)
                {
                faces.push_back(make_face(my_rank, ib, i, j, k, EntityAxis::Xi));
                }

            // FaceEt
            for (int i = 0; i < blk.mx; ++i)
            for (int j = 0; j <= blk.my; ++j)
                for (int k = 0; k < blk.mz; ++k)
                {
                faces.push_back(make_face(my_rank, ib, i, j, k, EntityAxis::Eta));
                }

            // FaceZe
            for (int i = 0; i < blk.mx; ++i)
            for (int j = 0; j < blk.my; ++j)
                for (int k = 0; k <= blk.mz; ++k)
                {
                faces.push_back(make_face(my_rank, ib, i, j, k, EntityAxis::Zeta));
                }
        }
        else if (dimension >= 2)
        {
            // 2D uses xi/eta face locations as line faces in the single k plane.
            const int k = 0;

            for (int i = 0; i <= blk.mx; ++i)
            for (int j = 0; j < blk.my; ++j)
            {
                faces.push_back(make_face(my_rank, ib, i, j, k, EntityAxis::Xi));
            }

            for (int i = 0; i < blk.mx; ++i)
            for (int j = 0; j <= blk.my; ++j)
            {
                faces.push_back(make_face(my_rank, ib, i, j, k, EntityAxis::Eta));
            }
        }
        }
        return faces;
    }

    
    bool face_key_is_shared(
        const FaceKey &key,
        const std::unordered_map<EntityKey, int, EntityKey::Hash> &rep_count)
    {
        return node_equiv_count(rep_count, key.a) > 1 ||
           node_equiv_count(rep_count, key.b) > 1 ||
           node_equiv_count(rep_count, key.c) > 1 ||
           node_equiv_count(rep_count, key.d) > 1;
    }

    
    void build_face_equivalence_from_nodes(
        const std::vector<EntityKey> &all_local_faces,
        Topology &equiv,
        FaceBuildScratch &scratch)
    {
        equiv.faces.local_to_qkey.clear();
        equiv.faces.local_to_qsign.clear();
        scratch.qkey_to_local_members.clear();

        for (const auto &f : all_local_faces)
        {
        int8_t sign_to_canonical = 0;
        FaceKey key = make_face_key(f, equiv.nodes.local_to_rep, sign_to_canonical);

        equiv.faces.local_to_qkey[f] = key;
        equiv.faces.local_to_qsign[f] = sign_to_canonical;
        scratch.qkey_to_local_members[key].push_back(f);
        }
    }

    
    void pack_face_local(std::vector<int> &buf, const EntityKey &f)
    {
        buf.push_back(f.rank);
        buf.push_back(f.block);
        buf.push_back(f.i);
        buf.push_back(f.j);
        buf.push_back(f.k);
        buf.push_back(axis_number(f.axis));
    }

    
    EntityKey unpack_face_local(const int *p)
    {
        return make_face(p[0], p[1], p[2], p[3], p[4], entity_axis(p[5]));
    }

    
    void pack_face_key(std::vector<int> &buf, const FaceKey &key)
    {
        pack_node(buf, key.a);
        pack_node(buf, key.b);
        pack_node(buf, key.c);
        pack_node(buf, key.d);
    }

    
    FaceKey unpack_face_key(const int *p)
    {
        EntityKey a = unpack_node(p + 0);
        EntityKey b = unpack_node(p + 5);
        EntityKey c = unpack_node(p + 10);
        EntityKey d = unpack_node(p + 15);
        return FaceKey{a, b, c, d};
    }

    
    void build_face_entity_ids(Topology &equiv)
    {
        std::vector<FaceKey> local_keys;
        local_keys.reserve(equiv.faces.local_to_qkey.size());
        for (const auto &[face, key] : equiv.faces.local_to_qkey)
        {
        (void)face;
        local_keys.push_back(key);
        }
        std::sort(local_keys.begin(), local_keys.end());
        local_keys.erase(std::unique(local_keys.begin(), local_keys.end()), local_keys.end());

        std::vector<int> send_buf;
        send_buf.reserve(local_keys.size() * 20);
        for (const FaceKey &key : local_keys)
        pack_face_key(send_buf, key);

        const std::vector<int> recv_buf =
        allgather_packed_records(send_buf, 20, "build_topology face EntityId");
        std::vector<FaceKey> global_keys;
        global_keys.reserve(recv_buf.size() / 20);
        for (std::size_t n = 0; n < recv_buf.size(); n += 20)
        global_keys.push_back(unpack_face_key(recv_buf.data() + n));

        std::sort(global_keys.begin(), global_keys.end());
        global_keys.erase(std::unique(global_keys.begin(), global_keys.end()), global_keys.end());
        equiv.faces.qkey_to_qid.clear();
        for (std::size_t n = 0; n < global_keys.size(); ++n)
        equiv.faces.qkey_to_qid[global_keys[n]] = static_cast<int>(n);
    }

    
    void add_stencil_term(FaceEdgeStencil &stencil, const EdgeKey &edge, int sign)
    {
        const int coefficient = (stencil[edge] += sign);
        if (coefficient == 0)
        stencil.erase(edge);
    }

    
    FaceEdgeStencil face_boundary_stencil_to_canonical_edges(
        const Topology &equiv,
        const EntityKey &face,
        const char *context)
    {
        FaceEdgeStencil stencil;
        for (const IncidenceEntry &entry : boundary_of_face(face))
        {
        int8_t edge_sign = 0;
        EdgeKey edge_key{};
        try
        {
            edge_key = make_edge_key(entry.entity, equiv.nodes.local_to_rep, edge_sign);
        }
        catch (const std::exception &error)
        {
            std::ostringstream oss;
            oss << context << ": cannot build boundary edge key"
            << " face=(" << face.rank << "," << face.block << ","
            << face.i << "," << face.j << "," << face.k
            << ",axis=" << axis_number(face.axis) << ")"
            << " edge=(" << entry.entity.rank << "," << entry.entity.block << ","
            << entry.entity.i << "," << entry.entity.j << "," << entry.entity.k
            << ",axis=" << axis_number(entry.entity.axis) << ")"
            << " reason=" << error.what();
            throw std::runtime_error(oss.str());
        }

        add_stencil_term(
            stencil,
            edge_key,
            entry.sign * static_cast<int>(edge_sign));
        }
        return stencil;
    }

    
    void pack_face_stencil(std::vector<int> &buf, const FaceEdgeStencil &stencil)
    {
        if (stencil.size() > kFaceStencilMaxEntries)
        throw std::runtime_error("build_topology: face stencil has more than four edge entries.");

        buf.push_back(static_cast<int>(stencil.size()));

        int packed_entries = 0;
        for (const auto &[edge_key, coefficient] : stencil)
        {
        pack_edge_key(buf, edge_key);
        buf.push_back(coefficient);
        ++packed_entries;
        }

        for (; packed_entries < kFaceStencilMaxEntries; ++packed_entries)
        {
        for (int n = 0; n < kPackedFaceStencilEntrySize; ++n)
            buf.push_back(0);
        }
    }

    
    FaceEdgeStencil unpack_face_stencil(const int *p)
    {
        const int count = p[0];
        if (count < 0 || count > kFaceStencilMaxEntries)
        throw std::runtime_error("build_topology: packed face stencil entry count is invalid.");

        FaceEdgeStencil stencil;
        const int *entry = p + 1;
        for (int n = 0; n < count; ++n)
        {
        const EdgeKey edge_key = unpack_edge_key(entry);
        const int coefficient = entry[10];
        add_stencil_term(stencil, edge_key, coefficient);
        entry += kPackedFaceStencilEntrySize;
        }
        return stencil;
    }

    
    void pack_face_owner_candidate(
        std::vector<int> &buf,
        const FaceKey &key,
        const EntityKey &f,
        const FaceEdgeStencil &stencil)
    {
        // FaceKey(20) + EntityKey(6) + stencil_count(1)
        // + up to four {EdgeKey(10), coefficient(1)} entries = 71 ints.
        pack_face_key(buf, key);
        pack_face_local(buf, f);
        pack_face_stencil(buf, stencil);
    }

    
    bool opposite_stencil(const FaceEdgeStencil &lhs, const FaceEdgeStencil &rhs)
    {
        if (lhs.size() != rhs.size())
        return false;

        for (const auto &[edge, coefficient] : lhs)
        {
        const auto it = rhs.find(edge);
        if (it == rhs.end() || it->second != -coefficient)
            return false;
        }
        return true;
    }

    
    void orient_face_members_to_owner(
        Topology &equiv,
        const FaceKey &key,
        const EntityKey &owner,
        const std::vector<EntityKey> &members,
        const std::unordered_map<EntityKey, FaceEdgeStencil, EntityKey::Hash> &stencils)
    {
        const auto owner_stencil_it = stencils.find(owner);
        if (owner_stencil_it == stencils.end())
        throw std::runtime_error("build_topology: owner face stencil is missing.");
        const FaceEdgeStencil &owner_stencil = owner_stencil_it->second;

        equiv.faces.local_to_qsign[owner] = +1;

        for (const EntityKey &member : members)
        {
        if (member == owner)
            continue;

        const auto member_stencil_it = stencils.find(member);
        if (member_stencil_it == stencils.end())
            throw std::runtime_error("build_topology: member face stencil is missing.");
        const FaceEdgeStencil &member_stencil = member_stencil_it->second;

        if (member_stencil == owner_stencil)
        {
            equiv.faces.local_to_qsign[member] = +1;
        }
        else if (opposite_stencil(member_stencil, owner_stencil))
        {
            equiv.faces.local_to_qsign[member] = -1;
        }
        else
        {
            std::ostringstream oss;
            oss << "build_topology: face member boundary does not match owner"
            << " key_first_corner=(" << key.a.rank << "," << key.a.block << ","
            << key.a.i << "," << key.a.j << "," << key.a.k << ")"
            << " owner=(" << owner.rank << "," << owner.block << ","
            << owner.i << "," << owner.j << "," << owner.k
            << ",axis=" << axis_number(owner.axis) << ")"
            << " member=(" << member.rank << "," << member.block << ","
            << member.i << "," << member.j << "," << member.k
            << ",axis=" << axis_number(member.axis) << ")";
            throw std::runtime_error(oss.str());
        }
        }
    }

    
    void select_face_owner_parallel(
        Topology &equiv,
        const FaceBuildScratch &scratch,
        const std::unordered_map<EntityKey, int, EntityKey::Hash> &rep_count)
    {
        std::vector<int> send_buf;
        send_buf.reserve(1024);

        for (const auto &[key, members] : scratch.qkey_to_local_members)
        {
        if (!face_key_is_shared(key, rep_count) && members.size() <= 1)
            continue;

        for (const auto &f : members)
        {
            const FaceEdgeStencil stencil =
            face_boundary_stencil_to_canonical_edges(
                equiv,
                f,
                "build_topology: pack face owner candidate");
            pack_face_owner_candidate(send_buf, key, f, stencil);
        }
        }

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

        if (total_recv % kPackedFaceOwnerCandidateSize != 0)
        {
        throw std::runtime_error(
            "build_topology: gathered face-candidate int count has invalid record size.");
        }

        std::unordered_map<FaceKey, EntityKey, FaceKey::Hash> global_owner;
        std::unordered_map<FaceKey, std::vector<EntityKey>, FaceKey::Hash> global_members;
        std::unordered_map<EntityKey, FaceEdgeStencil, EntityKey::Hash> global_stencils;

        const int ncand = total_recv / kPackedFaceOwnerCandidateSize;
        for (int c = 0; c < ncand; ++c)
        {
        const int *base = recv_buf.data() + kPackedFaceOwnerCandidateSize * c;

        FaceKey key = unpack_face_key(base + 0);
        EntityKey f = unpack_face_local(base + 20);
        FaceEdgeStencil stencil = unpack_face_stencil(base + 26);

        equiv.faces.local_to_qkey[f] = key;
        global_members[key].push_back(f);
        global_stencils[f] = std::move(stencil);

        auto it = global_owner.find(key);
        if (it == global_owner.end() || f < it->second)
        {
            global_owner[key] = f;
        }
        }

        equiv.faces.qkey_to_owner.clear();
        equiv.faces.local_is_owner.clear();
        equiv.faces.qkey_to_members.clear();

        for (auto &[key, members] : global_members)
        {
        std::sort(members.begin(), members.end());

        auto it = global_owner.find(key);
        if (it == global_owner.end())
        {
            throw std::runtime_error(
            "build_topology: local face key missing in global_owner.");
        }

        const EntityKey &owner = it->second;
        orient_face_members_to_owner(equiv, key, owner, members, global_stencils);

        equiv.faces.qkey_to_owner[key] = owner;
        equiv.faces.qkey_to_members[key] = members;

        for (const auto &f : members)
        {
            equiv.faces.local_is_owner[f] = (f == owner);
        }
        }

    }

    
    void build_face_owner_gid(int my_rank, Topology &equiv)
    {
        std::vector<EntityKey> local_owner_faces;
        local_owner_faces.reserve(equiv.faces.local_is_owner.size());

        for (const auto &[f, is_owner] : equiv.faces.local_is_owner)
        {
        if (is_owner && f.rank == my_rank)
            local_owner_faces.push_back(f);
        }

        std::sort(local_owner_faces.begin(), local_owner_faces.end());

        equiv.faces.n_local_owner = static_cast<int>(local_owner_faces.size());

        int nrank = 1;
        PARALLEL::mpi_size(&nrank);

        std::vector<int> counts(nrank, 0);
        PARALLEL::mpi_gather(&equiv.faces.n_local_owner, 1, counts.data(), 1, 0);
        PARALLEL::mpi_bcast(counts.data(), nrank, 0);

        equiv.faces.owner_gid_begin = 0;
        for (int r = 0; r < my_rank; ++r)
        {
        equiv.faces.owner_gid_begin += counts[r];
        }

        equiv.faces.owner_gid_end = equiv.faces.owner_gid_begin + equiv.faces.n_local_owner;

        equiv.faces.n_global_owner = 0;
        for (int r = 0; r < nrank; ++r)
        {
        equiv.faces.n_global_owner += counts[r];
        }

        equiv.faces.owner_to_gid.clear();
        equiv.faces.gid_to_owner.clear();

        std::vector<int> send_buf;
        send_buf.reserve(local_owner_faces.size() * 7);

        for (int n = 0; n < equiv.faces.n_local_owner; ++n)
        {
        const EntityKey &f = local_owner_faces[n];
        int gid = equiv.faces.owner_gid_begin + n;

        pack_face_local(send_buf, f);
        send_buf.push_back(gid);
        }

        std::vector<int> recv_counts(nrank, 0);
        std::vector<int> displs(nrank, 0);
        int send_count = static_cast<int>(send_buf.size());

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

        if (total_recv % 7 != 0)
        {
        throw std::runtime_error(
            "build_topology: gathered face-owner-gid int count is not multiple of 7.");
        }

        const int nowners = total_recv / 7;
        for (int n = 0; n < nowners; ++n)
        {
        const int *base = recv_buf.data() + 7 * n;
        EntityKey f = unpack_face_local(base);
        int gid = base[6];

        equiv.faces.owner_to_gid[f] = gid;
        equiv.faces.gid_to_owner[gid] = f;
        }
    }
    } // namespace TOPO::detail
