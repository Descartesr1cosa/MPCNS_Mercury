#include "4_halo/HaloEdgeOwner.h"
#include "0_basic/MPI_WRAPPER.h"

namespace HALO_OWNER
{
    namespace
    {
        constexpr int EDGE_OWNER_SYNC_FLAG_1FORM = 83011;
        constexpr int EDGE_OWNER_SYNC_FLAG_VEC = 83012;

        inline FieldBlock &edge_block(Field &fld,
                                      const IdTriplet &fid,
                                      const TOPO::EntityKey &e)
        {
            return fld.field(fid.at(TOPO::axis_number(e.axis)), e.block);
        }

        inline const FieldBlock &edge_block_const(const Field &fld,
                                                  const IdTriplet &fid,
                                                  const TOPO::EntityKey &e)
        {
            // Field::field 没有 const 重载的话，就别用这个 helper
            return const_cast<Field &>(fld).field(fid.at(TOPO::axis_number(e.axis)), e.block);
        }

        inline int check_edge_triplet_and_get_ncomp(Field &fld, const IdTriplet &fid)
        {
            fid.require_all("edge owner sync");

            const auto &dx = fld.descriptor(fid.xi);
            const auto &dy = fld.descriptor(fid.eta);
            const auto &dz = fld.descriptor(fid.zeta);

            if (dx.ncomp != dy.ncomp || dx.ncomp != dz.ncomp)
            {
                throw std::runtime_error("sync_edge_*: ncomp mismatch among xi/eta/zeta fields.");
            }

            // 这里不强制检查 location；如果你想更严，可以打开
            // if (dx.location != StaggerLocation::EdgeXi ||
            //     dy.location != StaggerLocation::EdgeEt ||
            //     dz.location != StaggerLocation::EdgeZe)
            // {
            //     throw std::runtime_error("sync_edge_*: field locations are not EdgeXi/EdgeEt/EdgeZe.");
            // }

            return dx.ncomp;
        }

        inline void copy_owner_to_rep_local(Field &fld,
                                            const IdTriplet &fid,
                                            const HALO_OWNER::EdgeOwnerLocalAliasItem &it,
                                            int ncomp,
                                            bool use_sign)
        {
            FieldBlock &fb_owner = edge_block(fld, fid, it.owner);
            FieldBlock &fb_rep = edge_block(fld, fid, it.rep);

            const double factor = (use_sign ? static_cast<double>(it.sign) : 1.0);

            for (int m = 0; m < ncomp; ++m)
            {
                fb_rep(it.rep.i, it.rep.j, it.rep.k, m) =
                    factor * fb_owner(it.owner.i, it.owner.j, it.owner.k, m);
            }
        }

        inline void pack_send_buffers(Field &fld,
                                      const IdTriplet &fid,
                                      const HALO_OWNER::EdgeOwnerSyncPattern &pattern,
                                      int ncomp,
                                      std::vector<std::vector<double>> &send_buf)
        {
            const int nrank = static_cast<int>(pattern.send_counts.size());

            if (static_cast<int>(send_buf.size()) != nrank)
                throw std::runtime_error("pack_send_buffers: send_buf_cache size mismatch.");

            for (int r = 0; r < nrank; ++r)
            {
                const int nitem = pattern.send_counts[r];
                const std::size_t need =
                    static_cast<std::size_t>(nitem) * static_cast<std::size_t>(ncomp);

                if (send_buf[r].size() < need)
                    send_buf[r].resize(need);
            }

            for (int r = 0; r < nrank; ++r)
            {
                const int begin = pattern.send_displs[r];
                const int count = pattern.send_counts[r];

                for (int t = 0; t < count; ++t)
                {
                    const auto &it = pattern.send_items[begin + t];
                    FieldBlock &fb_owner = edge_block(fld, fid, it.owner);

                    const int base = t * ncomp;
                    for (int m = 0; m < ncomp; ++m)
                    {
                        send_buf[r][base + m] =
                            fb_owner(it.owner.i, it.owner.j, it.owner.k, m);
                    }
                }
            }
        }

        inline void unpack_recv_buffers(Field &fld,
                                        const IdTriplet &fid,
                                        const HALO_OWNER::EdgeOwnerSyncPattern &pattern,
                                        int ncomp,
                                        const std::vector<std::vector<double>> &recv_buf,
                                        bool use_sign)
        {
            int nrank = static_cast<int>(pattern.recv_counts.size());

            for (int r = 0; r < nrank; ++r)
            {
                const int begin = pattern.recv_displs[r];
                const int count = pattern.recv_counts[r];

                for (int t = 0; t < count; ++t)
                {
                    const auto &it = pattern.recv_items[begin + t];
                    FieldBlock &fb_rep = edge_block(fld, fid, it.rep);

                    const double factor = (use_sign ? static_cast<double>(it.sign) : 1.0);
                    const int base = t * ncomp;

                    for (int m = 0; m < ncomp; ++m)
                    {
                        fb_rep(it.rep.i, it.rep.j, it.rep.k, m) =
                            factor * recv_buf[r][base + m];
                    }
                }
            }
        }

        void sync_edge_common(Field &fld,
                              const IdTriplet &field_id,
                              HALO_OWNER::EdgeOwnerSyncPattern &pattern,
                              bool use_sign,
                              int comm_flag)
        {
            const int ncomp = check_edge_triplet_and_get_ncomp(fld, field_id);

            // ------------------------------------------------------------
            // 1) local alias
            // ------------------------------------------------------------
            for (const auto &it : pattern.local_alias)
            {
                copy_owner_to_rep_local(fld, field_id, it, ncomp, use_sign);
            }

            // ------------------------------------------------------------
            // 2) prepare send / recv buffers
            // ------------------------------------------------------------
            int nrank = 1;
            PARALLEL::mpi_size(&nrank);

            auto &send_buf = pattern.send_buf_cache;
            auto &recv_buf = pattern.recv_buf_cache;

            if (static_cast<int>(recv_buf.size()) != nrank ||
                static_cast<int>(send_buf.size()) != nrank)
            {
                throw std::runtime_error("sync_edge_common: cache rank size mismatch.");
            }

            pack_send_buffers(fld, field_id, pattern, ncomp, send_buf);

            for (int r = 0; r < nrank; ++r)
            {
                const int nitem = pattern.recv_counts[r];
                const std::size_t need =
                    static_cast<std::size_t>(nitem) * static_cast<std::size_t>(ncomp);

                if (recv_buf[r].size() < need)
                    recv_buf[r].resize(need);
            }

            // ------------------------------------------------------------
            // 3) post recv
            // ------------------------------------------------------------
            auto &req_recv = pattern.req_recv_cache;
            auto &stat_recv = pattern.stat_recv_cache;

            int irecv = 0;
            for (int r = 0; r < nrank; ++r)
            {
                const int len = pattern.recv_counts[r] * ncomp;
                if (len <= 0)
                    continue;

                PARALLEL::mpi_data_recv(r, comm_flag, recv_buf[r].data(), len, &req_recv[irecv]);
                ++irecv;
            }

            // ------------------------------------------------------------
            // 4) post send
            // ------------------------------------------------------------
            auto &req_send = pattern.req_send_cache;
            auto &stat_send = pattern.stat_send_cache;

            int isend = 0;
            for (int r = 0; r < nrank; ++r)
            {
                const int len = pattern.send_counts[r] * ncomp;
                if (len <= 0)
                    continue;

                PARALLEL::mpi_data_send(r, comm_flag, send_buf[r].data(), len, &req_send[isend]);
                ++isend;
            }

            // ------------------------------------------------------------
            // 5) wait
            // ------------------------------------------------------------
            int nrecv = irecv;
            int nsend = isend;

            if (nrecv > 0)
                PARALLEL::mpi_wait(nrecv, req_recv.data(), stat_recv.data());
            if (nsend > 0)
                PARALLEL::mpi_wait(nsend, req_send.data(), stat_send.data());

            // PARALLEL::mpi_barrier();

            // ------------------------------------------------------------
            // 6) unpack recv buffers
            // ------------------------------------------------------------
            unpack_recv_buffers(fld, field_id, pattern, ncomp, recv_buf, use_sign);
        }

        //=================================================================
        //=================================================================

        inline int local_owner_index_from_gid(
            const TOPO::Topology &equiv,
            int gid)
        {
            if (gid < equiv.edge_owner_gid_begin || gid >= equiv.edge_owner_gid_end)
            {
                throw std::runtime_error(
                    "local_owner_index_from_gid: gid not owned by this rank.");
            }
            return gid - equiv.edge_owner_gid_begin;
        }

    }

    void gather_local_owner_edges_sorted(
        const TOPO::Topology &equiv,
        std::vector<TOPO::EntityKey> &owner_edges_sorted)
    {
        owner_edges_sorted.clear();
        owner_edges_sorted.reserve(equiv.n_local_edge_owner);

        for (const auto &[e, gid] : equiv.edge_owner_gid)
        {
            owner_edges_sorted.push_back(e);
        }

        std::sort(owner_edges_sorted.begin(), owner_edges_sorted.end(),
                  [&](const TOPO::EntityKey &a, const TOPO::EntityKey &b)
                  {
                      return equiv.edge_owner_gid.at(a) < equiv.edge_owner_gid.at(b);
                  });
    }

    void sync_edge_1form(
        Field &fld,
        const IdTriplet &field_id,
        EdgeOwnerSyncPattern &pattern)
    {
        sync_edge_common(
            fld,
            field_id,
            pattern,
            true, // use sign
            EDGE_OWNER_SYNC_FLAG_1FORM);
    }

    void sync_edge_vec(
        Field &fld,
        const IdTriplet &field_id,
        EdgeOwnerSyncPattern &pattern)
    {
        sync_edge_common(
            fld,
            field_id,
            pattern,
            false, // ignore sign
            EDGE_OWNER_SYNC_FLAG_VEC);
    }

    void pack_owner_edge_1form_local(
        Field &fld,
        const IdTriplet &field_id,
        const TOPO::Topology &equiv,
        const std::vector<TOPO::EntityKey> &owner_edges_sorted,
        std::vector<double> &buf_local)
    {
        const int ncomp = check_edge_triplet_and_get_ncomp(fld, field_id);

        if (equiv.n_local_edge_owner < 0)
        {
            throw std::runtime_error(
                "pack_owner_edge_1form_local: invalid equiv.n_local_edge_owner.");
        }

        if (static_cast<int>(owner_edges_sorted.size()) != equiv.n_local_edge_owner)
        {
            throw std::runtime_error(
                "pack_owner_edge_1form_local: owner edge count mismatch.");
        }

        buf_local.resize(
            static_cast<std::size_t>(equiv.n_local_edge_owner) * ncomp);

        for (int lid = 0; lid < equiv.n_local_edge_owner; ++lid)
        {
            const TOPO::EntityKey &e = owner_edges_sorted[lid];
            FieldBlock &fb = edge_block(fld, field_id, e);

            const int base = lid * ncomp;
            for (int m = 0; m < ncomp; ++m)
            {
                buf_local[base + m] = fb(e.i, e.j, e.k, m);
            }
        }
    }

    void unpack_owner_edge_1form_local(
        const std::vector<double> &buf_local,
        Field &fld,
        const IdTriplet &field_id,
        const TOPO::Topology &equiv,
        const std::vector<TOPO::EntityKey> &owner_edges_sorted,
        EdgeOwnerSyncPattern &pattern)
    {
        const int ncomp = check_edge_triplet_and_get_ncomp(fld, field_id);

        const std::size_t expect_size =
            static_cast<std::size_t>(equiv.n_local_edge_owner) * ncomp;

        if (buf_local.size() != expect_size)
        {
            throw std::runtime_error(
                "unpack_owner_edge_1form_local: buf_local size mismatch.");
        }

        if (static_cast<int>(owner_edges_sorted.size()) != equiv.n_local_edge_owner)
        {
            throw std::runtime_error(
                "unpack_owner_edge_1form_local: owner edge count mismatch.");
        }

        // 1) 写回本 rank owner edges
        for (int lid = 0; lid < equiv.n_local_edge_owner; ++lid)
        {
            const TOPO::EntityKey &e = owner_edges_sorted[lid];
            FieldBlock &fb = edge_block(fld, field_id, e);

            const int base = lid * ncomp;
            for (int m = 0; m < ncomp; ++m)
            {
                fb(e.i, e.j, e.k, m) = buf_local[base + m];
            }
        }

        // 2) owner -> non-owner sync
        sync_edge_1form(fld, field_id, pattern);

        // 3) 普通 halo 先不在这里做，留给外层统一调用
    }

} // namespace HALO_OWNER
