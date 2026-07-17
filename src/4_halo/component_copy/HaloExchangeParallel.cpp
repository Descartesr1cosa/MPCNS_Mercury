#include "4_halo/Halo.h"
#include "0_basic/MPI_WRAPPER.h"
#include "0_basic/Error.h"

void Halo::exchange_parallel(const std::string &field_name)
{
    //=========================================================================
    // 1. 找到这个 field 的 descriptor 和 fid
    const FieldDescriptor *desc = nullptr;
    int fid = -1;

    for (int id = 0; id < fld_->num_fields(); ++id)
    {
        if (field_name == fld_->descriptor(id).name)
        {
            desc = &(fld_->descriptor(id));
            fid = id;
            break;
        }
    }
    if (desc == nullptr)
    {
        std::cout << "Fatal Error!!! Can not find the field:\t" << field_name << std::endl;
        std::exit(-1);
    }

    PatternKey key{desc->location, desc->nghost};

    //=========================================================================
    // 2. 找到对应的 parallel pattern
    auto it = parallel_patterns_.find(key);
    if (it == parallel_patterns_.end())
    {
        std::cout << "Fatal Error!!! Can not find Parallel Halo pattern of field:\t"
                  << field_name << std::endl;
        std::exit(-1);
    }

    const HaloPattern &pat = it->second;
    const int ncomp = desc->ncomp;

    //=========================================================================
    // 3. 遍历所有 Parallel HaloRegion：每个 region 对应一次 “本块 <-> 邻居 rank” 的通信，打包
    //-------------------------------------------------------------------------
    // 检测缓冲空间是否足够
    const int num_face = pat.regions.size();
    std::vector<int> local_active(num_face, 0), peer_active(num_face, 0);
    std::vector<MPI_Request> active_send_req(num_face, MPI_REQUEST_NULL);
    std::vector<MPI_Request> active_recv_req(num_face, MPI_REQUEST_NULL);
    for (int ir = 0; ir < num_face; ++ir)
    {
        const HaloRegion &r = pat.regions[ir];
        local_active[ir] = fld_->field(fid, r.this_block).is_allocated() ? 1 : 0;
        MPI_Irecv(&peer_active[ir], 1, MPI_INT, r.neighbor_rank, r.recv_flag,
                  MPI_COMM_WORLD, &active_recv_req[ir]);
        MPI_Isend(&local_active[ir], 1, MPI_INT, r.neighbor_rank, r.send_flag,
                  MPI_COMM_WORLD, &active_send_req[ir]);
    }
    if (num_face > 0)
    {
        MPI_Waitall(num_face, active_send_req.data(), MPI_STATUSES_IGNORE);
        MPI_Waitall(num_face, active_recv_req.data(), MPI_STATUSES_IGNORE);
    }
    if (send_buf.size() < num_face)
        send_buf.resize(num_face);
    if (recv_buf.size() < num_face)
        recv_buf.resize(num_face);
    if (req_send.size() < num_face)
        req_send.resize(num_face);
    if (req_recv.size() < num_face)
        req_recv.resize(num_face);
    if (stat_send.size() < num_face)
        stat_send.resize(num_face);
    if (stat_recv.size() < num_face)
        stat_recv.resize(num_face);
    if (length.size() < num_face)
        length.resize(num_face);
    std::vector<int32_t> recv_length(num_face, 0);
    //-------------------------------------------------------------------------
    // 打包
    int index = 0;
    for (const HaloRegion &r : pat.regions)
    {
        if (!local_active[index] || !peer_active[index])
        {
            length[index] = 0;
            recv_length[index] = 0;
            ++index;
            continue;
        }
        FieldBlock &fb = fld_->field(fid, r.this_block); // 本 rank 上的块

        const Box3 &sb = r.send_box; // 本块 inner strip
        const Box3 &rb = r.recv_box; // 本块 ghost strip

        const int32_t send_total = (sb.hi.i - sb.lo.i) *
                                   (sb.hi.j - sb.lo.j) *
                                   (sb.hi.k - sb.lo.k) *
                                   ncomp;
        const int32_t recv_total = (rb.hi.i - rb.lo.i) *
                                   (rb.hi.j - rb.lo.j) *
                                   (rb.hi.k - rb.lo.k) *
                                   ncomp;
        length[index] = send_total;
        recv_length[index] = recv_total;
        // 4. 确保缓冲区足够大（复用 send_buf_ / recv_buf_）
        if (send_buf[index].size() < static_cast<std::size_t>(send_total))
            send_buf[index].resize(send_total);
        if (recv_buf[index].size() < static_cast<std::size_t>(recv_total))
            recv_buf[index].resize(recv_total);

        // 5. 打包：本块 inner strip -> send_buf_

        const TOPO::IndexTransform &transform = r.trans;

        // 1) 本地 send_box 的尺寸
        int loc_lo[3] = {sb.lo.i, sb.lo.j, sb.lo.k};
        int loc_hi[3] = {sb.hi.i - 1, sb.hi.j - 1, sb.hi.k - 1}; // 形成闭区间
        int len_loc[3] = {sb.hi.i - sb.lo.i, sb.hi.j - sb.lo.j, sb.hi.k - sb.lo.k};

        // 2) 获取目标块对应区域的范围，这里暂时用闭区间，主要获得最小值neighbor的lo.ijk, 存入tar_ref
        // tar_ref作为起点 用于编码
        int offset[3] = {r.trans.offset.i,
                         r.trans.offset.j,
                         r.trans.offset.k};
        int tar1[3], tar2[3], tar_ref[3];
        for (int d = 0; d < 3; ++d)
            tar1[transform.perm[d]] = transform.sign[d] * loc_lo[d] + offset[d];
        for (int d = 0; d < 3; ++d)
            tar2[transform.perm[d]] = transform.sign[d] * loc_hi[d] + offset[d];

        tar_ref[0] = (tar1[0] <= tar2[0]) ? tar1[0] : tar2[0];
        tar_ref[1] = (tar1[1] <= tar2[1]) ? tar1[1] : tar2[1];
        tar_ref[2] = (tar1[2] <= tar2[2]) ? tar1[2] : tar2[2];

        // 3) 对应到“邻居坐标系”下的尺寸：
        //    邻居 axis_nb = perm[d] 这一维的长度就是本地 d 维的长度
        int len_nb[3] = {0, 0, 0};
        for (int d = 0; d < 3; ++d)
            len_nb[transform.perm[d]] = len_loc[d];

        // 4) 真正打包：本地 (i,j,k) -> 邻居坐标 nb -> 相对坐标 (ri,rj,rk) -> buffer index
        int ijk[3], tar_ijk[3];
        for (ijk[0] = sb.lo.i; ijk[0] < sb.hi.i; ijk[0]++)
            for (ijk[1] = sb.lo.j; ijk[1] < sb.hi.j; ijk[1]++)
                for (ijk[2] = sb.lo.k; ijk[2] < sb.hi.k; ijk[2]++)
                {
                    for (int ii = 0; ii < 3; ii++)
                        tar_ijk[transform.perm[ii]] = transform.sign[ii] * ijk[ii] + offset[ii];
                    // 获得了目标tar_ijk[]坐标，现在设定编码方式：
                    //  按for i j k comp的顺序排列，以较小值tar_ref为起点
                    // 以 nb_ref 为原点的“邻居相对坐标”
                    int ri = tar_ijk[0] - tar_ref[0];
                    int rj = tar_ijk[1] - tar_ref[1];
                    int rk = tar_ijk[2] - tar_ref[2];
                    int32_t base = ((ri * len_nb[1] + rj) * len_nb[2] + rk) * ncomp;
                    // 真正的数据拷贝
                    for (int m = 0; m < ncomp; ++m)
                        send_buf[index][base + m] = fb(ijk[0], ijk[1], ijk[2], m);
                }

        // 5) 记录regions的个数
        index++;
    }
    //-------------------------------------------------------------------------
    // 等待
    // PARALLEL::mpi_barrier();
    //-------------------------------------------------------------------------
    // Post receives before sends.  The pattern fixes the message sizes, so a
    // blocking Probe here only serialized otherwise independent interfaces.
    // MPI_Get_count after Waitall retains the nonconforming-pattern check.
    std::vector<int32_t> actual_recv_length(num_face, 0);
    index = 0;
    for (const HaloRegion &r : pat.regions)
    {
        if (local_active[index] && peer_active[index])
            PARALLEL::mpi_data_recv(r.neighbor_rank, r.recv_flag,
                                    recv_buf[index].data(), recv_length[index], &(req_recv[index]));
        else
            req_recv[index] = MPI_REQUEST_NULL;
        index++;
    }
    index = 0;
    for (const HaloRegion &r : pat.regions)
    {
        if (local_active[index] && peer_active[index])
            PARALLEL::mpi_data_send(r.neighbor_rank, r.send_flag,
                                    send_buf[index].data(), length[index], &(req_send[index]));
        else
            req_send[index] = MPI_REQUEST_NULL;
        index++;
    }
    //----------------------------------------------------------------------
    // 等待完成
    int num_face_comm = num_face;
    PARALLEL::mpi_wait(num_face_comm, req_send.data(), stat_send.data());
    PARALLEL::mpi_wait(num_face_comm, req_recv.data(), stat_recv.data());

    for (int ir = 0; ir < num_face; ++ir)
    {
        if (local_active[ir] && peer_active[ir])
        {
            int incoming = 0;
            MPI_Get_count(&stat_recv[ir], MPI_DOUBLE, &incoming);
            actual_recv_length[ir] = incoming;
        }
    }

    bool length_mismatch = false;
    index = 0;
    for (const HaloRegion &r : pat.regions)
    {
        if (!local_active[index] || !peer_active[index])
        {
            ++index;
            continue;
        }
        if (actual_recv_length[index] != recv_length[index])
        {
            int myid = 0;
            PARALLEL::mpi_rank(&myid);
            const Box3 &rb = r.recv_box;
            std::cerr << "[Halo] nonconforming parallel face: field=" << field_name
                      << " rank=" << myid << " block=" << r.this_block
                      << " peer=" << r.neighbor_rank << " recv_tag=" << r.recv_flag
                      << " expected=" << recv_length[index]
                      << " incoming=" << actual_recv_length[index]
                      << " recv_box=[" << rb.lo.i << ',' << rb.lo.j << ',' << rb.lo.k
                      << "]-[" << rb.hi.i << ',' << rb.hi.j << ',' << rb.hi.k << "]\n";
            length_mismatch = true;
        }
        ++index;
    }
    if (length_mismatch)
        ERROR::Abort("Halo::exchange_parallel: nonconforming interface message length");
    //----------------------------------------------------------------------

    //----------------------------------------------------------------------
    // 6. 解包：recv_buf -> 本块的 recv_box（ghost cells）
    //----------------------------------------------------------------------
    index = 0;
    for (const HaloRegion &r : pat.regions)
    {
        if (!local_active[index] || !peer_active[index])
        {
            ++index;
            continue;
        }
        FieldBlock &fb = fld_->field(fid, r.this_block); // 本 rank 上的块
        const Box3 &rb = r.recv_box;                     // 要填充的 halo 区域

        const int ni = rb.hi.i - rb.lo.i;
        const int nj = rb.hi.j - rb.lo.j;
        const int nk = rb.hi.k - rb.lo.k;

#if if_Debug_Field_Array == 1
        // 防御式检查：与打包时记录的长度一致

        const int32_t n_total = ni * nj * nk * ncomp;
        if (n_total != recv_length[index])
        {
            std::cout << "Fatal Error!!! Parallel Halo unpack n_total mismatch "
                      << "(field=" << field_name << ", block=" << r.this_block << ")\n";
            std::exit(-1);
        }
#endif

        const std::vector<double> &buf = recv_buf[index];

        // 这里我们假设发送端按照“邻居坐标系（也就是当前块坐标系）”
        // 以 comp , k 最快、再 j、再 i 的顺序排布：
        //
        // base = (((i_rel) * nj + j_rel) * nk + k_rel) * ncomp
        //
        // 其中 i_rel = i - rb.lo.i，j_rel = j - rb.lo.j，k_rel = k - rb.lo.k
        for (int ii = 0; ii < ni; ++ii)
        {
            int i_loc = rb.lo.i + ii;
            for (int jj = 0; jj < nj; ++jj)
            {
                int j_loc = rb.lo.j + jj;
                for (int kk = 0; kk < nk; ++kk)
                {
                    int k_loc = rb.lo.k + kk;

                    int32_t base = (((ii * nj + jj) * nk) + kk) * ncomp;

                    for (int m = 0; m < ncomp; ++m)
                    {
                        fb(i_loc, j_loc, k_loc, m) = buf[base + m];
                    }
                }
            }
        }

        ++index;
    }
}

void Halo::exchange_parallel_vertex(const std::string &field_name)
{
    //=========================================================================
    // 1. 找到这个 field 的 descriptor 和 fid
    const FieldDescriptor *desc = nullptr;
    int fid = -1;

    for (int id = 0; id < fld_->num_fields(); ++id)
    {
        if (field_name == fld_->descriptor(id).name)
        {
            desc = &(fld_->descriptor(id));
            fid = id;
            break;
        }
    }
    if (desc == nullptr)
    {
        std::cout << "Fatal Error!!! Can not find the field:\t" << field_name << std::endl;
        std::exit(-1);
    }

    PatternKey key{desc->location, desc->nghost};

    //=========================================================================
    // 2. 找到对应的 parallel pattern
    auto it_send = parallel_vertex_patterns_send.find(key);
    if (it_send == parallel_vertex_patterns_send.end())
    {
        std::cout << "Fatal Error!!! Can not find Parallel Halo pattern of field:\t"
                  << field_name << std::endl;
        std::exit(-1);
    }

    auto it_recv = parallel_vertex_patterns_recv.find(key);
    if (it_recv == parallel_vertex_patterns_recv.end())
    {
        std::cout << "Fatal Error!!! Can not find Parallel Halo pattern of field:\t"
                  << field_name << std::endl;
        std::exit(-1);
    }

    const HaloPattern &pat_send = it_send->second;
    const HaloPattern &pat_recv = it_recv->second;
    const int ncomp = desc->ncomp;

    //=========================================================================
    // 3. 遍历所有 Parallel HaloRegion：每个 region 对应一次 “本块 <-> 邻居 rank” 的通信，打包
    //-------------------------------------------------------------------------
    // 检测缓冲空间是否足够
    const int num_face_send = pat_send.regions.size();
    // Parallel corner patterns exclude material-coupling patches.  Therefore
    // both ends of every retained region have the same block physics and a
    // physics-scoped field is active (or inactive) on both ranks.  Skip the
    // inactive pair locally instead of dereferencing its empty FieldBlock.
    std::vector<unsigned char> send_active(num_face_send, 0);
    if (send_buf.size() < num_face_send)
        send_buf.resize(num_face_send);
    if (req_send.size() < num_face_send)
        req_send.resize(num_face_send);
    if (stat_send.size() < num_face_send)
        stat_send.resize(num_face_send);
    if (length.size() < num_face_send)
        length.resize(num_face_send);
    //-------------------------------------------------------------------------
    // 打包
    int index = 0;
    for (const HaloRegion &r : pat_send.regions)
    {
        FieldBlock &fb = fld_->field(fid, r.this_block); // 本 rank 上的块
        send_active[index] = fb.is_allocated() ? 1 : 0;
        if (!send_active[index])
        {
            length[index] = 0;
            ++index;
            continue;
        }

        const Box3 &sb = r.send_box; // 本块 inner strip

        const int32_t n_total = (sb.hi.i - sb.lo.i) *
                                (sb.hi.j - sb.lo.j) *
                                (sb.hi.k - sb.lo.k) *
                                ncomp;
        length[index] = n_total;
        // 4. 确保缓冲区足够大（复用 send_buf_ ）
        if (send_buf[index].size() < n_total)
            send_buf[index].resize(n_total);

        // 5. 打包：本块 inner strip -> send_buf_
        const TOPO::IndexTransform &transform = r.trans;

        // 1) 本地 send_box 的尺寸
        int loc_lo[3] = {sb.lo.i, sb.lo.j, sb.lo.k};
        int loc_hi[3] = {sb.hi.i - 1, sb.hi.j - 1, sb.hi.k - 1}; // 形成闭区间
        int len_loc[3] = {sb.hi.i - sb.lo.i, sb.hi.j - sb.lo.j, sb.hi.k - sb.lo.k};

        // 2) 获取目标块对应区域的范围，这里暂时用闭区间，主要获得最小值neighbor的lo.ijk, 存入tar_ref
        // tar_ref作为起点 用于编码
        int offset[3] = {r.trans.offset.i,
                         r.trans.offset.j,
                         r.trans.offset.k};
        int tar1[3], tar2[3], tar_ref[3];
        for (int d = 0; d < 3; ++d)
            tar1[transform.perm[d]] = transform.sign[d] * loc_lo[d] + offset[d];
        for (int d = 0; d < 3; ++d)
            tar2[transform.perm[d]] = transform.sign[d] * loc_hi[d] + offset[d];

        tar_ref[0] = (tar1[0] <= tar2[0]) ? tar1[0] : tar2[0];
        tar_ref[1] = (tar1[1] <= tar2[1]) ? tar1[1] : tar2[1];
        tar_ref[2] = (tar1[2] <= tar2[2]) ? tar1[2] : tar2[2];

        // 3) 对应到“邻居坐标系”下的尺寸：
        //    邻居 axis_nb = perm[d] 这一维的长度就是本地 d 维的长度
        int len_nb[3] = {0, 0, 0};
        for (int d = 0; d < 3; ++d)
            len_nb[transform.perm[d]] = len_loc[d];

        // 4) 真正打包：本地 (i,j,k) -> 邻居坐标 nb -> 相对坐标 (ri,rj,rk) -> buffer index
        int ijk[3], tar_ijk[3];
        for (ijk[0] = sb.lo.i; ijk[0] < sb.hi.i; ijk[0]++)
            for (ijk[1] = sb.lo.j; ijk[1] < sb.hi.j; ijk[1]++)
                for (ijk[2] = sb.lo.k; ijk[2] < sb.hi.k; ijk[2]++)
                {
                    for (int ii = 0; ii < 3; ii++)
                        tar_ijk[transform.perm[ii]] = transform.sign[ii] * ijk[ii] + offset[ii];
                    // 获得了目标tar_ijk[]坐标，现在设定编码方式：
                    //  按for i j k comp的顺序排列，以较小值tar_ref为起点
                    // 以 nb_ref 为原点的“邻居相对坐标”
                    int ri = tar_ijk[0] - tar_ref[0];
                    int rj = tar_ijk[1] - tar_ref[1];
                    int rk = tar_ijk[2] - tar_ref[2];
                    int32_t base = ((ri * len_nb[1] + rj) * len_nb[2] + rk) * ncomp;
                    // 真正的数据拷贝
                    for (int m = 0; m < ncomp; ++m)
                        send_buf[index][base + m] = fb(ijk[0], ijk[1], ijk[2], m);
                }

        // 5) 记录regions的个数
        index++;
    }

    //=========================================================================
    // 4. 作为接受块还需要单独开接受的空间
    //-------------------------------------------------------------------------
    // 检测缓冲空间是否足够
    const int num_face_recv = pat_recv.regions.size();
    std::vector<unsigned char> recv_active(num_face_recv, 0);
    if (recv_buf.size() < num_face_recv)
        recv_buf.resize(num_face_recv);
    if (req_recv.size() < num_face_recv)
        req_recv.resize(num_face_recv);
    if (stat_recv.size() < num_face_recv)
        stat_recv.resize(num_face_recv);
    if (length_corner_recv.size() < num_face_recv)
        length_corner_recv.resize(num_face_recv);
    //-------------------------------------------------------------------------
    // 统计空间
    int index_recv = 0;
    for (const HaloRegion &r : pat_recv.regions)
    {
        FieldBlock &fb = fld_->field(fid, r.this_block); // 本 rank 上的块
        recv_active[index_recv] = fb.is_allocated() ? 1 : 0;
        if (!recv_active[index_recv])
        {
            length_corner_recv[index_recv] = 0;
            ++index_recv;
            continue;
        }

        const Box3 &rb = r.recv_box; // 本块 ghost strip

        const int32_t n_total = (rb.hi.i - rb.lo.i) *
                                (rb.hi.j - rb.lo.j) *
                                (rb.hi.k - rb.lo.k) *
                                ncomp;
        length_corner_recv[index_recv] = n_total;
        // 确保缓冲区足够大（复用 recv_buf_ ）
        if (recv_buf[index_recv].size() < n_total)
            recv_buf[index_recv].resize(n_total);

        // 5) 记录regions的个数
        index_recv++;
    }

    //=========================================================================
    // 5. MPI
    //-------------------------------------------------------------------------
    // 等待
    // PARALLEL::mpi_barrier();
    //-------------------------------------------------------------------------
    // mpi发送接收
    index = 0;
    for (const HaloRegion &r : pat_send.regions)
    {
        if (send_active[index])
            PARALLEL::mpi_data_send(r.neighbor_rank, r.send_flag, send_buf[index].data(), length[index], &(req_send[index]));
        else
            req_send[index] = MPI_REQUEST_NULL;
        index++;
    }
    index_recv = 0;
    for (const HaloRegion &r : pat_recv.regions)
    {
        if (recv_active[index_recv])
            PARALLEL::mpi_data_recv(r.neighbor_rank, r.recv_flag, recv_buf[index_recv].data(), length_corner_recv[index_recv], &(req_recv[index_recv]));
        else
            req_recv[index_recv] = MPI_REQUEST_NULL;
        index_recv++;
    }
    //----------------------------------------------------------------------
    // 等待完成
    int num_face_comm_send = num_face_send;
    int num_face_comm_recv = num_face_recv;
    PARALLEL::mpi_wait(num_face_comm_send, req_send.data(), stat_send.data());
    PARALLEL::mpi_wait(num_face_comm_recv, req_recv.data(), stat_recv.data());
    //----------------------------------------------------------------------

    //=========================================================================
    // 6. 解包：recv_buf -> 本块的 recv_box（ghost cells）
    //----------------------------------------------------------------------
    index_recv = 0;
    for (const HaloRegion &r : pat_recv.regions)
    {
        if (!recv_active[index_recv])
        {
            ++index_recv;
            continue;
        }
        FieldBlock &fb = fld_->field(fid, r.this_block); // 本 rank 上的块
        const Box3 &rb = r.recv_box;                     // 要填充的 halo 区域

        const int ni = rb.hi.i - rb.lo.i;
        const int nj = rb.hi.j - rb.lo.j;
        const int nk = rb.hi.k - rb.lo.k;

        const std::vector<double> &buf = recv_buf[index_recv];

        // 这里我们假设发送端按照“邻居坐标系（也就是当前块坐标系）”
        // 以 comp , k 最快、再 j、再 i 的顺序排布：
        //
        // base = (((i_rel) * nj + j_rel) * nk + k_rel) * ncomp
        //
        // 其中 i_rel = i - rb.lo.i，j_rel = j - rb.lo.j，k_rel = k - rb.lo.k
        for (int ii = 0; ii < ni; ++ii)
        {
            int i_loc = rb.lo.i + ii;
            for (int jj = 0; jj < nj; ++jj)
            {
                int j_loc = rb.lo.j + jj;
                for (int kk = 0; kk < nk; ++kk)
                {
                    int k_loc = rb.lo.k + kk;

                    int32_t base = (((ii * nj + jj) * nk) + kk) * ncomp;

                    for (int m = 0; m < ncomp; ++m)
                    {
                        fb(i_loc, j_loc, k_loc, m) = buf[base + m];
                    }
                }
            }
        }

        ++index_recv;
    }
}

void Halo::exchange_parallel_edge(const std::string &field_name)
{
    //=========================================================================
    // 1. 找到这个 field 的 descriptor 和 fid
    const FieldDescriptor *desc = nullptr;
    int fid = -1;

    for (int id = 0; id < fld_->num_fields(); ++id)
    {
        if (field_name == fld_->descriptor(id).name)
        {
            desc = &(fld_->descriptor(id));
            fid = id;
            break;
        }
    }
    if (desc == nullptr)
    {
        std::cout << "Fatal Error!!! Can not find the field:\t" << field_name << std::endl;
        std::exit(-1);
    }

    PatternKey key{desc->location, desc->nghost};

    //=========================================================================
    // 2. 找到对应的 parallel pattern
    auto it_send = parallel_edge_patterns_send.find(key);
    if (it_send == parallel_edge_patterns_send.end())
    {
        std::cout << "Fatal Error!!! Can not find Parallel Halo pattern of field:\t"
                  << field_name << std::endl;
        std::exit(-1);
    }

    auto it_recv = parallel_edge_patterns_recv.find(key);
    if (it_recv == parallel_edge_patterns_recv.end())
    {
        std::cout << "Fatal Error!!! Can not find Parallel Halo pattern of field:\t"
                  << field_name << std::endl;
        std::exit(-1);
    }

    const HaloPattern &pat_send = it_send->second;
    const HaloPattern &pat_recv = it_recv->second;
    const int ncomp = desc->ncomp;

    //=========================================================================
    // 3. 遍历所有 Parallel HaloRegion：每个 region 对应一次 “本块 <-> 邻居 rank” 的通信，打包
    //-------------------------------------------------------------------------
    // 检测缓冲空间是否足够
    const int num_face_send = pat_send.regions.size();
    // See exchange_parallel_vertex(): non-coupling edge peers have matching
    // field allocation, so inactive regions must post no MPI request.
    std::vector<unsigned char> send_active(num_face_send, 0);
    if (send_buf.size() < num_face_send)
        send_buf.resize(num_face_send);
    if (req_send.size() < num_face_send)
        req_send.resize(num_face_send);
    if (stat_send.size() < num_face_send)
        stat_send.resize(num_face_send);
    if (length.size() < num_face_send)
        length.resize(num_face_send);
    //-------------------------------------------------------------------------
    // 打包
    int index = 0;
    for (const HaloRegion &r : pat_send.regions)
    {
        FieldBlock &fb = fld_->field(fid, r.this_block); // 本 rank 上的块
        send_active[index] = fb.is_allocated() ? 1 : 0;
        if (!send_active[index])
        {
            length[index] = 0;
            ++index;
            continue;
        }

        const Box3 &sb = r.send_box; // 本块 inner strip

        const int32_t n_total = (sb.hi.i - sb.lo.i) *
                                (sb.hi.j - sb.lo.j) *
                                (sb.hi.k - sb.lo.k) *
                                ncomp;
        length[index] = n_total;
        // 4. 确保缓冲区足够大（复用 send_buf_ ）
        if (send_buf[index].size() < n_total)
            send_buf[index].resize(n_total);

        // 5. 打包：本块 inner strip -> send_buf_
        const TOPO::IndexTransform &transform = r.trans;

        // 1) 本地 send_box 的尺寸
        int loc_lo[3] = {sb.lo.i, sb.lo.j, sb.lo.k};
        int loc_hi[3] = {sb.hi.i - 1, sb.hi.j - 1, sb.hi.k - 1}; // 形成闭区间
        int len_loc[3] = {sb.hi.i - sb.lo.i, sb.hi.j - sb.lo.j, sb.hi.k - sb.lo.k};

        // 2) 获取目标块对应区域的范围，这里暂时用闭区间，主要获得最小值neighbor的lo.ijk, 存入tar_ref
        // tar_ref作为起点 用于编码
        int offset[3] = {r.trans.offset.i,
                         r.trans.offset.j,
                         r.trans.offset.k};
        int tar1[3], tar2[3], tar_ref[3];
        for (int d = 0; d < 3; ++d)
            tar1[transform.perm[d]] = transform.sign[d] * loc_lo[d] + offset[d];
        for (int d = 0; d < 3; ++d)
            tar2[transform.perm[d]] = transform.sign[d] * loc_hi[d] + offset[d];

        tar_ref[0] = (tar1[0] <= tar2[0]) ? tar1[0] : tar2[0];
        tar_ref[1] = (tar1[1] <= tar2[1]) ? tar1[1] : tar2[1];
        tar_ref[2] = (tar1[2] <= tar2[2]) ? tar1[2] : tar2[2];

        // 3) 对应到“邻居坐标系”下的尺寸：
        //    邻居 axis_nb = perm[d] 这一维的长度就是本地 d 维的长度
        int len_nb[3] = {0, 0, 0};
        for (int d = 0; d < 3; ++d)
            len_nb[transform.perm[d]] = len_loc[d];

        // 4) 真正打包：本地 (i,j,k) -> 邻居坐标 nb -> 相对坐标 (ri,rj,rk) -> buffer index
        int ijk[3], tar_ijk[3];
        for (ijk[0] = sb.lo.i; ijk[0] < sb.hi.i; ijk[0]++)
            for (ijk[1] = sb.lo.j; ijk[1] < sb.hi.j; ijk[1]++)
                for (ijk[2] = sb.lo.k; ijk[2] < sb.hi.k; ijk[2]++)
                {
                    for (int ii = 0; ii < 3; ii++)
                        tar_ijk[transform.perm[ii]] = transform.sign[ii] * ijk[ii] + offset[ii];
                    // 获得了目标tar_ijk[]坐标，现在设定编码方式：
                    //  按for i j k comp的顺序排列，以较小值tar_ref为起点
                    // 以 nb_ref 为原点的“邻居相对坐标”
                    int ri = tar_ijk[0] - tar_ref[0];
                    int rj = tar_ijk[1] - tar_ref[1];
                    int rk = tar_ijk[2] - tar_ref[2];
                    int32_t base = ((ri * len_nb[1] + rj) * len_nb[2] + rk) * ncomp;
                    // 真正的数据拷贝
                    for (int m = 0; m < ncomp; ++m)
                        send_buf[index][base + m] = fb(ijk[0], ijk[1], ijk[2], m);
                }

        // 5) 记录regions的个数
        index++;
    }

    //=========================================================================
    // 4. 作为接受块还需要单独开接受的空间
    //-------------------------------------------------------------------------
    // 检测缓冲空间是否足够
    const int num_face_recv = pat_recv.regions.size();
    std::vector<unsigned char> recv_active(num_face_recv, 0);
    if (recv_buf.size() < num_face_recv)
        recv_buf.resize(num_face_recv);
    if (req_recv.size() < num_face_recv)
        req_recv.resize(num_face_recv);
    if (stat_recv.size() < num_face_recv)
        stat_recv.resize(num_face_recv);
    if (length_corner_recv.size() < num_face_recv)
        length_corner_recv.resize(num_face_recv);
    //-------------------------------------------------------------------------
    // 统计空间
    int index_recv = 0;
    for (const HaloRegion &r : pat_recv.regions)
    {
        FieldBlock &fb = fld_->field(fid, r.this_block); // 本 rank 上的块
        recv_active[index_recv] = fb.is_allocated() ? 1 : 0;
        if (!recv_active[index_recv])
        {
            length_corner_recv[index_recv] = 0;
            ++index_recv;
            continue;
        }

        const Box3 &rb = r.recv_box; // 本块 ghost strip

        const int32_t n_total = (rb.hi.i - rb.lo.i) *
                                (rb.hi.j - rb.lo.j) *
                                (rb.hi.k - rb.lo.k) *
                                ncomp;
        length_corner_recv[index_recv] = n_total;
        // 确保缓冲区足够大（复用 recv_buf_ ）
        if (recv_buf[index_recv].size() < n_total)
            recv_buf[index_recv].resize(n_total);

        // 5) 记录regions的个数
        index_recv++;
    }

    //=========================================================================
    // 5. MPI
    //-------------------------------------------------------------------------
    // 等待
    // PARALLEL::mpi_barrier();
    //-------------------------------------------------------------------------
    // mpi发送接收
    index = 0;
    for (const HaloRegion &r : pat_send.regions)
    {
        if (send_active[index])
            PARALLEL::mpi_data_send(r.neighbor_rank, r.send_flag, send_buf[index].data(), length[index], &(req_send[index]));
        else
            req_send[index] = MPI_REQUEST_NULL;
        index++;
    }
    index_recv = 0;
    for (const HaloRegion &r : pat_recv.regions)
    {
        if (recv_active[index_recv])
            PARALLEL::mpi_data_recv(r.neighbor_rank, r.recv_flag, recv_buf[index_recv].data(), length_corner_recv[index_recv], &(req_recv[index_recv]));
        else
            req_recv[index_recv] = MPI_REQUEST_NULL;
        index_recv++;
    }
    //----------------------------------------------------------------------
    // 等待完成
    int num_face_comm_send = num_face_send;
    int num_face_comm_recv = num_face_recv;
    PARALLEL::mpi_wait(num_face_comm_send, req_send.data(), stat_send.data());
    PARALLEL::mpi_wait(num_face_comm_recv, req_recv.data(), stat_recv.data());
    //----------------------------------------------------------------------

    //=========================================================================
    // 6. 解包：recv_buf -> 本块的 recv_box（ghost cells）
    //----------------------------------------------------------------------
    index_recv = 0;
    for (const HaloRegion &r : pat_recv.regions)
    {
        if (!recv_active[index_recv])
        {
            ++index_recv;
            continue;
        }
        FieldBlock &fb = fld_->field(fid, r.this_block); // 本 rank 上的块
        const Box3 &rb = r.recv_box;                     // 要填充的 halo 区域

        const int ni = rb.hi.i - rb.lo.i;
        const int nj = rb.hi.j - rb.lo.j;
        const int nk = rb.hi.k - rb.lo.k;

        const std::vector<double> &buf = recv_buf[index_recv];

        // 这里我们假设发送端按照“邻居坐标系（也就是当前块坐标系）”
        // 以 comp , k 最快、再 j、再 i 的顺序排布：
        //
        // base = (((i_rel) * nj + j_rel) * nk + k_rel) * ncomp
        //
        // 其中 i_rel = i - rb.lo.i，j_rel = j - rb.lo.j，k_rel = k - rb.lo.k
        for (int ii = 0; ii < ni; ++ii)
        {
            int i_loc = rb.lo.i + ii;
            for (int jj = 0; jj < nj; ++jj)
            {
                int j_loc = rb.lo.j + jj;
                for (int kk = 0; kk < nk; ++kk)
                {
                    int k_loc = rb.lo.k + kk;

                    int32_t base = (((ii * nj + jj) * nk) + kk) * ncomp;

                    for (int m = 0; m < ncomp; ++m)
                    {
                        fb(i_loc, j_loc, k_loc, m) = buf[base + m];
                    }
                }
            }
        }

        ++index_recv;
    }
}
