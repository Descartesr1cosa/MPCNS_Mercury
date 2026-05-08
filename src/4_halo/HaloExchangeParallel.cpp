#include "4_halo/Halo.h"
#include "0_basic/MPI_WRAPPER.h"

void Halo::exchange_parallel(std::string field_name)
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
    //-------------------------------------------------------------------------
    // 打包
    int index = 0;
    for (const HaloRegion &r : pat.regions)
    {
        FieldBlock &fb = fld_->field(fid, r.this_block); // 本 rank 上的块

        const Box3 &sb = r.send_box; // 本块 inner strip

#if if_Debug_Field_Array == 1
        // send/recv 区域尺寸（应该相同）
        const int ni = sb.hi.i - sb.lo.i;
        const int nj = sb.hi.j - sb.lo.j;
        const int nk = sb.hi.k - sb.lo.k;

        const Box3 &rb = r.recv_box; // recv
        const int ni_r = rb.hi.i - rb.lo.i;
        const int nj_r = rb.hi.j - rb.lo.j;
        const int nk_r = rb.hi.k - rb.lo.k;
        // 防御式：确保 send_box 和 recv_box 的大小一致
        if (ni != ni_r || nj != nj_r || nk != nk_r)
        {
            std::cout << "Fatal Error!!! Parallel Halo send/recv box size mismatch "
                      << "(field=" << field_name << ", block=" << r.this_block << ")\n";
            std::exit(-1);
        }
#endif
        const int32_t n_total = (sb.hi.i - sb.lo.i) *
                                (sb.hi.j - sb.lo.j) *
                                (sb.hi.k - sb.lo.k) *
                                ncomp;
        length[index] = n_total;
        // 4. 确保缓冲区足够大（复用 send_buf_ / recv_buf_）
        if (send_buf[index].size() < n_total)
            send_buf[index].resize(n_total);
        if (recv_buf[index].size() < n_total)
            recv_buf[index].resize(n_total);

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
    // mpi发送接收
    index = 0;
    for (const HaloRegion &r : pat.regions)
    {
        PARALLEL::mpi_data_send(r.neighbor_rank, r.send_flag, send_buf[index].data(), length[index], &(req_send[index]));
        PARALLEL::mpi_data_recv(r.neighbor_rank, r.recv_flag, recv_buf[index].data(), length[index], &(req_recv[index]));
        index++;
    }
    //----------------------------------------------------------------------
    // 等待完成
    int num_face_comm = num_face;
    PARALLEL::mpi_wait(num_face_comm, req_send.data(), stat_send.data());
    PARALLEL::mpi_wait(num_face_comm, req_recv.data(), stat_recv.data());
    PARALLEL::mpi_barrier();
    //----------------------------------------------------------------------

    //----------------------------------------------------------------------
    // 6. 解包：recv_buf -> 本块的 recv_box（ghost cells）
    //----------------------------------------------------------------------
    index = 0;
    for (const HaloRegion &r : pat.regions)
    {
        FieldBlock &fb = fld_->field(fid, r.this_block); // 本 rank 上的块
        const Box3 &rb = r.recv_box;                     // 要填充的 halo 区域

        const int ni = rb.hi.i - rb.lo.i;
        const int nj = rb.hi.j - rb.lo.j;
        const int nk = rb.hi.k - rb.lo.k;

#if if_Debug_Field_Array == 1
        // 防御式检查：与打包时记录的长度一致

        const int32_t n_total = ni * nj * nk * ncomp;
        if (n_total != length[index])
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

void Halo::exchange_parallel_vertex(std::string field_name)
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
        PARALLEL::mpi_data_send(r.neighbor_rank, r.send_flag, send_buf[index].data(), length[index], &(req_send[index]));
        index++;
    }
    index_recv = 0;
    for (const HaloRegion &r : pat_recv.regions)
    {
        PARALLEL::mpi_data_recv(r.neighbor_rank, r.recv_flag, recv_buf[index_recv].data(), length_corner_recv[index_recv], &(req_recv[index_recv]));
        index_recv++;
    }
    //----------------------------------------------------------------------
    // 等待完成
    int num_face_comm_send = num_face_send;
    int num_face_comm_recv = num_face_recv;
    PARALLEL::mpi_wait(num_face_comm_send, req_send.data(), stat_send.data());
    PARALLEL::mpi_wait(num_face_comm_recv, req_recv.data(), stat_recv.data());
    PARALLEL::mpi_barrier();
    //----------------------------------------------------------------------

    //=========================================================================
    // 6. 解包：recv_buf -> 本块的 recv_box（ghost cells）
    //----------------------------------------------------------------------
    index_recv = 0;
    for (const HaloRegion &r : pat_recv.regions)
    {
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

void Halo::exchange_parallel_edge(std::string field_name)
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
        PARALLEL::mpi_data_send(r.neighbor_rank, r.send_flag, send_buf[index].data(), length[index], &(req_send[index]));
        index++;
    }
    index_recv = 0;
    for (const HaloRegion &r : pat_recv.regions)
    {
        PARALLEL::mpi_data_recv(r.neighbor_rank, r.recv_flag, recv_buf[index_recv].data(), length_corner_recv[index_recv], &(req_recv[index_recv]));
        index_recv++;
    }
    //----------------------------------------------------------------------
    // 等待完成
    int num_face_comm_send = num_face_send;
    int num_face_comm_recv = num_face_recv;
    PARALLEL::mpi_wait(num_face_comm_send, req_send.data(), stat_send.data());
    PARALLEL::mpi_wait(num_face_comm_recv, req_recv.data(), stat_recv.data());
    PARALLEL::mpi_barrier();
    //----------------------------------------------------------------------

    //=========================================================================
    // 6. 解包：recv_buf -> 本块的 recv_box（ghost cells）
    //----------------------------------------------------------------------
    index_recv = 0;
    for (const HaloRegion &r : pat_recv.regions)
    {
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
