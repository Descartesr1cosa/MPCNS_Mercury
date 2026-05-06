#include "4_halo/1_MPCNS_Halo.h"
#include "0_basic/MPI_WRAPPER.h"

void Halo::mpi_exchange_edge_meta(
    const std::map<int, std::vector<EdgeMeta>> &meta_to_send,
    std::vector<EdgeMeta> &recv_metas)
{
    int nrank;
    PARALLEL::mpi_size(&nrank);
    int myrank;
    PARALLEL::mpi_rank(&myrank);

    // 1) 每个 rank 要发给别人的数量数组 send_counts[nrank] / recv_counts[nrank]
    std::vector<int> send_counts(nrank, 0), recv_counts(nrank, 0);
    for (auto &kv : meta_to_send)
    {
        int dest = kv.first;
        send_counts[dest] = (int)kv.second.size();
    }

    // 用 Alltoall 或者一圈 Send/Recv 换得 recv_counts
    PARALLEL::mpi_alltoall(send_counts.data(), 1, recv_counts.data(), 1);

    // 2) 计算偏移、总数
    std::vector<int> sdispls(nrank, 0), rdispls(nrank, 0);
    int total_send = 0, total_recv = 0;
    for (int r = 0; r < nrank; ++r)
    {
        sdispls[r] = total_send;
        rdispls[r] = total_recv;
        total_send += send_counts[r];
        total_recv += recv_counts[r];
    }

    std::vector<EdgeMeta> send_buf(total_send);
    recv_metas.resize(total_recv);

    // 3) 把 meta_to_send 展平到 send_buf（按 rank 顺序）
    for (int r = 0; r < nrank; ++r)
    {
        auto it = meta_to_send.find(r);
        if (it == meta_to_send.end())
            continue;
        const auto &v = it->second;
        std::copy(v.begin(), v.end(), send_buf.begin() + sdispls[r]);
    }

    // 4) Alltoallv 交换元数据（用 MPI_BYTE 发送 struct）
    // 把“struct 个数”转换成“字节数”，调用 MPI_Alltoallv
    const int sz = static_cast<int>(sizeof(EdgeMeta));

    std::vector<int> send_counts_bytes(nrank), recv_counts_bytes(nrank);
    std::vector<int> sdispls_bytes(nrank), rdispls_bytes(nrank);
    for (int r = 0; r < nrank; ++r)
    {
        send_counts_bytes[r] = send_counts[r] * sz;
        recv_counts_bytes[r] = recv_counts[r] * sz;
        sdispls_bytes[r] = sdispls[r] * sz;
        rdispls_bytes[r] = rdispls[r] * sz;
    }

    // 5) 真正交换 meta 的 Alltoallv
    MPI_Alltoallv(
        // send
        reinterpret_cast<const void *>(send_buf.data()),
        send_counts_bytes.data(),
        sdispls_bytes.data(),
        MPI_BYTE,
        // recv
        reinterpret_cast<void *>(recv_metas.data()),
        recv_counts_bytes.data(),
        rdispls_bytes.data(),
        MPI_BYTE,
        // comm
        MPI_COMM_WORLD);
}

void Halo::mpi_exchange_vertex_meta(
    const std::map<int, std::vector<VertexMeta>> &meta_to_send,
    std::vector<VertexMeta> &recv_metas)
{
    int nrank;
    PARALLEL::mpi_size(&nrank);
    int myrank;
    PARALLEL::mpi_rank(&myrank);

    // 1) 每个 rank 要发给别人的数量数组 send_counts[nrank] / recv_counts[nrank]
    std::vector<int> send_counts(nrank, 0), recv_counts(nrank, 0);
    for (auto &kv : meta_to_send)
    {
        int dest = kv.first;
        send_counts[dest] = (int)kv.second.size();
    }

    // 用 Alltoall 或者一圈 Send/Recv 换得 recv_counts
    PARALLEL::mpi_alltoall(send_counts.data(), 1, recv_counts.data(), 1);

    // 2) 计算偏移、总数
    std::vector<int> sdispls(nrank, 0), rdispls(nrank, 0);
    int total_send = 0, total_recv = 0;
    for (int r = 0; r < nrank; ++r)
    {
        sdispls[r] = total_send;
        rdispls[r] = total_recv;
        total_send += send_counts[r];
        total_recv += recv_counts[r];
    }

    std::vector<VertexMeta> send_buf(total_send);
    recv_metas.resize(total_recv);

    // 3) 把 meta_to_send 展平到 send_buf（按 rank 顺序）
    for (int r = 0; r < nrank; ++r)
    {
        auto it = meta_to_send.find(r);
        if (it == meta_to_send.end())
            continue;
        const auto &v = it->second;
        std::copy(v.begin(), v.end(), send_buf.begin() + sdispls[r]);
    }

    // 4) Alltoallv 交换元数据（用 MPI_BYTE 发送 struct）
    // 把“struct 个数”转换成“字节数”，调用 MPI_Alltoallv
    const int sz = static_cast<int>(sizeof(VertexMeta));

    std::vector<int> send_counts_bytes(nrank), recv_counts_bytes(nrank);
    std::vector<int> sdispls_bytes(nrank), rdispls_bytes(nrank);
    for (int r = 0; r < nrank; ++r)
    {
        send_counts_bytes[r] = send_counts[r] * sz;
        recv_counts_bytes[r] = recv_counts[r] * sz;
        sdispls_bytes[r] = sdispls[r] * sz;
        rdispls_bytes[r] = rdispls[r] * sz;
    }

    // 5) 真正交换 meta 的 Alltoallv
    MPI_Alltoallv(
        // send
        reinterpret_cast<const void *>(send_buf.data()),
        send_counts_bytes.data(),
        sdispls_bytes.data(),
        MPI_BYTE,
        // recv
        reinterpret_cast<void *>(recv_metas.data()),
        recv_counts_bytes.data(),
        rdispls_bytes.data(),
        MPI_BYTE,
        // comm
        MPI_COMM_WORLD);
}
