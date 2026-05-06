#include "0_basic/MPI_WRAPPER.h"

namespace PARALLEL
{
    void mpi_initial(int arg, char **argv)
    {
        MPI_Init(&arg, &argv);
    }
    void mpi_finalize()
    {
        MPI_Finalize();
    }
    void mpi_barrier()
    {
        MPI_Barrier(MPI_COMM_WORLD);
    }
    void mpi_max(double *number, double *maximun, int num)
    {
        MPI_Allreduce(number, maximun, num, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    }
    void mpi_min(double *number, double *minimun, int num)
    {
        MPI_Allreduce(number, minimun, num, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    }
    void mpi_data_send(int tar_id, int send_flag, double *data, int length, MPI_Request *status)
    {
        MPI_Isend(data, length, MPI_DOUBLE, tar_id, send_flag, MPI_COMM_WORLD, status);
    }
    void mpi_data_recv(int tar_id, int recv_flag, double *data, int length, MPI_Request *status)
    {
        MPI_Irecv(data, length, MPI_DOUBLE, tar_id, recv_flag, MPI_COMM_WORLD, status);
    }
    void mpi_rank(int *rank)
    {
        MPI_Comm_rank(MPI_COMM_WORLD, rank);
    }
    void mpi_size(int *size)
    {
        MPI_Comm_size(MPI_COMM_WORLD, size);
    }
    void mpi_wait(int &count, MPI_Request *request_s, MPI_Status *status_s)
    {
        MPI_Waitall(count, request_s, status_s);
    }
    void mpi_allgatherv(double *data_send, int send_count, double *data_recv, int *recv_count, int *displs)
    {
        MPI_Allgatherv(data_send, send_count, MPI_DOUBLE, data_recv, recv_count, displs, MPI_DOUBLE, MPI_COMM_WORLD);
    }
    void mpi_allgatherv(int *data_send, int send_count, int *data_recv, int *recv_count, int *displs)
    {
        MPI_Allgatherv(data_send, send_count, MPI_INT, data_recv, recv_count, displs, MPI_INT, MPI_COMM_WORLD);
    }
    void mpi_gather(int *data_send, int send_count, int *data_recv, int recv_count, int root)
    {
        MPI_Gather(data_send, send_count, MPI_INT, data_recv, recv_count, MPI_INT, root, MPI_COMM_WORLD);
    }
    void mpi_bcast(int *data_send, int send_count, int root)
    {
        MPI_Bcast(data_send, send_count, MPI_INT, root, MPI_COMM_WORLD);
    }
    void mpi_sum(double *number, double *sum, int num)
    {
        MPI_Allreduce(number, sum, num, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }

    void mpi_alltoall(int *send, int send_num, int *recv, int recv_num)
    {
        MPI_Alltoall(send, send_num, MPI_INT, recv, recv_num, MPI_INT, MPI_COMM_WORLD);
    }
}