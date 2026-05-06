#pragma once
#include <mpi.h>

namespace PARALLEL
{
    // 封装的mpi函数，不同的MPI可以只需更改这里封装的函数即可
    void mpi_initial(int argc, char **argv);
    void mpi_finalize();
    void mpi_rank(int *rank);
    void mpi_size(int *size);
    void mpi_barrier();
    void mpi_wait(int &count, MPI_Request *request_s, MPI_Status *status_s);

    void mpi_max(double *number, double *maximun, int num);
    void mpi_min(double *number, double *minimun, int num);
    void mpi_allgatherv(double *data_send, int send_count, double *data_recv, int *recv_count, int *displs); // DOUBLE
    void mpi_allgatherv(int *data_send, int send_count, int *data_recv, int *recv_count, int *displs);       // INT
    void mpi_gather(int *data_send, int send_count, int *data_recv, int recv_count, int root);
    void mpi_bcast(int *data_send, int send_count, int root);
    void mpi_sum(double *number, double *sum, int num);

    void mpi_data_send(int tar_id, int send_flag, double *data, int length, MPI_Request *status);
    void mpi_data_recv(int tar_id, int recv_flag, double *data, int length, MPI_Request *status);

    void mpi_alltoall(int *send, int send_num, int *recv, int recv_num);
}