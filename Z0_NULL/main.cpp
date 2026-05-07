#include "Z0_NullRunner.h"

#include "0_basic/MPI_WRAPPER.h"

int main(int argc, char **argv)
{
    PARALLEL::mpi_initial(argc, argv);

    const int rc = Z0_NULL::run(argc, argv);

    PARALLEL::mpi_finalize();
    return rc;
}
