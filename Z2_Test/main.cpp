#include "Z0_Runner.h"

#include "0_basic/MPI_WRAPPER.h"

int main(int argc, char **argv)
{
    PARALLEL::mpi_initial(argc, argv);
    const int rc = Z0::run(argc, argv);
    PARALLEL::mpi_finalize();
    return rc;
}
