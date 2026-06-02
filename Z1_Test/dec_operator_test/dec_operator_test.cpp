#include "Z1_Test_Common.h"

#include "8_dec/DecDiagnostics.h"
#include "8_dec/DecOps.h"

#include <exception>
#include <sstream>

int main(int argc, char **argv)
{
    PARALLEL::mpi_initial(argc, argv);
    int exit_code = 0;
    try
    {
        Z1_TEST::CaseContext ctx = Z1_TEST::load_case_and_topology();
        const int nghost = ctx.param.GetInt("ngg");
        Field field(&ctx.grid, &ctx.param, nghost);
        Z1_TEST::register_dec_test_fields(field, nghost);
        Z1_TEST::fill_edge_test_form(field);

        const DEC::EdgeFormNames e_names{"E_xi", "E_eta", "E_zeta"};
        const DEC::FaceFormNames b_names{"B_xi", "B_eta", "B_zeta"};
        DEC::d1_edge_to_face(field, e_names, b_names);
        DEC::d2_face_to_cell(field, b_names, "divB");
        const double max_div_curl = DEC::max_abs_div_curl(field, e_names);

        if (ctx.myid == 0)
            std::cout << "DEC operator tests\n";

        std::ostringstream detail;
        detail << "max=" << max_div_curl;
        bool passed = true;
        passed &= Z1_TEST::print_pass_fail(ctx.myid,
                                           "9. exterior_derivative_1_to_2 then div(curl(E)) == 0",
                                           max_div_curl < 1.0e-10,
                                           detail.str());
        passed &= Z1_TEST::print_result(ctx.myid, "SKIP",
                                        "10. codifferential_2_to_1 adjointness",
                                        "(framework placeholder: codifferential_2_to_1 API not present yet)");
        exit_code = passed ? 0 : 1;
    }
    catch (const std::exception &ex)
    {
        int myid = 0;
        PARALLEL::mpi_rank(&myid);
        if (myid == 0)
            std::cerr << "dec_operator_test failed with exception: " << ex.what() << "\n";
        exit_code = 1;
    }
    PARALLEL::mpi_finalize();
    return exit_code;
}
