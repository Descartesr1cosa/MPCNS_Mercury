#include "Z1_Test_Common.h"

#include "7_metric/Metric.h"

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
        const METRIC::MetricDiagnostics diag = METRIC::diagnose_metric_fields(field);

        if (ctx.myid == 0)
            std::cout << "Metric/Hodge tests\n";

        const bool hodge_finite =
            diag.Hodge_star_inverse_2form_to_1form_nonfinite == 0 &&
            diag.Hodge_star_2form_to_1form_nonfinite == 0;

        const bool hodge_positive =
            Z1_TEST::finite_positive_field(field, "Hodge_star_2form_to_1form_face_xi_lumped", true) &&
            Z1_TEST::finite_positive_field(field, "Hodge_star_2form_to_1form_face_eta_lumped", true) &&
            Z1_TEST::finite_positive_field(field, "Hodge_star_2form_to_1form_face_zeta_lumped", true);

        const bool axis_finite =
            diag.jac_nonpositive >= 0 &&
            diag.area_nonpositive >= 0 &&
            diag.dl_nonpositive >= 0 &&
            hodge_finite;

        std::ostringstream finite_detail;
        finite_detail << "hodge2to1_nonfinite=" << diag.Hodge_star_2form_to_1form_nonfinite
                      << " inv2to1_nonfinite=" << diag.Hodge_star_inverse_2form_to_1form_nonfinite;

        std::ostringstream positive_detail;
        positive_detail << "near_axis_singular=" << diag.near_axis_singular
                        << " near_axis_capped=" << diag.near_axis_capped;

        bool passed = true;
        passed &= Z1_TEST::print_pass_fail(ctx.myid, "11. Hodge scale finite",
                                           hodge_finite, finite_detail.str());
        passed &= Z1_TEST::print_pass_fail(ctx.myid,
                                           "12. Hodge positive except masked singular DOFs",
                                           hodge_positive, positive_detail.str());
        passed &= Z1_TEST::print_pass_fail(ctx.myid,
                                           "13. axis/collapsed entity does not produce NaN/Inf",
                                           axis_finite, finite_detail.str());

        exit_code = passed ? 0 : 1;
    }
    catch (const std::exception &ex)
    {
        int myid = 0;
        PARALLEL::mpi_rank(&myid);
        if (myid == 0)
            std::cerr << "metric_hodge_test failed with exception: " << ex.what() << "\n";
        exit_code = 1;
    }
    PARALLEL::mpi_finalize();
    return exit_code;
}
