#include "Z1_Test_Common.h"

#include "8_dec/DecOps.h"

#include <exception>
#include <sstream>

namespace
{
    bool check_edge_owner_alias_signs(const TOPO::Topology &topology)
    {
        for (const TOPO::EquivClass &cls : topology.edge_classes)
        {
            const int owner_sign = cls.owner.orient_sign;
            for (const TOPO::EquivMember &member : cls.members)
            {
                if (topology.owner_of(member.entity) != cls.owner.entity)
                    return false;
                if (topology.sign_to_owner(member.entity) != member.orient_sign * owner_sign)
                    return false;
            }
        }
        return true;
    }

    bool check_face_owner_alias_signs(const TOPO::Topology &topology)
    {
        for (const TOPO::EquivClass &cls : topology.face_classes)
        {
            const int owner_sign = cls.owner.orient_sign;
            for (const TOPO::EquivMember &member : cls.members)
            {
                if (topology.owner_of(member.entity) != cls.owner.entity)
                    return false;
                if (topology.sign_to_owner(member.entity) != member.orient_sign * owner_sign)
                    return false;
            }
        }
        return true;
    }
}

int main(int argc, char **argv)
{
    PARALLEL::mpi_initial(argc, argv);
    int exit_code = 0;
    try
    {
        Z1_TEST::CaseContext ctx = Z1_TEST::load_case_and_topology();
        if (ctx.myid == 0)
            std::cout << "Halo/Form tests\n";

        bool passed = true;
        passed &= Z1_TEST::print_pass_fail(ctx.myid, "5. edge 1-form owner/alias sign correct",
                                           check_edge_owner_alias_signs(ctx.topology));
        passed &= Z1_TEST::print_pass_fail(ctx.myid, "6. face 2-form owner/alias sign correct",
                                           check_face_owner_alias_signs(ctx.topology));

        passed &= Z1_TEST::print_result(ctx.myid, "SKIP", "7. block interface d1(E) consistent",
                                        "(framework placeholder: needs interface sampling oracle)");
        passed &= Z1_TEST::print_result(ctx.myid, "SKIP", "8. block interface d2(B) consistent",
                                        "(framework placeholder: needs interface sampling oracle)");

        exit_code = passed ? 0 : 1;
    }
    catch (const std::exception &ex)
    {
        int myid = 0;
        PARALLEL::mpi_rank(&myid);
        if (myid == 0)
            std::cerr << "halo_form_test failed with exception: " << ex.what() << "\n";
        exit_code = 1;
    }
    PARALLEL::mpi_finalize();
    return exit_code;
}
