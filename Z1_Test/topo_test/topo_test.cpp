#include "0_basic/1_MPCNS_Parameter.h"
#include "0_basic/MPI_WRAPPER.h"
#include "1_grid/1_MPCNS_Grid.h"
#include "2_topology/GlobalIncidence.h"
#include "2_topology/LocalIncidence.h"
#include "2_topology/TopologyBuilder.h"

#include <cstdint>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace
{
    struct TestSummary
    {
        long long checked = 0;
        long long failed = 0;
    };

    void mark_failure(TestSummary &summary)
    {
        ++summary.failed;
    }

    TestSummary check_local_boundary_boundary_face(const TOPO::Topology &topology)
    {
        TestSummary result;
        for (const auto &entry : topology.face2key)
        {
            ++result.checked;
            if (!TOPO::check_boundary_boundary_face_zero(entry.first))
                ++result.failed;
        }
        return result;
    }

    TestSummary check_local_boundary_boundary_cell(const TOPO::Topology &topology)
    {
        TestSummary result;
        for (const auto &entry : topology.cell_to_id)
        {
            ++result.checked;
            if (!TOPO::check_boundary_boundary_cell_zero(entry.first))
                ++result.failed;
        }
        return result;
    }

    bool valid_sign(int value)
    {
        return value == +1 || value == -1;
    }

    template <class Key, class KeyHash>
    bool contains_entity(const std::unordered_map<Key, std::vector<TOPO::EntityKey>, KeyHash> &members_by_key,
                         const Key &key,
                         const TOPO::EntityKey &entity)
    {
        const auto members_it = members_by_key.find(key);
        if (members_it == members_by_key.end())
            return false;
        for (const TOPO::EntityKey &member : members_it->second)
            if (member == entity)
                return true;
        return false;
    }

    template <class Key, class KeyHash>
    TestSummary check_owner_alias_maps(
        const std::unordered_map<TOPO::EntityKey, Key, TOPO::EntityKey::Hash> &entity_to_key,
        const std::unordered_map<TOPO::EntityKey, int8_t, TOPO::EntityKey::Hash> &entity_to_sign,
        const std::unordered_map<Key, std::vector<TOPO::EntityKey>, KeyHash> &members_by_key,
        const std::unordered_map<Key, TOPO::EntityKey, KeyHash> &owner_by_key,
        const std::unordered_map<TOPO::EntityKey, bool, TOPO::EntityKey::Hash> &is_owner_by_entity,
        const TOPO::Topology &topology)
    {
        TestSummary result;

        for (const auto &entry : entity_to_key)
        {
            const TOPO::EntityKey &entity = entry.first;
            const Key &key = entry.second;
            ++result.checked;

            const auto sign_it = entity_to_sign.find(entity);
            if (sign_it == entity_to_sign.end() || !valid_sign(sign_it->second))
            {
                mark_failure(result);
                continue;
            }

            const auto owner_it = owner_by_key.find(key);
            if (owner_it == owner_by_key.end())
                continue;

            if (!contains_entity(members_by_key, key, entity))
                mark_failure(result);

            const TOPO::EntityKey &owner = owner_it->second;
            const auto owner_sign_it = entity_to_sign.find(owner);
            if (owner_sign_it == entity_to_sign.end() || !valid_sign(owner_sign_it->second))
            {
                mark_failure(result);
                continue;
            }

            if (topology.owner_of(entity) != owner)
                mark_failure(result);

            const int expected_sign_to_owner =
                static_cast<int>(sign_it->second) * static_cast<int>(owner_sign_it->second);
            if (topology.sign_to_owner(entity) != expected_sign_to_owner)
                mark_failure(result);

            const auto is_owner_it = is_owner_by_entity.find(entity);
            if (is_owner_it == is_owner_by_entity.end() ||
                is_owner_it->second != (entity == owner) ||
                topology.is_owner(entity) != (entity == owner))
                mark_failure(result);
        }

        for (const auto &entry : members_by_key)
        {
            const Key &key = entry.first;
            const std::vector<TOPO::EntityKey> &members = entry.second;
            std::set<TOPO::EntityKey> unique_members;
            int owner_count = 0;

            const auto owner_it = owner_by_key.find(key);
            for (const TOPO::EntityKey &member : members)
            {
                ++result.checked;
                if (!unique_members.insert(member).second)
                    mark_failure(result);

                const auto key_it = entity_to_key.find(member);
                if (key_it == entity_to_key.end() || !(key_it->second == key))
                    mark_failure(result);

                if (owner_it != owner_by_key.end() && member == owner_it->second)
                    ++owner_count;
            }

            if (owner_it != owner_by_key.end() && owner_count != 1)
                mark_failure(result);
        }

        for (const auto &entry : owner_by_key)
        {
            const Key &key = entry.first;
            const TOPO::EntityKey &owner = entry.second;
            ++result.checked;

            const auto owner_key_it = entity_to_key.find(owner);
            if (owner_key_it == entity_to_key.end() || !(owner_key_it->second == key))
                mark_failure(result);

            const auto is_owner_it = is_owner_by_entity.find(owner);
            if (is_owner_it == is_owner_by_entity.end() || !is_owner_it->second)
                mark_failure(result);

            int owner_flag_count = 0;
            const auto members_it = members_by_key.find(key);
            if (members_it == members_by_key.end())
            {
                mark_failure(result);
                continue;
            }

            for (const TOPO::EntityKey &member : members_it->second)
            {
                const auto is_owner_it_for_member = is_owner_by_entity.find(member);
                if (is_owner_it_for_member != is_owner_by_entity.end() &&
                    is_owner_it_for_member->second)
                    ++owner_flag_count;
            }
            if (owner_flag_count != 1)
                mark_failure(result);
        }

        return result;
    }

    TestSummary check_owner_alias_uniqueness_and_signs(const TOPO::Topology &topology)
    {
        TestSummary result;

        for (const auto &entry : topology.node2eq)
        {
            ++result.checked;
            if (topology.owner_of(entry.first) != entry.second)
                mark_failure(result);
            if (topology.sign_to_owner(entry.first) != +1)
                mark_failure(result);
        }

        TestSummary edge_result = check_owner_alias_maps<TOPO::EdgeKey, TOPO::EdgeKey::Hash>(
            topology.edge2key,
            topology.edge2sign,
            topology.edge_members,
            topology.edge_owner,
            topology.edge_is_owner,
            topology);

        TestSummary face_result = check_owner_alias_maps<TOPO::FaceKey, TOPO::FaceKey::Hash>(
            topology.face2key,
            topology.face2sign,
            topology.face_members,
            topology.face_owner,
            topology.face_is_owner,
            topology);

        result.checked += edge_result.checked + face_result.checked;
        result.failed += edge_result.failed + face_result.failed;
        return result;
    }

    bool print_test_result(const std::string &name, bool passed, int myid)
    {
        if (myid == 0)
            std::cout << "  [" << (passed ? "PASS" : "FAIL") << "] " << name << "\n";
        return passed;
    }

    bool print_counted_test_result(const std::string &name, const TestSummary &summary, int myid)
    {
        const bool passed = summary.failed == 0;
        if (myid == 0)
        {
            std::cout << "  [" << (passed ? "PASS" : "FAIL") << "] " << name
                      << " checked=" << summary.checked
                      << " failed=" << summary.failed << "\n";
        }
        return passed;
    }
}

int main(int argc, char **argv)
{
    PARALLEL::mpi_initial(argc, argv);

    int exit_code = 0;
    try
    {
        int myid = 0;
        PARALLEL::mpi_rank(&myid);

        Param param;
        param.ReadParam(myid);

        Grid grid;
        grid.Grid_Preprocess(&param);

        const int dimension = param.GetInt("dimension");
        TOPO::Topology topology = TOPO::build_topology(grid, myid, dimension);
        TOPO::GlobalIncidence incidence(topology);

        if (myid == 0)
            std::cout << "Topology tests\n";

        bool passed = true;
        passed &= print_counted_test_result(
            "local boundary(boundary(face)) == 0",
            check_local_boundary_boundary_face(topology),
            myid);
        passed &= print_counted_test_result(
            "local boundary(boundary(cell)) == 0",
            check_local_boundary_boundary_cell(topology),
            myid);
        passed &= print_test_result("global D1 * D0 == 0",
                                    incidence.check_d1_d0_zero(),
                                    myid);
        passed &= print_test_result("global D2 * D1 == 0",
                                    incidence.check_d2_d1_zero(),
                                    myid);
        passed &= print_counted_test_result(
            "owner/alias uniqueness and sign consistency",
            check_owner_alias_uniqueness_and_signs(topology),
            myid);

        exit_code = passed ? 0 : 1;
    }
    catch (const std::exception &ex)
    {
        int myid = 0;
        PARALLEL::mpi_rank(&myid);
        if (myid == 0)
            std::cerr << "topo_test failed with exception: " << ex.what() << "\n";
        exit_code = 1;
    }

    PARALLEL::mpi_finalize();
    return exit_code;
}
