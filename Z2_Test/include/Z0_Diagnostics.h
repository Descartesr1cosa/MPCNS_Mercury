#pragma once

#include "2_topology/TopologyBuilder.h"
#include "2_topology/Topology.h"
#include "3_field/FieldDescriptor.h"

#include <iosfwd>
#include <string>

class Field;
class Halo;

namespace Z0
{
    struct TestResult
    {
        bool pass = true;
        double max_error = 0.0;
    };

    void print_banner();
    void print_diagnostics(const Field &fields,
                           const TOPO::Topology &topology,
                           const TOPO::Topology &equiv,
                           int dimension,
                           int nghost,
                           std::ostream &os);
    void dump_field_catalog(const Field &fields, int my_rank, std::ostream &os);
    void dump_topology_equiv_summary(const TOPO::Topology &equiv, int my_rank, std::ostream &os);
    void report_test(const std::string &name, const TestResult &result, std::ostream &os);
    void report_failure(const char *test,
                        const Field &fields,
                        const FieldDescriptor &desc,
                        int rank,
                        int block,
                        int i,
                        int j,
                        int k,
                        int comp,
                        double expected,
                        double actual,
                        std::ostream &os,
                        const TOPO::IndexTransform *tr = nullptr);
}
