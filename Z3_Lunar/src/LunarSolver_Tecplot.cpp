#include "LunarSolver.h"

#include "0_basic/MPI_WRAPPER.h"
#include "2_topology/Topology.h"
#include "3_field/Field.h"

#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
struct TecNodeFieldPair
{
    const char *cell_name;
    const char *node_name;
};
}

void LunarSolver::UpdateTecplotNodeFields_()
{
    if (!topo_ || !fld_)
        throw std::runtime_error("Lunar Tecplot node reconstruction needs topology and fields");

    const std::array<TecNodeFieldPair, 4> pairs{{
        {"PV_H", "PV_H_tecnode"},
        {"B_cell", "B_cell_tecnode"},
        {"Bind_cell", "Bind_cell_tecnode"},
        {"J_cell", "J_cell_tecnode"},
    }};

    const std::size_t node_count = topo_->nodes.rep_to_qid.size();
    if (node_count == 0)
        throw std::runtime_error("Lunar Tecplot node reconstruction found no quotient nodes");

    int rank = 0;
    PARALLEL::mpi_rank(&rank);

    for (const TecNodeFieldPair &pair : pairs)
    {
        const int cell_fid = fld_->field_id(pair.cell_name);
        const int node_fid = fld_->field_id(pair.node_name);
        const int ncomp = fld_->descriptor(cell_fid).ncomp;
        if (fld_->descriptor(node_fid).ncomp != ncomp)
            throw std::runtime_error(std::string("Tecplot node component mismatch for ") + pair.cell_name);

        std::vector<double> local_sum(node_count * static_cast<std::size_t>(ncomp), 0.0);
        std::vector<double> global_sum(local_sum.size(), 0.0);
        std::vector<double> local_weight(node_count, 0.0);
        std::vector<double> global_weight(node_count, 0.0);

        // Accumulate only real cells.  Each cell contributes to its eight
        // quotient nodes; shared panel aliases collapse to the same global id.
        for (int ib = 0; ib < fld_->num_blocks(); ++ib)
        {
            FieldBlock &cell = fld_->field(cell_fid, ib);
            FieldBlock &jac = fld_->field("Jac", ib);
            if (!cell.is_allocated() || !jac.is_allocated())
                continue;

            const Int3 lo = cell.inner_lo();
            const Int3 hi = cell.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const double weight = std::abs(jac(i, j, k, 0));
                        if (!std::isfinite(weight) || weight <= 0.0)
                            continue;

                        for (int di = 0; di <= 1; ++di)
                            for (int dj = 0; dj <= 1; ++dj)
                                for (int dk = 0; dk <= 1; ++dk)
                                {
                                    const TOPO::EntityKey node =
                                        TOPO::make_node(rank, ib, i + di, j + dj, k + dk);
                                    const int gid = topo_->id_of(node).id;
                                    if (gid < 0 || static_cast<std::size_t>(gid) >= node_count)
                                        throw std::runtime_error("Lunar Tecplot quotient-node id is out of range");
                                    local_weight[gid] += weight;
                                    const std::size_t base =
                                        static_cast<std::size_t>(gid) * static_cast<std::size_t>(ncomp);
                                    for (int comp = 0; comp < ncomp; ++comp)
                                        local_sum[base + comp] += weight * cell(i, j, k, comp);
                                }
                    }
        }

        MPI_Allreduce(local_sum.data(), global_sum.data(),
                      static_cast<int>(local_sum.size()), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(local_weight.data(), global_weight.data(),
                      static_cast<int>(local_weight.size()), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        // Materialize the same physical-node value into every local alias.
        for (int ib = 0; ib < fld_->num_blocks(); ++ib)
        {
            FieldBlock &node_field = fld_->field(node_fid, ib);
            if (!node_field.is_allocated())
                continue;
            const Int3 lo = node_field.inner_lo();
            const Int3 hi = node_field.inner_hi();
            for (int i = lo.i; i < hi.i; ++i)
                for (int j = lo.j; j < hi.j; ++j)
                    for (int k = lo.k; k < hi.k; ++k)
                    {
                        const TOPO::EntityKey node = TOPO::make_node(rank, ib, i, j, k);
                        const int gid = topo_->id_of(node).id;
                        const double weight = global_weight[gid];
                        if (!(weight > 0.0))
                            throw std::runtime_error("Lunar Tecplot node has no incident real cell");
                        const std::size_t base =
                            static_cast<std::size_t>(gid) * static_cast<std::size_t>(ncomp);
                        for (int comp = 0; comp < ncomp; ++comp)
                            node_field(i, j, k, comp) = global_sum[base + comp] / weight;
                    }
        }
    }
}
