#include "2_topology/TopologyEquivDetail.h"

#include "1_grid/1_MPCNS_Grid.h"

#include <algorithm>

namespace TOPO::detail
{
    std::vector<EntityKey> collect_all_local_cells(
        Grid &grid,
        int my_rank,
        int dimension)
    {
        std::vector<EntityKey> cells;

        if (dimension < 3)
        return cells;

        for (int ib = 0; ib < grid.nblock; ++ib)
        {
        const auto &blk = grid.grids(ib);

        for (int i = 0; i < blk.mx; ++i)
            for (int j = 0; j < blk.my; ++j)
            for (int k = 0; k < blk.mz; ++k)
            {
                cells.push_back(make_cell(my_rank, ib, i, j, k));
            }
        }
        return cells;
    }

    
    void build_cell_entity_ids(
        const std::vector<EntityKey> &all_local_cells,
        Topology &equiv)
    {
        std::vector<EntityKey> local_keys = all_local_cells;
        std::sort(local_keys.begin(), local_keys.end());
        local_keys.erase(std::unique(local_keys.begin(), local_keys.end()), local_keys.end());

        std::vector<int> send_buf;
        send_buf.reserve(local_keys.size() * 5);
        for (const EntityKey &cell : local_keys)
        pack_node(send_buf, make_node(cell.rank, cell.block, cell.i, cell.j, cell.k));

        const std::vector<int> recv_buf =
        allgather_packed_records(send_buf, 5, "build_topology cell EntityId");
        std::vector<EntityKey> global_keys;
        global_keys.reserve(recv_buf.size() / 5);
        for (std::size_t n = 0; n < recv_buf.size(); n += 5)
        {
        const EntityKey node = unpack_node(recv_buf.data() + n);
        global_keys.push_back(make_cell(node.rank, node.block, node.i, node.j, node.k));
        }

        std::sort(global_keys.begin(), global_keys.end());
        global_keys.erase(std::unique(global_keys.begin(), global_keys.end()), global_keys.end());
        equiv.cells.local_to_qid.clear();
        for (std::size_t n = 0; n < global_keys.size(); ++n)
        equiv.cells.local_to_qid[global_keys[n]] = static_cast<int>(n);
    }

    } // namespace TOPO::detail
