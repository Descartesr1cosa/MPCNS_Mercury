#include "mpcns/mesh/mesh_partition.hpp"

namespace mpcns::mesh
{
namespace
{
void validate_patch(const Patch &patch, int size, int local_blocks)
{
    if (!patch.this_block.valid() || patch.this_block.rank >= size ||
        patch.this_block.local_block >= local_blocks)
        throw std::invalid_argument("MeshPartition patch has invalid this block");
    if (patch.has_neighbor)
    {
        if (!patch.neighbor_block.valid() || patch.neighbor_block.rank >= size ||
            !patch.transform.has_value())
            throw std::invalid_argument("MeshPartition neighbor patch is incomplete");
    }
    else if (patch.neighbor_block.valid() || patch.transform.has_value() ||
             patch.kind != PatchKind::Physical)
    {
        throw std::invalid_argument("MeshPartition physical patch must have no neighbor");
    }
}
}
MeshPartition::MeshPartition(int rank, int size, std::vector<LocalBlock> blocks,
                             std::vector<Patch> faces, std::vector<Patch> edges,
                             std::vector<Patch> vertices)
    : rank_(rank), size_(size), blocks_(std::move(blocks)), faces_(std::move(faces)),
      edges_(std::move(edges)), vertices_(std::move(vertices))
{
    if (rank < 0 || size < 1 || rank >= size)
        throw std::invalid_argument("MeshPartition requires a valid MPI rank snapshot");
    for (int i = 0; i < static_cast<int>(blocks_.size()); ++i)
    {
        if (blocks_[i].address != BlockAddress{rank, i})
            throw std::invalid_argument("MeshPartition local blocks must be rank-local and ordered");
    }
    for (const Patch &p : faces_) validate_patch(p, size, static_cast<int>(blocks_.size()));
    for (const Patch &p : edges_) validate_patch(p, size, static_cast<int>(blocks_.size()));
    for (const Patch &p : vertices_) validate_patch(p, size, static_cast<int>(blocks_.size()));
}
} // namespace mpcns::mesh
