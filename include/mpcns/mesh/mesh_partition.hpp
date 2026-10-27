#pragma once
#include "mpcns/core/logical_index.hpp"

#include <array>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace mpcns::mesh
{
struct BlockAddress
{
    int rank = -1;
    int local_block = -1;
    bool valid() const noexcept { return rank >= 0 && local_block >= 0; }
    friend bool operator==(const BlockAddress &, const BlockAddress &) = default;
};

struct IndexTransform
{
    std::array<int, 3> permutation{};
    std::array<int, 3> sign{};
    core::LogicalIndex offset{};
    friend bool operator==(const IndexTransform &, const IndexTransform &) = default;
};

enum class PatchKind { Inner, Parallel, Physical };

struct LocalBlock
{
    BlockAddress address;
    std::string name;
    core::LogicalBox node_extent;
};

struct Patch
{
    PatchKind kind = PatchKind::Physical;
    BlockAddress this_block;
    BlockAddress neighbor_block; // invalid when has_neighbor is false
    std::string this_block_name;
    std::string neighbor_block_name;
    core::LogicalBox this_node_box;
    core::LogicalBox neighbor_node_box;
    // A transform exists only when a topological neighbor exists.  Physical
    // patches deliberately have std::nullopt rather than a misleading zero map.
    std::optional<IndexTransform> transform;
    bool has_neighbor = false;
    // Legacy coupling convention: data direction is neighbor -> this_block.
    // This is topology metadata only; it does not imply an MPI send/recv role.
    bool is_coupling = false;
    int direction = 0;
    int neighbor_direction = 0;
    std::array<int, 3> directions{};
    int boundary_id = 0;
    std::string boundary_name;
};

class MeshPartition
{
public:
    MeshPartition(int rank, int size, std::vector<LocalBlock> blocks,
                  std::vector<Patch> faces, std::vector<Patch> edges,
                  std::vector<Patch> vertices);

    int rank() const noexcept { return rank_; }
    int size() const noexcept { return size_; }
    const std::vector<LocalBlock> &local_blocks() const noexcept { return blocks_; }
    const std::vector<Patch> &face_patches() const noexcept { return faces_; }
    const std::vector<Patch> &edge_patches() const noexcept { return edges_; }
    const std::vector<Patch> &vertex_patches() const noexcept { return vertices_; }

private:
    int rank_;
    int size_;
    std::vector<LocalBlock> blocks_;
    std::vector<Patch> faces_;
    std::vector<Patch> edges_;
    std::vector<Patch> vertices_;
};
} // namespace mpcns::mesh
