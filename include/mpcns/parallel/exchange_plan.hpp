#pragma once

#include "mpcns/core/logical_index.hpp"
#include "mpcns/mesh/mesh_partition.hpp"

#include <cstddef>
#include <optional>
#include <vector>

namespace mpcns::parallel
{
enum class FieldLocation
{
    Cell, Node, FaceXi, FaceEta, FaceZeta, EdgeXi, EdgeEta, EdgeZeta
};
enum class HaloStage { Face, Edge, Vertex };
enum class ExchangeKind
{
    InnerCopy,
    ParallelFace,
    ParallelEdgeSend,
    ParallelEdgeRecv,
    ParallelVertexSend,
    ParallelVertexRecv
};

struct ExchangeRegion
{
    FieldLocation location{};
    int ghost_layers = 0;
    HaloStage stage{};
    ExchangeKind kind{};
    mesh::BlockAddress local_block;
    mesh::BlockAddress peer_block;
    std::optional<core::LogicalBox> send_box;
    std::optional<core::LogicalBox> recv_box;
    std::optional<mesh::IndexTransform> transform;
    std::optional<int> send_tag;
    std::optional<int> recv_tag;
    std::size_t pattern_order = 0;
};

class ExchangePlan
{
public:
    ExchangePlan(FieldLocation location, int ghost_layers,
                 std::vector<ExchangeRegion> regions);

    FieldLocation location() const noexcept { return location_; }
    int ghost_layers() const noexcept { return ghost_layers_; }
    const std::vector<ExchangeRegion> &regions() const noexcept { return regions_; }

private:
    FieldLocation location_;
    int ghost_layers_;
    std::vector<ExchangeRegion> regions_;
};
} // namespace mpcns::parallel
