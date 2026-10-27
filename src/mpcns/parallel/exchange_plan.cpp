#include "mpcns/parallel/exchange_plan.hpp"

#include <stdexcept>
#include <utility>

namespace
{
using mpcns::parallel::ExchangeKind;
using mpcns::parallel::HaloStage;

bool stage_matches(HaloStage stage, ExchangeKind kind)
{
    if (kind == ExchangeKind::InnerCopy)
        return true;
    if (stage == HaloStage::Face)
        return kind == ExchangeKind::ParallelFace;
    if (stage == HaloStage::Edge)
        return kind == ExchangeKind::ParallelEdgeSend ||
               kind == ExchangeKind::ParallelEdgeRecv;
    return kind == ExchangeKind::ParallelVertexSend ||
           kind == ExchangeKind::ParallelVertexRecv;
}
} // namespace

namespace mpcns::parallel
{
ExchangePlan::ExchangePlan(FieldLocation location, int ghost_layers,
                           std::vector<ExchangeRegion> regions)
    : location_(location), ghost_layers_(ghost_layers), regions_(std::move(regions))
{
    if (ghost_layers < 0)
        throw std::invalid_argument("ExchangePlan ghost layers must be nonnegative");

    for (const ExchangeRegion &region : regions_)
    {
        if (region.location != location_ || region.ghost_layers != ghost_layers_)
            throw std::invalid_argument("ExchangePlan region key mismatch");
        if (!region.local_block.valid() || !region.peer_block.valid())
            throw std::invalid_argument("ExchangePlan requires valid block addresses");
        if (!region.transform)
            throw std::invalid_argument("ExchangePlan region requires an index transform");
        if (!stage_matches(region.stage, region.kind))
            throw std::invalid_argument("ExchangePlan stage/kind mismatch");

        const bool both = region.kind == ExchangeKind::InnerCopy ||
                          region.kind == ExchangeKind::ParallelFace;
        const bool send = region.kind == ExchangeKind::ParallelEdgeSend ||
                          region.kind == ExchangeKind::ParallelVertexSend;
        const bool recv = region.kind == ExchangeKind::ParallelEdgeRecv ||
                          region.kind == ExchangeKind::ParallelVertexRecv;

        if (both && (!region.send_box || !region.recv_box))
            throw std::invalid_argument("two-sided region requires both boxes");
        if (send && (!region.send_box || region.recv_box))
            throw std::invalid_argument("send region box roles invalid");
        if (recv && (region.send_box || !region.recv_box))
            throw std::invalid_argument("receive region box roles invalid");
        if (region.kind == ExchangeKind::InnerCopy &&
            (region.send_tag || region.recv_tag))
            throw std::invalid_argument("inner region cannot carry MPI tags");
        if (region.kind == ExchangeKind::ParallelFace &&
            (!region.send_tag || !region.recv_tag))
            throw std::invalid_argument("parallel face requires send and receive tags");
        if (send && (!region.send_tag || region.recv_tag))
            throw std::invalid_argument("send region tag roles invalid");
        if (recv && (region.send_tag || !region.recv_tag))
            throw std::invalid_argument("receive region tag roles invalid");
    }
}
} // namespace mpcns::parallel
