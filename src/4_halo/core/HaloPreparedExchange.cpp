#include "4_halo/Halo.h"

#include "0_basic/Error.h"
#include "4_halo/detail/HaloBuildTools.h"

#include <algorithm>
#include <cstring>
#include <sstream>
#include <vector>

namespace
{
int box_length(const Box3 &box, int ncomp)
{
    const int ni = std::max(0, box.hi.i - box.lo.i);
    const int nj = std::max(0, box.hi.j - box.lo.j);
    const int nk = std::max(0, box.hi.k - box.lo.k);
    return ni * nj * nk * ncomp;
}

struct MatchRecord
{
    int rank;
    int peer;
    int tag;
    int length;
    int direction; // 0 send, 1 receive
    int ordinal;
    int enabled;
};

int message_ordinal(const std::vector<HALO_EXEC::Message> &messages,
                    std::size_t index)
{
    int ordinal = 0;
    for (std::size_t i = 0; i < index; ++i)
        if (messages[i].peer == messages[index].peer &&
            messages[i].tag == messages[index].tag)
            ++ordinal;
    return ordinal;
}

void collectively_match_plan(HALO_EXEC::Plan &plan)
{
    int rank = 0;
    int size = 1;
    PARALLEL::mpi_rank(&rank);
    PARALLEL::mpi_size(&size);

    std::vector<MatchRecord> local;
    local.reserve(plan.sends.size() + plan.recvs.size());
    for (std::size_t i = 0; i < plan.sends.size(); ++i)
        local.push_back({rank, plan.sends[i].peer, plan.sends[i].tag,
                         plan.sends[i].length, 0,
                         message_ordinal(plan.sends, i),
                         plan.sends[i].enabled ? 1 : 0});
    for (std::size_t i = 0; i < plan.recvs.size(); ++i)
        local.push_back({rank, plan.recvs[i].peer, plan.recvs[i].tag,
                         plan.recvs[i].length, 1,
                         message_ordinal(plan.recvs, i),
                         plan.recvs[i].enabled ? 1 : 0});

    const int local_bytes =
        static_cast<int>(local.size() * sizeof(MatchRecord));
    std::vector<int> byte_counts(size, 0);
    MPI_Allgather(&local_bytes, 1, MPI_INT, byte_counts.data(), 1,
                  MPI_INT, MPI_COMM_WORLD);
    std::vector<int> displacements(size, 0);
    int total_bytes = 0;
    for (int peer = 0; peer < size; ++peer)
    {
        displacements[peer] = total_bytes;
        total_bytes += byte_counts[peer];
    }
    std::vector<unsigned char> gathered(
        static_cast<std::size_t>(std::max(1, total_bytes)), 0);
    MPI_Allgatherv(local.empty() ? nullptr : local.data(), local_bytes,
                   MPI_BYTE, gathered.data(), byte_counts.data(),
                   displacements.data(), MPI_BYTE, MPI_COMM_WORLD);

    if (total_bytes % static_cast<int>(sizeof(MatchRecord)) != 0)
        ERROR::Abort("[Halo] corrupt prepared exchange schedule: " +
                     plan.name);
    const MatchRecord *records =
        reinterpret_cast<const MatchRecord *>(gathered.data());
    const std::size_t record_count =
        static_cast<std::size_t>(total_bytes) / sizeof(MatchRecord);

    const auto match = [&](HALO_EXEC::Message &message, int direction,
                           int ordinal)
    {
        const int opposite = direction == 0 ? 1 : 0;
        const MatchRecord *found = nullptr;
        for (std::size_t i = 0; i < record_count; ++i)
        {
            const MatchRecord &candidate = records[i];
            if (candidate.rank == message.peer &&
                candidate.peer == rank &&
                candidate.tag == message.tag &&
                candidate.direction == opposite &&
                candidate.ordinal == ordinal)
            {
                if (found)
                    ERROR::Abort("[Halo] duplicate prepared peer message: " +
                                 plan.name);
                found = &candidate;
            }
        }
        if (!found || found->length != message.length)
        {
            std::ostringstream stream;
            stream << "[Halo] prepared message mismatch plan=" << plan.name
                   << " rank=" << rank << " peer=" << message.peer
                   << " tag=" << message.tag
                   << " length=" << message.length;
            ERROR::Abort(stream.str());
        }
        message.enabled = message.enabled && found->enabled != 0;
    };

    for (std::size_t i = 0; i < plan.sends.size(); ++i)
        match(plan.sends[i], 0, message_ordinal(plan.sends, i));
    for (std::size_t i = 0; i < plan.recvs.size(); ++i)
        match(plan.recvs[i], 1, message_ordinal(plan.recvs, i));

    HALO_EXEC::finalize_plan(plan);
}

void pack_transformed(FieldBlock &field, const HaloRegion &region,
                      int ncomp, double *buffer, int length)
{
    if (!field.is_allocated() || length == 0)
        return;

    const Box3 &box = region.send_box;
    const int len_local[3] = {
        box.hi.i - box.lo.i,
        box.hi.j - box.lo.j,
        box.hi.k - box.lo.k};
    if (len_local[0] <= 0 || len_local[1] <= 0 || len_local[2] <= 0)
        return;

    const int local_lo[3] = {box.lo.i, box.lo.j, box.lo.k};
    const int local_hi[3] = {box.hi.i - 1, box.hi.j - 1,
                             box.hi.k - 1};
    const int offset[3] = {region.trans.offset.i,
                           region.trans.offset.j,
                           region.trans.offset.k};
    int target1[3], target2[3], target_ref[3], target_length[3] = {0,0,0};
    for (int axis = 0; axis < 3; ++axis)
    {
        target1[region.trans.perm[axis]] =
            region.trans.sign[axis] * local_lo[axis] + offset[axis];
        target2[region.trans.perm[axis]] =
            region.trans.sign[axis] * local_hi[axis] + offset[axis];
        target_length[region.trans.perm[axis]] = len_local[axis];
    }
    for (int axis = 0; axis < 3; ++axis)
        target_ref[axis] = std::min(target1[axis], target2[axis]);

    int source[3], target[3];
    for (source[0] = box.lo.i; source[0] < box.hi.i; ++source[0])
        for (source[1] = box.lo.j; source[1] < box.hi.j; ++source[1])
            for (source[2] = box.lo.k; source[2] < box.hi.k; ++source[2])
            {
                for (int axis = 0; axis < 3; ++axis)
                    target[region.trans.perm[axis]] =
                        region.trans.sign[axis] * source[axis] +
                        offset[axis];
                const int ri = target[0] - target_ref[0];
                const int rj = target[1] - target_ref[1];
                const int rk = target[2] - target_ref[2];
                const int base =
                    ((ri * target_length[1] + rj) * target_length[2] +
                     rk) * ncomp;
                if (base < 0 || base + ncomp > length)
                    ERROR::Abort("[Halo] prepared component pack is out of bounds");
                for (int component = 0; component < ncomp; ++component)
                    buffer[base + component] =
                        field(source[0], source[1], source[2], component);
            }
}

void unpack_box(FieldBlock &field, const HaloRegion &region, int ncomp,
                const double *buffer, int length)
{
    if (!field.is_allocated() || length == 0)
        return;
    const Box3 &box = region.recv_box;
    const int ni = box.hi.i - box.lo.i;
    const int nj = box.hi.j - box.lo.j;
    const int nk = box.hi.k - box.lo.k;
    if (ni * nj * nk * ncomp != length)
        ERROR::Abort("[Halo] prepared component unpack length mismatch");

    const Int3 lo = field.get_lo();
    const Int3 hi = field.get_hi();
    for (int ii = 0; ii < ni; ++ii)
        for (int jj = 0; jj < nj; ++jj)
            for (int kk = 0; kk < nk; ++kk)
            {
                const int i = box.lo.i + ii;
                const int j = box.lo.j + jj;
                const int k = box.lo.k + kk;
                if (i < lo.i || i >= hi.i || j < lo.j || j >= hi.j ||
                    k < lo.k || k >= hi.k)
                    ERROR::Abort("[Halo] prepared component recv box is out of bounds");
                const int base = ((ii * nj + jj) * nk + kk) * ncomp;
                for (int component = 0; component < ncomp; ++component)
                    field(i, j, k, component) =
                        buffer[base + component];
            }
}
}

void Halo::build_component_copy_exchange_plans_(std::uint64_t generation)
{
    component_exchange_plans_.clear();
    for (const std::string &field_name : component_copy_fields_)
    {
        const FieldHaloRequest &request = halo_request_(field_name);
        const int fid = fld_->field_id(field_name);
        for (int stage_value = static_cast<int>(HaloLevel::Corner1D);
             stage_value <= static_cast<int>(request.level); ++stage_value)
        {
            const HaloLevel stage = static_cast<HaloLevel>(stage_value);
            const PatternKey key{request.location, request.nghost};
            const HaloPattern *send_pattern = nullptr;
            const HaloPattern *recv_pattern = nullptr;
            if (stage == HaloLevel::Corner1D)
            {
                const auto it = parallel_patterns_.find(key);
                if (it != parallel_patterns_.end())
                    send_pattern = recv_pattern = &it->second;
            }
            else if (stage == HaloLevel::Corner2D)
            {
                const auto send = parallel_edge_patterns_send.find(key);
                const auto recv = parallel_edge_patterns_recv.find(key);
                if (send != parallel_edge_patterns_send.end())
                    send_pattern = &send->second;
                if (recv != parallel_edge_patterns_recv.end())
                    recv_pattern = &recv->second;
            }
            else
            {
                const auto send = parallel_vertex_patterns_send.find(key);
                const auto recv = parallel_vertex_patterns_recv.find(key);
                if (send != parallel_vertex_patterns_send.end())
                    send_pattern = &send->second;
                if (recv != parallel_vertex_patterns_recv.end())
                    recv_pattern = &recv->second;
            }

            ComponentExchangePlan prepared;
            prepared.fid = fid;
            prepared.ncomp = request.ncomp;
            prepared.stage = stage;
            prepared.exchange_plan.name = "component:" + field_name +
                                          ":" + std::to_string(stage_value);
            prepared.exchange_plan.generation = generation;
            if (send_pattern)
                prepared.send_regions = send_pattern->regions;
            if (recv_pattern)
                prepared.recv_regions = recv_pattern->regions;

            for (const HaloRegion &region : prepared.send_regions)
                prepared.exchange_plan.sends.push_back({
                    region.neighbor_rank, region.send_flag,
                    box_length(region.send_box, request.ncomp), 0,
                    field_active_(fid, region.this_block)});
            for (const HaloRegion &region : prepared.recv_regions)
                prepared.exchange_plan.recvs.push_back({
                    region.neighbor_rank, region.recv_flag,
                    box_length(region.recv_box, request.ncomp), 0,
                    field_active_(fid, region.this_block)});

            collectively_match_plan(prepared.exchange_plan);
            prepared_exchange_executor_.reserve_plan(
                prepared.exchange_plan);
            component_exchange_plans_.emplace(
                ComponentExchangeKey{field_name, stage_value},
                std::move(prepared));
        }
    }
}

const Halo::ComponentExchangePlan &Halo::component_exchange_plan_(
    const std::string &field_name, HaloLevel stage) const
{
    const auto found = component_exchange_plans_.find(
        ComponentExchangeKey{field_name, static_cast<int>(stage)});
    if (found == component_exchange_plans_.end())
        ERROR::Abort("[Halo] component exchange plan was not built: " +
                     field_name);
    return found->second;
}

void Halo::execute_component_exchange_plan_(
    const ComponentExchangePlan &prepared)
{
    prepared_exchange_executor_.execute(
        prepared.exchange_plan,
        [&](std::size_t index, double *buffer, int length)
        {
            const HaloRegion &region = prepared.send_regions[index];
            FieldBlock &field =
                fld_->field(prepared.fid, region.this_block);
            pack_transformed(field, region, prepared.ncomp, buffer, length);
        },
        [&](std::size_t index, const double *buffer, int length)
        {
            const HaloRegion &region = prepared.recv_regions[index];
            FieldBlock &field =
                fld_->field(prepared.fid, region.this_block);
            unpack_box(field, region, prepared.ncomp, buffer, length);
        });
}

std::string Halo::oriented_face_plan_key_(
    const std::string &group, FieldValueKind value_kind) const
{
    return (value_kind == FieldValueKind::EdgeCovariant1Form
                ? "edge1:"
                : "face2:") + group;
}

void Halo::build_oriented_face_exchange_plans_(std::uint64_t generation)
{
    oriented_face_exchange_plans_.clear();
    std::vector<std::string> edge_groups;
    std::vector<std::string> face_groups;
    for (const auto &entry : edge_1form_triplets_)
        edge_groups.push_back(entry.first);
    for (const auto &entry : face_2form_triplets_)
        face_groups.push_back(entry.first);
    std::sort(edge_groups.begin(), edge_groups.end());
    std::sort(face_groups.begin(), face_groups.end());

    const auto build_group =
        [&](const HaloTripletRequest &triplet, FieldValueKind kind)
    {
        const std::string names[3] = {
            triplet.xi, triplet.eta, triplet.zeta};
        const int fids[3] = {
            fld_->field_id(names[0]), fld_->field_id(names[1]),
            fld_->field_id(names[2])};
        const int ncomp = fld_->descriptor(fids[0]).ncomp;
        std::vector<OrientedFaceExchangePlan> plans;
        plans.reserve(3);

        for (int destination_axis = 0; destination_axis < 3;
             ++destination_axis)
        {
            const StaggerLocation destination_location =
                kind == FieldValueKind::EdgeCovariant1Form
                    ? (destination_axis == 0
                           ? StaggerLocation::EdgeXi
                           : destination_axis == 1
                                 ? StaggerLocation::EdgeEt
                                 : StaggerLocation::EdgeZe)
                    : (destination_axis == 0
                           ? StaggerLocation::FaceXi
                           : destination_axis == 1
                                 ? StaggerLocation::FaceEt
                                 : StaggerLocation::FaceZe);
            const PatternKey destination_key{
                destination_location, triplet.nghost};
            const auto destination =
                parallel_patterns_.find(destination_key);

            OrientedFaceExchangePlan prepared;
            prepared.ncomp = ncomp;
            prepared.exchange_plan.name =
                oriented_face_plan_key_(triplet.group_name, kind) +
                ":" + std::to_string(destination_axis);
            prepared.exchange_plan.generation = generation;

            if (destination != parallel_patterns_.end())
            {
                const std::vector<HaloRegion> &regions =
                    destination->second.regions;
                for (std::size_t index = 0; index < regions.size();
                     ++index)
                {
                    const HaloRegion &destination_region =
                        regions[index];
                    int source_axis = -1;
                    for (int axis = 0; axis < 3; ++axis)
                        if (destination_region.trans.perm[axis] ==
                            destination_axis)
                            source_axis = axis;
                    if (source_axis < 0)
                        ERROR::Abort("[Halo] invalid oriented Face-stage permutation");

                    const StaggerLocation source_location =
                        kind == FieldValueKind::EdgeCovariant1Form
                            ? (source_axis == 0
                                   ? StaggerLocation::EdgeXi
                                   : source_axis == 1
                                         ? StaggerLocation::EdgeEt
                                         : StaggerLocation::EdgeZe)
                            : (source_axis == 0
                                   ? StaggerLocation::FaceXi
                                   : source_axis == 1
                                         ? StaggerLocation::FaceEt
                                         : StaggerLocation::FaceZe);
                    const auto source = parallel_patterns_.find(
                        PatternKey{source_location, triplet.nghost});
                    if (source == parallel_patterns_.end() ||
                        index >= source->second.regions.size())
                        ERROR::Abort("[Halo] oriented Face-stage source pattern is incomplete");

                    OrientedFaceTransfer transfer;
                    transfer.destination_region = destination_region;
                    transfer.source_region = source->second.regions[index];
                    transfer.source_region.trans =
                        destination_region.trans;
                    transfer.source_fid = fids[source_axis];
                    transfer.destination_fid = fids[destination_axis];
                    transfer.sign =
                        kind == FieldValueKind::EdgeCovariant1Form
                            ? destination_region.trans.sign[source_axis]
                            : HALO_TOOLS::face_2form_orientation_sign(
                                  destination_region.trans, source_axis);
                    prepared.send_transfers.push_back(transfer);
                    prepared.recv_transfers.push_back(transfer);

                    HALO_EXEC::Message send_message;
                    send_message.peer =
                        destination_region.neighbor_rank;
                    send_message.tag = destination_region.send_flag;
                    send_message.length = box_length(
                        transfer.source_region.send_box, ncomp);
                    send_message.enabled = field_active_(
                        transfer.source_fid,
                        transfer.source_region.this_block);
                    prepared.exchange_plan.sends.push_back(send_message);

                    HALO_EXEC::Message recv_message;
                    recv_message.peer =
                        destination_region.neighbor_rank;
                    recv_message.tag = destination_region.recv_flag;
                    recv_message.length = box_length(
                        destination_region.recv_box, ncomp);
                    recv_message.enabled = field_active_(
                        transfer.destination_fid,
                        destination_region.this_block);
                    prepared.exchange_plan.recvs.push_back(recv_message);
                }
            }
            collectively_match_plan(prepared.exchange_plan);
            prepared_exchange_executor_.reserve_plan(
                prepared.exchange_plan);
            plans.push_back(std::move(prepared));
        }
        oriented_face_exchange_plans_.emplace(
            oriented_face_plan_key_(triplet.group_name, kind),
            std::move(plans));
    };

    for (const std::string &group : edge_groups)
        build_group(edge_1form_triplets_.at(group),
                    FieldValueKind::EdgeCovariant1Form);
    for (const std::string &group : face_groups)
        build_group(face_2form_triplets_.at(group),
                    FieldValueKind::FaceContravariant2Form);
}

void Halo::execute_oriented_face_exchange_plans_(
    const std::vector<std::string> &fields, FieldValueKind value_kind)
{
    if (fields.empty())
        ERROR::Abort("[Halo] empty oriented Face-stage group");
    const auto group = field_to_group_.find(fields[0]);
    if (group == field_to_group_.end())
        ERROR::Abort("[Halo] oriented Face-stage group is not registered");
    const auto found = oriented_face_exchange_plans_.find(
        oriented_face_plan_key_(group->second, value_kind));
    if (found == oriented_face_exchange_plans_.end())
        ERROR::Abort("[Halo] oriented Face-stage plan was not built");

    for (const OrientedFaceExchangePlan &prepared : found->second)
        prepared_exchange_executor_.execute(
            prepared.exchange_plan,
            [&](std::size_t index, double *buffer, int length)
            {
                const OrientedFaceTransfer &transfer =
                    prepared.send_transfers[index];
                FieldBlock &field = fld_->field(
                    transfer.source_fid,
                    transfer.source_region.this_block);
                pack_transformed(field, transfer.source_region,
                                 prepared.ncomp, buffer, length);
                if (transfer.sign < 0)
                    for (int value = 0; value < length; ++value)
                        buffer[value] = -buffer[value];
            },
            [&](std::size_t index, const double *buffer, int length)
            {
                const OrientedFaceTransfer &transfer =
                    prepared.recv_transfers[index];
                FieldBlock &field = fld_->field(
                    transfer.destination_fid,
                    transfer.destination_region.this_block);
                unpack_box(field, transfer.destination_region,
                           prepared.ncomp, buffer, length);
            });
}

HALO_EXEC::Statistics
Halo::debug_prepared_exchange_statistics() const noexcept
{
    return prepared_exchange_executor_.statistics();
}

void Halo::configure_owner_sync_tags_()
{
    owner_sync_tags_.clear();
    int *tag_upper_bound = nullptr;
    int present = 0;
    MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &tag_upper_bound,
                      &present);
    if (!present || tag_upper_bound == nullptr)
        ERROR::Abort("[Halo] MPI_TAG_UB is unavailable");

    const int required = static_cast<int>(owner_sync_requests_.size());
    if (required > *tag_upper_bound)
        ERROR::Abort("[Halo] insufficient MPI tags for Owner/Alias plans");

    // The managed Owner/Alias namespace is allocated from MPI_TAG_UB down.
    // Existing legacy halo tags are built from small local ordinals, so the
    // two namespaces cannot collide.
    for (int index = 0; index < required; ++index)
        owner_sync_tags_.emplace(owner_sync_requests_[index].field_name,
                                 *tag_upper_bound - index);
}

void Halo::validate_owner_sync_pattern_lengths_(OwnerSyncPattern &pat) const
{
    int size = 1;
    int rank = 0;
    PARALLEL::mpi_size(&size);
    PARALLEL::mpi_rank(&rank);
    std::vector<int> send_lengths(size, 0), recv_lengths(size, 0);
    std::vector<int> peer_send_lengths(size, 0);
    std::vector<int> send_counts(size, 0), recv_counts(size, 0);

    for (const OwnerSyncSendOp &op : pat.send_ops)
    {
        send_lengths[op.dst_rank] += op.ncomp;
        ++send_counts[op.dst_rank];
    }
    for (const OwnerSyncRecvOp &op : pat.recv_ops)
    {
        recv_lengths[op.src_rank] += op.ncomp;
        ++recv_counts[op.src_rank];
    }

    PARALLEL::mpi_alltoall(send_lengths.data(), 1,
                           peer_send_lengths.data(), 1);
    for (int peer = 0; peer < size; ++peer)
        if (recv_lengths[peer] != peer_send_lengths[peer])
        {
            std::ostringstream stream;
            stream << "[Halo] owner sync length mismatch field="
                   << pat.field_name << " rank=" << rank
                   << " peer=" << peer
                   << " recv=" << recv_lengths[peer]
                   << " send=" << peer_send_lengths[peer];
            ERROR::Abort(stream.str());
        }

    std::vector<int> expected_recv_counts(size, 0);
    PARALLEL::mpi_alltoall(send_counts.data(), 1,
                           expected_recv_counts.data(), 1);
    for (int peer = 0; peer < size; ++peer)
        if (recv_counts[peer] != expected_recv_counts[peer])
            ERROR::Abort("[Halo] owner sync operation-count mismatch: " +
                         pat.field_name);

    std::vector<int> send_displacements(size, 0);
    std::vector<int> recv_displacements(size, 0);
    int total_send = 0, total_recv = 0;
    for (int peer = 0; peer < size; ++peer)
    {
        send_displacements[peer] = total_send;
        recv_displacements[peer] = total_recv;
        total_send += send_counts[peer];
        total_recv += recv_counts[peer];
    }
    std::vector<int> send_active(total_send, 0);
    std::vector<int> recv_active(total_recv, 0);
    std::vector<int> cursor = send_displacements;
    for (const OwnerSyncSendOp &op : pat.send_ops)
        send_active[cursor[op.dst_rank]++] = op.active ? 1 : 0;

    MPI_Alltoallv(send_active.empty() ? nullptr : send_active.data(),
                  send_counts.data(), send_displacements.data(), MPI_INT,
                  recv_active.empty() ? nullptr : recv_active.data(),
                  recv_counts.data(), recv_displacements.data(), MPI_INT,
                  MPI_COMM_WORLD);
    cursor = recv_displacements;
    for (OwnerSyncRecvOp &op : pat.recv_ops)
        op.active = op.active &&
                    recv_active[cursor[op.src_rank]++] != 0;

    pat.mpi_lengths_validated = true;
}
