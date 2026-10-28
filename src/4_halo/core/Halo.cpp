#include "4_halo/Halo.h"
#include "4_halo/detail/HaloBuildTools.h"

#include <algorithm>
#include <array>
#include <sstream>

namespace
{
int edge1_box_volume(const Box3 &box)
{
    return (box.hi.i - box.lo.i) *
           (box.hi.j - box.lo.j) *
           (box.hi.k - box.lo.k);
}

bool edge1_valid_box(const Box3 &box)
{
    return box.hi.i >= box.lo.i &&
           box.hi.j >= box.lo.j &&
           box.hi.k >= box.lo.k;
}

bool edge1_valid_transform(const TOPO::IndexTransform &transform)
{
    std::array<bool, 3> seen = {{false, false, false}};
    for (int axis = 0; axis < 3; ++axis)
    {
        const int target = transform.perm[axis];
        if (target < 0 || target >= 3 || seen[target] ||
            (transform.sign[axis] != -1 &&
             transform.sign[axis] != 1))
            return false;
        seen[target] = true;
    }
    return true;
}

bool edge1_valid_transfer(
    const LegacyOrientedCornerTransferSnapshot &transfer,
    LegacyOrientedCornerRole expected_role,
    int rank, int size, int ncomp, int nghost)
{
    if (transfer.role != expected_role ||
        transfer.recv_rank < 0 || transfer.recv_rank >= size ||
        transfer.send_rank < 0 || transfer.send_rank >= size ||
        transfer.recv_block < 0 || transfer.send_block < 0 ||
        transfer.recv_axis < 0 || transfer.recv_axis >= 3 ||
        transfer.send_axis < 0 || transfer.send_axis >= 3 ||
        (transfer.sign != -1 && transfer.sign != 1) ||
        transfer.tag < 0 ||
        transfer.ncomp != ncomp ||
        transfer.nghost != nghost ||
        !edge1_valid_box(transfer.recv_box) ||
        !edge1_valid_box(transfer.node_box_on_send) ||
        !edge1_valid_transform(transfer.trans_recv_to_send) ||
        transfer.message_length !=
            edge1_box_volume(transfer.recv_box) * ncomp)
        return false;

    if (expected_role == LegacyOrientedCornerRole::Receive)
    {
        if (transfer.recv_rank != rank ||
            transfer.peer_rank != transfer.send_rank ||
            transfer.send_box.has_value())
            return false;
    }
    else
    {
        if (transfer.send_rank != rank ||
            transfer.peer_rank != transfer.recv_rank ||
            !transfer.send_box ||
            !edge1_valid_box(*transfer.send_box))
            return false;
    }
    return true;
}

bool edge1_valid_vertex_direction(int direction)
{
    return direction >= -3 && direction <= 3 && direction != 0;
}

bool edge1_valid_vertex_transfer(
    const LegacyOrientedVertexTransferSnapshot &transfer,
    LegacyOrientedCornerRole expected_role,
    int rank, int size, int ncomp, int nghost)
{
    bool seen_direction[3] = {false, false, false};
    const int directions[3] = {
        transfer.dir1_send, transfer.dir2_send,
        transfer.dir3_send};
    for (int direction : directions)
    {
        if (!edge1_valid_vertex_direction(direction))
            return false;
        const int axis = direction < 0 ? -direction - 1 : direction - 1;
        if (seen_direction[axis])
            return false;
        seen_direction[axis] = true;
    }

    if (transfer.role != expected_role ||
        transfer.recv_rank < 0 || transfer.recv_rank >= size ||
        transfer.send_rank < 0 || transfer.send_rank >= size ||
        transfer.recv_block < 0 || transfer.send_block < 0 ||
        transfer.recv_axis < 0 || transfer.recv_axis >= 3 ||
        transfer.send_axis < 0 || transfer.send_axis >= 3 ||
        transfer.send_axis !=
            transfer.trans_recv_to_send.perm[transfer.recv_axis] ||
        transfer.sign !=
            transfer.trans_recv_to_send.sign[transfer.recv_axis] ||
        (transfer.sign != -1 && transfer.sign != 1) ||
        transfer.tag < 0 ||
        transfer.ncomp != ncomp ||
        transfer.nghost != nghost ||
        !edge1_valid_box(transfer.recv_box) ||
        !edge1_valid_box(transfer.node_box_on_send) ||
        !edge1_valid_transform(transfer.trans_recv_to_send) ||
        transfer.message_length !=
            edge1_box_volume(transfer.recv_box) * ncomp)
        return false;

    if (expected_role == LegacyOrientedCornerRole::Receive)
    {
        return transfer.recv_rank == rank &&
               transfer.peer_rank == transfer.send_rank &&
               !transfer.send_box;
    }
    return transfer.send_rank == rank &&
           transfer.peer_rank == transfer.recv_rank &&
           transfer.send_box &&
           edge1_valid_box(*transfer.send_box);
}

template <class Transfer>
bool face2_valid_transfer_common(
    const Transfer &transfer,
    LegacyOrientedCornerRole expected_role,
    int rank, int size, int ncomp, int nghost)
{
    if (transfer.role != expected_role ||
        transfer.recv_rank < 0 || transfer.recv_rank >= size ||
        transfer.send_rank < 0 || transfer.send_rank >= size ||
        transfer.recv_block < 0 || transfer.send_block < 0 ||
        transfer.recv_axis < 0 || transfer.recv_axis >= 3 ||
        transfer.send_axis < 0 || transfer.send_axis >= 3 ||
        (transfer.sign != -1 && transfer.sign != 1) ||
        transfer.tag < 0 || transfer.ncomp != ncomp ||
        transfer.nghost != nghost ||
        !edge1_valid_box(transfer.recv_box) ||
        !edge1_valid_box(transfer.node_box_on_send) ||
        !edge1_valid_transform(transfer.trans_recv_to_send) ||
        transfer.message_length !=
            edge1_box_volume(transfer.recv_box) * ncomp)
        return false;

    if (transfer.trans_recv_to_send.perm[transfer.send_axis] !=
            transfer.recv_axis ||
        transfer.sign != HALO_TOOLS::face_2form_orientation_sign(
            transfer.trans_recv_to_send, transfer.send_axis))
        return false;

    if (expected_role == LegacyOrientedCornerRole::Receive)
        return transfer.recv_rank == rank &&
               transfer.peer_rank == transfer.send_rank &&
               !transfer.send_box;
    return transfer.send_rank == rank &&
           transfer.peer_rank == transfer.recv_rank &&
           transfer.send_box &&
           edge1_valid_box(*transfer.send_box);
}

bool face2_valid_edge_transfer(
    const LegacyFace2FormEdgeTransferSnapshot &transfer,
    LegacyOrientedCornerRole expected_role,
    int rank, int size, int ncomp, int nghost)
{
    if (!edge1_valid_vertex_direction(transfer.dir1_send) ||
        !edge1_valid_vertex_direction(transfer.dir2_send) ||
        std::abs(transfer.dir1_send) == std::abs(transfer.dir2_send))
        return false;
    return face2_valid_transfer_common(
        transfer, expected_role, rank, size, ncomp, nghost);
}

bool face2_valid_vertex_transfer(
    const LegacyFace2FormVertexTransferSnapshot &transfer,
    LegacyOrientedCornerRole expected_role,
    int rank, int size, int ncomp, int nghost)
{
    const int directions[3] = {
        transfer.dir1_send, transfer.dir2_send, transfer.dir3_send};
    bool seen[3] = {false, false, false};
    for (int direction : directions)
    {
        if (!edge1_valid_vertex_direction(direction))
            return false;
        const int axis = std::abs(direction) - 1;
        if (seen[axis])
            return false;
        seen[axis] = true;
    }
    return face2_valid_transfer_common(
        transfer, expected_role, rank, size, ncomp, nghost);
}


bool edge1_all_ranks_same_string(
    const std::string &local, int size)
{
    const int local_size = static_cast<int>(local.size());
    std::vector<int> sizes(size, 0);
    MPI_Allgather(
        &local_size, 1, MPI_INT, sizes.data(), 1, MPI_INT,
        MPI_COMM_WORLD);

    std::vector<int> displacements(size, 0);
    int total = 0;
    for (int rank = 0; rank < size; ++rank)
    {
        displacements[rank] = total;
        total += sizes[rank];
    }

    std::vector<char> gathered(total);
    MPI_Allgatherv(
        local.empty() ? nullptr : local.data(), local_size, MPI_CHAR,
        gathered.empty() ? nullptr : gathered.data(), sizes.data(),
        displacements.data(), MPI_CHAR, MPI_COMM_WORLD);

    for (int rank = 0; rank < size; ++rank)
    {
        const std::string remote(
            gathered.data() + displacements[rank],
            static_cast<std::size_t>(sizes[rank]));
        if (remote != local)
            return false;
    }
    return true;
}

bool edge1_all_ranks_true(bool local_valid)
{
    const int local = local_valid ? 1 : 0;
    int global = 0;
    MPI_Allreduce(
        &local, &global, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    return global == 1;
}
} // namespace


Halo::Halo(Field *field, TOPO::Topology *topo)
    : fld_(field), topo_(topo)
{
}

void Halo::set_topology_equiv(const TOPO::Topology *equiv)
{
    equiv_ = equiv;
}

const TOPO::Topology *Halo::topology_equiv() const
{
    return equiv_;
}

std::vector<HaloRegion> Halo::debug_halo_regions(StaggerLocation location,
                                                 int nghost,
                                                 HaloLevel stage) const
{
    const PatternKey key = {location, nghost};
    std::vector<HaloRegion> regions;

    auto append = [&regions](const std::map<PatternKey, HaloPattern> &patterns,
                             const PatternKey &pattern_key)
    {
        auto it = patterns.find(pattern_key);
        if (it == patterns.end())
            return;
        regions.insert(regions.end(), it->second.regions.begin(), it->second.regions.end());
    };

    if (stage == HaloLevel::Corner1D)
    {
        append(inner_patterns_, key);
        append(parallel_patterns_, key);
    }
    else if (stage == HaloLevel::Corner2D)
    {
        append(inner_edge_patterns_, key);
        append(parallel_edge_patterns_recv, key);
    }
    else if (stage == HaloLevel::Corner3D)
    {
        append(inner_vertex_patterns_, key);
        append(parallel_vertex_patterns_recv, key);
    }

    return regions;
}

std::vector<HaloRegion> Halo::debug_halo_send_regions(StaggerLocation location,
                                                       int nghost,
                                                       HaloLevel stage) const
{
    const PatternKey key = {location, nghost};
    std::vector<HaloRegion> regions;

    auto append = [&regions](const std::map<PatternKey, HaloPattern> &patterns,
                             const PatternKey &pattern_key)
    {
        auto it = patterns.find(pattern_key);
        if (it == patterns.end())
            return;
        regions.insert(regions.end(), it->second.regions.begin(), it->second.regions.end());
    };

    if (stage == HaloLevel::Corner1D)
    {
        append(parallel_patterns_, key);
    }
    else if (stage == HaloLevel::Corner2D)
    {
        append(parallel_edge_patterns_send, key);
    }
    else if (stage == HaloLevel::Corner3D)
    {
        append(parallel_vertex_patterns_send, key);
    }

    return regions;
}

std::vector<LegacyHaloPatternSnapshot>
Halo::debug_component_copy_patterns(StaggerLocation location, int nghost) const
{
    const PatternKey key{location, nghost};
    std::vector<LegacyHaloPatternSnapshot> result;
    const auto append = [&](const auto &patterns, LegacyHaloPatternKind kind, HaloLevel stage)
    {
        const auto it = patterns.find(key);
        if (it != patterns.end())
            result.push_back({kind, stage, it->second.location, it->second.nghost,
                              it->second.regions});
    };
    append(inner_patterns_, LegacyHaloPatternKind::InnerCopy, HaloLevel::Corner1D);
    append(parallel_patterns_, LegacyHaloPatternKind::ParallelFace, HaloLevel::Corner1D);
    append(inner_edge_patterns_, LegacyHaloPatternKind::InnerCopy, HaloLevel::Corner2D);
    append(parallel_edge_patterns_send, LegacyHaloPatternKind::ParallelEdgeSend, HaloLevel::Corner2D);
    append(parallel_edge_patterns_recv, LegacyHaloPatternKind::ParallelEdgeRecv, HaloLevel::Corner2D);
    append(inner_vertex_patterns_, LegacyHaloPatternKind::InnerCopy, HaloLevel::Corner3D);
    append(parallel_vertex_patterns_send, LegacyHaloPatternKind::ParallelVertexSend, HaloLevel::Corner3D);
    append(parallel_vertex_patterns_recv, LegacyHaloPatternKind::ParallelVertexRecv, HaloLevel::Corner3D);
    return result;
}

std::vector<LegacyOrientedTripletSnapshot>
Halo::debug_oriented_form_triplets() const
{
    std::vector<LegacyOrientedTripletSnapshot> result;
    const auto append = [&](const auto &groups, LegacyOrientedFormKind kind)
    {
        for (const auto &entry : groups)
        {
            const HaloTripletRequest &triplet = entry.second;
            LegacyOrientedTripletSnapshot snapshot;
            snapshot.group_name = triplet.group_name;
            snapshot.form_kind = kind;
            snapshot.field_names = {triplet.xi, triplet.eta, triplet.zeta};
            snapshot.level = triplet.level;
            snapshot.orientation_aware = triplet.orientation_aware;
            for (int axis = 0; axis < 3; ++axis)
            {
                const FieldDescriptor &descriptor =
                    fld_->descriptor(snapshot.field_names[axis]);
                snapshot.locations[axis] = descriptor.location;
                snapshot.value_kinds[axis] = descriptor.value_kind;
                snapshot.ncomp[axis] = descriptor.ncomp;
                snapshot.nghost[axis] = descriptor.nghost;
            }
            result.push_back(std::move(snapshot));
        }
    };

    append(edge_1form_triplets_, LegacyOrientedFormKind::Edge1Form);
    append(face_2form_triplets_, LegacyOrientedFormKind::Face2Form);
    return result;
}

void Halo::enable_edge_1form_corner_plan_capture(bool enabled)
{
    if (enabled)
        edge_1form_corner_plan_snapshots_.clear();
    capture_edge_1form_corner_plans_ = enabled;
}

std::optional<LegacyOrientedCornerPlanSnapshot>
Halo::debug_edge_1form_corner_plan(
    const std::string &group_name) const
{
    const auto it =
        edge_1form_corner_plan_snapshots_.find(group_name);
    if (it == edge_1form_corner_plan_snapshots_.end())
        return std::nullopt;
    return it->second;
}

void Halo::invalidate_prepared_edge_1form_corner_metadata_()
{
    invalidate_prepared_edge_1form_edge_execution_();
    prepared_edge_1form_corner_plans_.clear();
    prepared_edge_1form_corner_metadata_ready_ = false;
    prepared_edge_1form_corner_registry_generation_ = 0;
    prepared_edge_1form_corner_pattern_generation_ = 0;
}

bool Halo::prepared_edge_1form_corner_metadata_valid_() const noexcept
{
    return prepared_edge_1form_corner_metadata_ready_ &&
           registered_patterns_built_ &&
           prepared_edge_1form_corner_registry_generation_ ==
               registry_generation_ &&
           prepared_edge_1form_corner_pattern_generation_ ==
               pattern_generation_;
}
void Halo::invalidate_prepared_edge_1form_edge_execution_()
{
    prepared_edge_1form_edge_execution_plans_.clear();
    edge_1form_edge_execution_registry_generation_ = 0;
    edge_1form_edge_execution_pattern_generation_ = 0;
}

bool Halo::prepared_edge_1form_edge_execution_valid_() const noexcept
{
    return edge_1form_edge_execution_mode_ ==
               Edge1FormEdgeExecutionMode::Prepared &&
           prepared_edge_1form_corner_metadata_valid_() &&
           edge_1form_edge_execution_registry_generation_ ==
               registry_generation_ &&
           edge_1form_edge_execution_pattern_generation_ ==
               pattern_generation_ &&
           !prepared_edge_1form_edge_execution_plans_.empty();
}



std::optional<LegacyOrientedCornerPlanSnapshot>
Halo::debug_prepared_edge_1form_corner_plan(
    const std::string &group_name) const
{
    if (!prepared_edge_1form_corner_metadata_valid_())
        return std::nullopt;

    const auto it =
        prepared_edge_1form_corner_plans_.find(group_name);
    if (it == prepared_edge_1form_corner_plans_.end() ||
        it->second.registry_generation != registry_generation_ ||
        it->second.pattern_generation != pattern_generation_)
        return std::nullopt;
    return it->second.snapshot;
}

std::optional<std::uint64_t>
Halo::debug_prepared_edge_1form_corner_plan_generation() const
{
    if (!prepared_edge_1form_corner_metadata_valid_())
        return std::nullopt;
    return prepared_edge_1form_corner_pattern_generation_;
}

std::uint64_t
Halo::debug_prepared_edge_1form_corner_metadata_exchange_count() const noexcept
{
    return prepared_edge_1form_corner_metadata_exchange_count_;
}

bool Halo::set_edge_1form_edge_execution_mode(
    Edge1FormEdgeExecutionMode mode)
{
    int rank = 0;
    int size = 1;
    PARALLEL::mpi_rank(&rank);
    PARALLEL::mpi_size(&size);

    const int requested_mode = static_cast<int>(mode);
    std::vector<int> requested_modes(size, 0);
    MPI_Allgather(
        &requested_mode, 1, MPI_INT,
        requested_modes.data(), 1, MPI_INT, MPI_COMM_WORLD);
    const bool same_mode = std::all_of(
        requested_modes.begin(), requested_modes.end(),
        [&](int remote)
        {
            return remote == requested_mode;
        });
    if (!same_mode)
    {
        invalidate_prepared_edge_1form_edge_execution_();
        edge_1form_edge_execution_mode_ =
            Edge1FormEdgeExecutionMode::Legacy;
        return false;
    }

    std::vector<std::string> groups;
    groups.reserve(edge_1form_triplets_.size());
    for (const auto &entry : edge_1form_triplets_)
    {
        if (static_cast<int>(entry.second.level) >=
            static_cast<int>(HaloLevel::Corner2D))
            groups.push_back(entry.first);
    }
    std::sort(groups.begin(), groups.end());

    std::ostringstream descriptor;
    descriptor << "edge1-edge-execution-v1";
    for (const std::string &group_name : groups)
    {
        const HaloTripletRequest &triplet =
            edge_1form_triplets_.at(group_name);
        descriptor << "|group=" << group_name
                   << "|stage=" <<
            static_cast<int>(HaloLevel::Corner2D)
                   << "|nghost=" << triplet.nghost;
        const std::string fields[3] = {
            triplet.xi, triplet.eta, triplet.zeta};
        for (int axis = 0; axis < 3; ++axis)
        {
            const FieldDescriptor &field =
                fld_->descriptor(fields[axis]);
            descriptor << "|field=" << fields[axis]
                       << "|location=" <<
                static_cast<int>(field.location)
                       << "|kind=" <<
                static_cast<int>(field.value_kind)
                       << "|ncomp=" << field.ncomp
                       << "|field_nghost=" << field.nghost
                       << "|halo=" << field.sync.do_halo
                       << "|oriented=" <<
                field.sync.orientation_aware;
        }
    }

    const bool same_descriptor =
        edge1_all_ranks_same_string(descriptor.str(), size);
    if (!same_descriptor)
    {
        invalidate_prepared_edge_1form_edge_execution_();
        edge_1form_edge_execution_mode_ =
            Edge1FormEdgeExecutionMode::Legacy;
        return false;
    }

    if (mode == Edge1FormEdgeExecutionMode::Legacy)
    {
        invalidate_prepared_edge_1form_edge_execution_();
        edge_1form_edge_execution_mode_ =
            Edge1FormEdgeExecutionMode::Legacy;
        return true;
    }

    bool local_valid =
        prepared_edge_1form_corner_metadata_valid_() &&
        !groups.empty();
    for (const std::string &group_name : groups)
    {
        const auto entry =
            prepared_edge_1form_corner_plans_.find(group_name);
        if (entry == prepared_edge_1form_corner_plans_.end() ||
            entry->second.registry_generation !=
                registry_generation_ ||
            entry->second.pattern_generation !=
                pattern_generation_)
        {
            local_valid = false;
            continue;
        }

        const LegacyOrientedCornerPlanSnapshot &snapshot =
            entry->second.snapshot;
        const HaloTripletRequest &triplet =
            edge_1form_triplets_.at(group_name);
        const std::array<std::string, 3> fields = {
            triplet.xi, triplet.eta, triplet.zeta};
        if (snapshot.group_name != group_name ||
            snapshot.field_names != fields ||
            snapshot.ncomp <= 0 ||
            snapshot.nghost != triplet.nghost)
        {
            local_valid = false;
            continue;
        }

        for (const auto &transfer : snapshot.recv_plan)
        {
            if (!edge1_valid_transfer(
                    transfer, LegacyOrientedCornerRole::Receive,
                    rank, size, snapshot.ncomp, snapshot.nghost))
                local_valid = false;
        }
        for (const auto &transfer : snapshot.send_plan)
        {
            if (!edge1_valid_transfer(
                    transfer, LegacyOrientedCornerRole::Send,
                    rank, size, snapshot.ncomp, snapshot.nghost))
                local_valid = false;
        }
    }

    if (!edge1_all_ranks_true(local_valid))
    {
        invalidate_prepared_edge_1form_edge_execution_();
        edge_1form_edge_execution_mode_ =
            Edge1FormEdgeExecutionMode::Legacy;
        return false;
    }

    invalidate_prepared_edge_1form_edge_execution_();
    for (const std::string &group_name : groups)
    {
        const auto entry =
            prepared_edge_1form_corner_plans_.find(group_name);
        prepared_edge_1form_edge_execution_plans_.emplace(
            group_name,
            make_prepared_edge_1form_edge_execution_view_(
                entry->second.snapshot));
    }
    reserve_prepared_exchange_views_();
    edge_1form_edge_execution_registry_generation_ =
        registry_generation_;
    edge_1form_edge_execution_pattern_generation_ =
        pattern_generation_;
    edge_1form_edge_execution_mode_ =
        Edge1FormEdgeExecutionMode::Prepared;
    return true;
}

bool Halo::activate_prepared_edge_1form_edge_execution()
{
    return set_edge_1form_edge_execution_mode(
        Edge1FormEdgeExecutionMode::Prepared);
}

bool Halo::deactivate_prepared_edge_1form_edge_execution()
{
    return set_edge_1form_edge_execution_mode(
        Edge1FormEdgeExecutionMode::Legacy);
}

Edge1FormEdgeExecutionMode
Halo::debug_edge_1form_edge_execution_mode() const noexcept
{
    return edge_1form_edge_execution_mode_;
}

bool Halo::debug_prepared_edge_1form_edge_execution_active() const noexcept
{
    return prepared_edge_1form_edge_execution_valid_();
}

Edge1FormEdgeExecutionCounters
Halo::debug_edge_1form_edge_execution_counters() const noexcept
{
    return edge_1form_edge_execution_counters_;
}


void Halo::invalidate_prepared_vertex_edge_1form_metadata_()
{
    invalidate_prepared_edge_1form_vertex_execution_();
    prepared_vertex_edge_1form_plans_.clear();
    prepared_vertex_edge_1form_metadata_ready_ = false;
    prepared_vertex_edge_1form_registry_generation_ = 0;
    prepared_vertex_edge_1form_pattern_generation_ = 0;
}

bool Halo::prepared_vertex_edge_1form_metadata_valid_() const noexcept
{
    return prepared_vertex_edge_1form_metadata_ready_ &&
           registered_patterns_built_ &&
           prepared_vertex_edge_1form_registry_generation_ ==
               registry_generation_ &&
           prepared_vertex_edge_1form_pattern_generation_ ==
               pattern_generation_;
}

std::optional<LegacyOrientedVertexPlanSnapshot>
Halo::debug_prepared_vertex_edge_1form_plan(
    const std::string &group_name) const
{
    if (!prepared_vertex_edge_1form_metadata_valid_())
        return std::nullopt;

    const auto it =
        prepared_vertex_edge_1form_plans_.find(group_name);
    if (it == prepared_vertex_edge_1form_plans_.end() ||
        it->second.registry_generation != registry_generation_ ||
        it->second.pattern_generation != pattern_generation_)
        return std::nullopt;
    return it->second.snapshot;
}

std::optional<std::uint64_t>
Halo::debug_prepared_vertex_edge_1form_plan_generation() const
{
    if (!prepared_vertex_edge_1form_metadata_valid_())
        return std::nullopt;
    return prepared_vertex_edge_1form_pattern_generation_;
}

std::uint64_t
Halo::debug_prepared_vertex_edge_1form_metadata_exchange_count() const noexcept
{
    return prepared_vertex_edge_1form_metadata_exchange_count_;
}

void Halo::invalidate_prepared_edge_1form_vertex_execution_()
{
    prepared_edge_1form_vertex_execution_plans_.clear();
    edge_1form_vertex_execution_registry_generation_ = 0;
    edge_1form_vertex_execution_pattern_generation_ = 0;
}

bool Halo::prepared_edge_1form_vertex_execution_valid_() const noexcept
{
    return edge_1form_vertex_execution_mode_ ==
               Edge1FormVertexExecutionMode::Prepared &&
           prepared_vertex_edge_1form_metadata_valid_() &&
           edge_1form_vertex_execution_registry_generation_ ==
               registry_generation_ &&
           edge_1form_vertex_execution_pattern_generation_ ==
               pattern_generation_ &&
           !prepared_edge_1form_vertex_execution_plans_.empty();
}

bool Halo::set_edge_1form_vertex_execution_mode(
    Edge1FormVertexExecutionMode mode)
{
    int rank = 0;
    int size = 1;
    PARALLEL::mpi_rank(&rank);
    PARALLEL::mpi_size(&size);

    const int requested_mode = static_cast<int>(mode);
    std::vector<int> requested_modes(size, 0);
    MPI_Allgather(
        &requested_mode, 1, MPI_INT,
        requested_modes.data(), 1, MPI_INT, MPI_COMM_WORLD);
    const bool same_mode = std::all_of(
        requested_modes.begin(), requested_modes.end(),
        [&](int remote)
        {
            return remote == requested_mode;
        });
    if (!same_mode)
    {
        invalidate_prepared_edge_1form_vertex_execution_();
        edge_1form_vertex_execution_mode_ =
            Edge1FormVertexExecutionMode::Legacy;
        return false;
    }

    std::vector<std::string> groups;
    groups.reserve(edge_1form_triplets_.size());
    for (const auto &entry : edge_1form_triplets_)
    {
        if (static_cast<int>(entry.second.level) >=
            static_cast<int>(HaloLevel::Corner3D))
            groups.push_back(entry.first);
    }
    std::sort(groups.begin(), groups.end());

    std::ostringstream descriptor;
    descriptor << "edge1-vertex-execution-v1";
    for (const std::string &group_name : groups)
    {
        const HaloTripletRequest &triplet =
            edge_1form_triplets_.at(group_name);
        descriptor << "|group=" << group_name
                   << "|level=" <<
            static_cast<int>(triplet.level)
                   << "|stage=" <<
            static_cast<int>(HaloLevel::Corner3D)
                   << "|nghost=" << triplet.nghost
                   << "|triplet_oriented=" <<
            triplet.orientation_aware;
        const std::string fields[3] = {
            triplet.xi, triplet.eta, triplet.zeta};
        for (int axis = 0; axis < 3; ++axis)
        {
            const FieldDescriptor &field =
                fld_->descriptor(fields[axis]);
            descriptor << "|field=" << fields[axis]
                       << "|location=" <<
                static_cast<int>(field.location)
                       << "|kind=" <<
                static_cast<int>(field.value_kind)
                       << "|ncomp=" << field.ncomp
                       << "|field_nghost=" << field.nghost
                       << "|halo=" << field.sync.do_halo
                       << "|oriented=" <<
                field.sync.orientation_aware;
        }
    }

    if (!edge1_all_ranks_same_string(descriptor.str(), size))
    {
        invalidate_prepared_edge_1form_vertex_execution_();
        edge_1form_vertex_execution_mode_ =
            Edge1FormVertexExecutionMode::Legacy;
        return false;
    }

    if (mode == Edge1FormVertexExecutionMode::Legacy)
    {
        invalidate_prepared_edge_1form_vertex_execution_();
        edge_1form_vertex_execution_mode_ =
            Edge1FormVertexExecutionMode::Legacy;
        return true;
    }

    const StaggerLocation expected_locations[3] = {
        StaggerLocation::EdgeXi,
        StaggerLocation::EdgeEt,
        StaggerLocation::EdgeZe};
    bool local_valid =
        prepared_vertex_edge_1form_metadata_valid_() &&
        !groups.empty();
    for (const std::string &group_name : groups)
    {
        const auto entry =
            prepared_vertex_edge_1form_plans_.find(group_name);
        if (entry == prepared_vertex_edge_1form_plans_.end() ||
            entry->second.registry_generation !=
                registry_generation_ ||
            entry->second.pattern_generation !=
                pattern_generation_)
        {
            local_valid = false;
            continue;
        }

        const LegacyOrientedVertexPlanSnapshot &snapshot =
            entry->second.snapshot;
        const HaloTripletRequest &triplet =
            edge_1form_triplets_.at(group_name);
        const std::array<std::string, 3> fields = {
            triplet.xi, triplet.eta, triplet.zeta};
        if (snapshot.group_name != group_name ||
            snapshot.field_names != fields ||
            snapshot.ncomp <= 0 ||
            snapshot.nghost != triplet.nghost ||
            static_cast<int>(triplet.level) <
                static_cast<int>(HaloLevel::Corner3D) ||
            !triplet.orientation_aware)
        {
            local_valid = false;
            continue;
        }

        for (int axis = 0; axis < 3; ++axis)
        {
            const FieldDescriptor &field =
                fld_->descriptor(fields[axis]);
            if (field.location != expected_locations[axis] ||
                field.value_kind !=
                    FieldValueKind::EdgeCovariant1Form ||
                field.ncomp != snapshot.ncomp ||
                field.nghost != snapshot.nghost ||
                !field.sync.do_halo ||
                !field.sync.orientation_aware)
                local_valid = false;
        }
        for (const auto &transfer : snapshot.recv_plan)
        {
            if (!edge1_valid_vertex_transfer(
                    transfer, LegacyOrientedCornerRole::Receive,
                    rank, size, snapshot.ncomp, snapshot.nghost))
                local_valid = false;
        }
        for (const auto &transfer : snapshot.send_plan)
        {
            if (!edge1_valid_vertex_transfer(
                    transfer, LegacyOrientedCornerRole::Send,
                    rank, size, snapshot.ncomp, snapshot.nghost))
                local_valid = false;
        }
    }

    if (!edge1_all_ranks_true(local_valid))
    {
        invalidate_prepared_edge_1form_vertex_execution_();
        edge_1form_vertex_execution_mode_ =
            Edge1FormVertexExecutionMode::Legacy;
        return false;
    }

    invalidate_prepared_edge_1form_vertex_execution_();
    for (const std::string &group_name : groups)
    {
        const auto entry =
            prepared_vertex_edge_1form_plans_.find(group_name);
        prepared_edge_1form_vertex_execution_plans_.emplace(
            group_name,
            make_prepared_edge_1form_vertex_execution_view_(
                entry->second.snapshot));
    }
    reserve_prepared_exchange_views_();
    edge_1form_vertex_execution_registry_generation_ =
        registry_generation_;
    edge_1form_vertex_execution_pattern_generation_ =
        pattern_generation_;
    edge_1form_vertex_execution_mode_ =
        Edge1FormVertexExecutionMode::Prepared;
    return true;
}

bool Halo::activate_prepared_edge_1form_vertex_execution()
{
    return set_edge_1form_vertex_execution_mode(
        Edge1FormVertexExecutionMode::Prepared);
}

bool Halo::deactivate_prepared_edge_1form_vertex_execution()
{
    return set_edge_1form_vertex_execution_mode(
        Edge1FormVertexExecutionMode::Legacy);
}

Edge1FormVertexExecutionMode
Halo::debug_edge_1form_vertex_execution_mode() const noexcept
{
    return edge_1form_vertex_execution_mode_;
}

bool Halo::debug_prepared_edge_1form_vertex_execution_active() const noexcept
{
    return prepared_edge_1form_vertex_execution_valid_();
}

Edge1FormVertexExecutionCounters
Halo::debug_edge_1form_vertex_execution_counters() const noexcept
{
    return edge_1form_vertex_execution_counters_;
}


void Halo::invalidate_prepared_edge_face_2form_metadata_()
{
    invalidate_prepared_face_2form_edge_execution_();
    prepared_edge_face_2form_plans_.clear();
    prepared_edge_face_2form_metadata_ready_ = false;
    prepared_edge_face_2form_registry_generation_ = 0;
    prepared_edge_face_2form_pattern_generation_ = 0;
}

bool Halo::prepared_edge_face_2form_metadata_valid_() const noexcept
{
    return prepared_edge_face_2form_metadata_ready_ &&
           registered_patterns_built_ &&
           prepared_edge_face_2form_registry_generation_ ==
               registry_generation_ &&
           prepared_edge_face_2form_pattern_generation_ ==
               pattern_generation_;
}

std::optional<LegacyFace2FormEdgePlanSnapshot>
Halo::debug_prepared_edge_face_2form_plan(
    const std::string &group_name) const
{
    if (!prepared_edge_face_2form_metadata_valid_())
        return std::nullopt;

    const auto it = prepared_edge_face_2form_plans_.find(group_name);
    if (it == prepared_edge_face_2form_plans_.end() ||
        it->second.registry_generation != registry_generation_ ||
        it->second.pattern_generation != pattern_generation_)
        return std::nullopt;
    return it->second.snapshot;
}

std::optional<std::uint64_t>
Halo::debug_prepared_edge_face_2form_plan_generation() const
{
    if (!prepared_edge_face_2form_metadata_valid_())
        return std::nullopt;
    return prepared_edge_face_2form_pattern_generation_;
}

std::uint64_t
Halo::debug_prepared_edge_face_2form_metadata_exchange_count() const noexcept
{
    return prepared_edge_face_2form_metadata_exchange_count_;
}
void Halo::invalidate_prepared_face_2form_edge_execution_()
{
    prepared_face_2form_edge_execution_plans_.clear();
    face_2form_edge_execution_registry_generation_ = 0;
    face_2form_edge_execution_pattern_generation_ = 0;
}

bool Halo::prepared_face_2form_edge_execution_valid_() const noexcept
{
    return face_2form_edge_execution_mode_ ==
               Face2FormEdgeExecutionMode::Prepared &&
           prepared_edge_face_2form_metadata_valid_() &&
           face_2form_edge_execution_registry_generation_ ==
               registry_generation_ &&
           face_2form_edge_execution_pattern_generation_ ==
               pattern_generation_ &&
           !prepared_face_2form_edge_execution_plans_.empty();
}

bool Halo::set_face_2form_edge_execution_mode(
    Face2FormEdgeExecutionMode mode)
{
    int rank = 0;
    int size = 1;
    PARALLEL::mpi_rank(&rank);
    PARALLEL::mpi_size(&size);

    const int requested_mode = static_cast<int>(mode);
    std::vector<int> requested_modes(size, 0);
    MPI_Allgather(
        &requested_mode, 1, MPI_INT, requested_modes.data(), 1,
        MPI_INT, MPI_COMM_WORLD);
    const bool same_mode = std::all_of(
        requested_modes.begin(), requested_modes.end(),
        [&](int remote) { return remote == requested_mode; });
    if (!same_mode)
    {
        invalidate_prepared_face_2form_edge_execution_();
        face_2form_edge_execution_mode_ =
            Face2FormEdgeExecutionMode::Legacy;
        return false;
    }

    std::vector<std::string> groups;
    for (const auto &entry : face_2form_triplets_)
    {
        if (entry.second.orientation_aware &&
            static_cast<int>(entry.second.level) >=
                static_cast<int>(HaloLevel::Corner2D))
            groups.push_back(entry.first);
    }
    std::sort(groups.begin(), groups.end());

    std::ostringstream descriptor;
    descriptor << "face2-edge-execution-v1";
    for (const std::string &group_name : groups)
    {
        const HaloTripletRequest &triplet =
            face_2form_triplets_.at(group_name);
        descriptor << "|group=" << group_name
                   << "|level=" << static_cast<int>(triplet.level)
                   << "|stage=" << static_cast<int>(HaloLevel::Corner2D)
                   << "|nghost=" << triplet.nghost;
        const std::string fields[3] = {
            triplet.xi, triplet.eta, triplet.zeta};
        for (int axis = 0; axis < 3; ++axis)
        {
            const FieldDescriptor &field = fld_->descriptor(fields[axis]);
            descriptor << "|field=" << fields[axis]
                       << "|location=" << static_cast<int>(field.location)
                       << "|kind=" << static_cast<int>(field.value_kind)
                       << "|ncomp=" << field.ncomp
                       << "|field_nghost=" << field.nghost
                       << "|halo=" << field.sync.do_halo
                       << "|oriented=" << field.sync.orientation_aware;
        }
    }
    if (!edge1_all_ranks_same_string(descriptor.str(), size))
    {
        invalidate_prepared_face_2form_edge_execution_();
        face_2form_edge_execution_mode_ =
            Face2FormEdgeExecutionMode::Legacy;
        return false;
    }

    if (mode == Face2FormEdgeExecutionMode::Legacy)
    {
        invalidate_prepared_face_2form_edge_execution_();
        face_2form_edge_execution_mode_ =
            Face2FormEdgeExecutionMode::Legacy;
        return true;
    }

    const StaggerLocation expected_locations[3] = {
        StaggerLocation::FaceXi, StaggerLocation::FaceEt,
        StaggerLocation::FaceZe};
    bool local_valid =
        prepared_edge_face_2form_metadata_valid_() && !groups.empty();
    for (const std::string &group_name : groups)
    {
        const auto entry = prepared_edge_face_2form_plans_.find(group_name);
        if (entry == prepared_edge_face_2form_plans_.end() ||
            entry->second.registry_generation != registry_generation_ ||
            entry->second.pattern_generation != pattern_generation_)
        {
            local_valid = false;
            continue;
        }

        const LegacyFace2FormEdgePlanSnapshot &snapshot =
            entry->second.snapshot;
        const HaloTripletRequest &triplet =
            face_2form_triplets_.at(group_name);
        const std::array<std::string, 3> fields = {
            triplet.xi, triplet.eta, triplet.zeta};
        if (snapshot.group_name != group_name ||
            snapshot.field_names != fields || snapshot.ncomp <= 0 ||
            snapshot.nghost != triplet.nghost ||
            !triplet.orientation_aware ||
            static_cast<int>(triplet.level) <
                static_cast<int>(HaloLevel::Corner2D))
        {
            local_valid = false;
            continue;
        }
        for (int axis = 0; axis < 3; ++axis)
        {
            const FieldDescriptor &field = fld_->descriptor(fields[axis]);
            if (field.location != expected_locations[axis] ||
                field.value_kind !=
                    FieldValueKind::FaceContravariant2Form ||
                field.ncomp != snapshot.ncomp ||
                field.nghost != snapshot.nghost ||
                !field.sync.do_halo || !field.sync.orientation_aware)
                local_valid = false;
        }
        for (std::size_t index = 0; index < snapshot.recv_plan.size();
             ++index)
        {
            const auto &transfer = snapshot.recv_plan[index];
            if (transfer.local_plan_order != index ||
                transfer.peer_message_order !=
                    static_cast<std::size_t>(transfer.tag) ||
                !face2_valid_edge_transfer(
                    transfer, LegacyOrientedCornerRole::Receive, rank,
                    size, snapshot.ncomp, snapshot.nghost))
                local_valid = false;
        }
        for (std::size_t index = 0; index < snapshot.send_plan.size();
             ++index)
        {
            const auto &transfer = snapshot.send_plan[index];
            if (transfer.local_plan_order != index ||
                transfer.peer_message_order !=
                    static_cast<std::size_t>(transfer.tag) ||
                !face2_valid_edge_transfer(
                    transfer, LegacyOrientedCornerRole::Send, rank,
                    size, snapshot.ncomp, snapshot.nghost))
                local_valid = false;
        }
    }

    if (!edge1_all_ranks_true(local_valid))
    {
        invalidate_prepared_face_2form_edge_execution_();
        face_2form_edge_execution_mode_ =
            Face2FormEdgeExecutionMode::Legacy;
        return false;
    }

    invalidate_prepared_face_2form_edge_execution_();
    for (const std::string &group_name : groups)
    {
        const auto entry = prepared_edge_face_2form_plans_.find(group_name);
        prepared_face_2form_edge_execution_plans_.emplace(
            group_name,
            make_prepared_face_2form_edge_execution_view_(
                entry->second.snapshot));
    }
    reserve_prepared_exchange_views_();
    face_2form_edge_execution_registry_generation_ =
        registry_generation_;
    face_2form_edge_execution_pattern_generation_ =
        pattern_generation_;
    face_2form_edge_execution_mode_ =
        Face2FormEdgeExecutionMode::Prepared;
    return true;
}

bool Halo::activate_prepared_face_2form_edge_execution()
{
    return set_face_2form_edge_execution_mode(
        Face2FormEdgeExecutionMode::Prepared);
}

bool Halo::deactivate_prepared_face_2form_edge_execution()
{
    return set_face_2form_edge_execution_mode(
        Face2FormEdgeExecutionMode::Legacy);
}

Face2FormEdgeExecutionMode
Halo::debug_face_2form_edge_execution_mode() const noexcept
{
    return face_2form_edge_execution_mode_;
}

bool Halo::debug_prepared_face_2form_edge_execution_active() const noexcept
{
    return prepared_face_2form_edge_execution_valid_();
}

Face2FormEdgeExecutionCounters
Halo::debug_face_2form_edge_execution_counters() const noexcept
{
    return face_2form_edge_execution_counters_;
}


void Halo::enable_vertex_edge_1form_plan_capture(bool enabled)
{
    if (enabled)
        vertex_edge_1form_plan_snapshots_.clear();
    capture_vertex_edge_1form_plans_ = enabled;
}

std::optional<LegacyOrientedVertexPlanSnapshot>
Halo::debug_vertex_edge_1form_plan(
    const std::string &group_name) const
{
    const auto it =
        vertex_edge_1form_plan_snapshots_.find(group_name);
    if (it == vertex_edge_1form_plan_snapshots_.end())
        return std::nullopt;
    return it->second;
}

void Halo::enable_edge_face_2form_plan_capture(bool enabled)
{
    if (enabled)
        edge_face_2form_plan_snapshots_.clear();
    capture_edge_face_2form_plans_ = enabled;
}

std::optional<LegacyFace2FormEdgePlanSnapshot>
Halo::debug_edge_face_2form_plan(
    const std::string &group_name) const
{
    const auto it =
        edge_face_2form_plan_snapshots_.find(group_name);
    if (it == edge_face_2form_plan_snapshots_.end())
        return std::nullopt;
    return it->second;
}

void Halo::enable_vertex_face_2form_plan_capture(bool enabled)
{
    if (enabled)
        vertex_face_2form_plan_snapshots_.clear();
    capture_vertex_face_2form_plans_ = enabled;
}

std::optional<LegacyFace2FormVertexPlanSnapshot>
Halo::debug_vertex_face_2form_plan(
    const std::string &group_name) const
{
    const auto it =
        vertex_face_2form_plan_snapshots_.find(group_name);
    if (it == vertex_face_2form_plan_snapshots_.end())
        return std::nullopt;
    return it->second;
}

void Halo::data_trans_1DCorner(std::string &field_name)
{
    exchange_inner(field_name);
    exchange_parallel(field_name);
}

void Halo::data_trans_2DCorner(std::string &field_name)
{
    exchange_inner_edge(field_name);
    exchange_parallel_edge(field_name);
}

void Halo::data_trans_3DCorner(std::string &field_name)
{
    exchange_inner_vertex(field_name);
    exchange_parallel_vertex(field_name);
}

void Halo::coupling_trans_1DCorner(std::string &src, std::string &dst)
{
    coupling_inner_face(src, dst);
    coupling_parallel_face(src, dst);
}

void Halo::coupling_trans_1DCorner(std::string &src,
                                   std::string &dst,
                                   std::vector<int32_t> &field_cids)
{
    coupling_inner_face(src, dst, field_cids);
    coupling_parallel_face(src, dst, field_cids);
}

void Halo::coupling_trans_2DCorner(std::string &src, std::string &dst)
{
    coupling_inner_edge(src, dst);
    coupling_parallel_edge(src, dst);
}

void Halo::coupling_trans_2DCorner(std::string &src,
                                   std::string &dst,
                                   std::vector<int32_t> &field_cids)
{
    coupling_inner_edge(src, dst, field_cids);
    coupling_parallel_edge(src, dst, field_cids);
}

void Halo::coupling_trans_3DCorner(std::string &src, std::string &dst)
{
    coupling_inner_vertex(src, dst);
    coupling_parallel_vertex(src, dst);
}

void Halo::coupling_trans_3DCorner(std::string &src,
                                   std::string &dst,
                                   std::vector<int32_t> &field_cids)
{
    coupling_inner_vertex(src, dst, field_cids);
    coupling_parallel_vertex(src, dst, field_cids);
}
void Halo::invalidate_prepared_vertex_face_2form_metadata_()
{
    invalidate_prepared_face_2form_vertex_execution_();
    prepared_vertex_face_2form_plans_.clear();
    prepared_vertex_face_2form_metadata_ready_ = false;
    prepared_vertex_face_2form_registry_generation_ = 0;
    prepared_vertex_face_2form_pattern_generation_ = 0;
}

bool Halo::prepared_vertex_face_2form_metadata_valid_() const noexcept
{
    return prepared_vertex_face_2form_metadata_ready_ &&
           registered_patterns_built_ &&
           prepared_vertex_face_2form_registry_generation_ ==
               registry_generation_ &&
           prepared_vertex_face_2form_pattern_generation_ ==
               pattern_generation_;
}

std::optional<LegacyFace2FormVertexPlanSnapshot>
Halo::debug_prepared_vertex_face_2form_plan(
    const std::string &group_name) const
{
    if (!prepared_vertex_face_2form_metadata_valid_())
        return std::nullopt;

    const auto it = prepared_vertex_face_2form_plans_.find(group_name);
    if (it == prepared_vertex_face_2form_plans_.end() ||
        it->second.registry_generation != registry_generation_ ||
        it->second.pattern_generation != pattern_generation_)
        return std::nullopt;
    return it->second.snapshot;
}

std::optional<std::uint64_t>
Halo::debug_prepared_vertex_face_2form_plan_generation() const
{
    if (!prepared_vertex_face_2form_metadata_valid_())
        return std::nullopt;
    return prepared_vertex_face_2form_pattern_generation_;
}

std::uint64_t
Halo::debug_prepared_vertex_face_2form_metadata_exchange_count() const noexcept
{
    return prepared_vertex_face_2form_metadata_exchange_count_;
}
void Halo::invalidate_prepared_face_2form_vertex_execution_()
{
    prepared_face_2form_vertex_execution_plans_.clear();
    face_2form_vertex_execution_registry_generation_ = 0;
    face_2form_vertex_execution_pattern_generation_ = 0;
}

bool Halo::prepared_face_2form_vertex_execution_valid_() const noexcept
{
    return face_2form_vertex_execution_mode_ ==
               Face2FormVertexExecutionMode::Prepared &&
           prepared_vertex_face_2form_metadata_valid_() &&
           face_2form_vertex_execution_registry_generation_ ==
               registry_generation_ &&
           face_2form_vertex_execution_pattern_generation_ ==
               pattern_generation_ &&
           !prepared_face_2form_vertex_execution_plans_.empty();
}

bool Halo::set_face_2form_vertex_execution_mode(
    Face2FormVertexExecutionMode mode)
{
    int rank = 0;
    int size = 1;
    PARALLEL::mpi_rank(&rank);
    PARALLEL::mpi_size(&size);

    const int requested_mode = static_cast<int>(mode);
    std::vector<int> requested_modes(size, 0);
    MPI_Allgather(
        &requested_mode, 1, MPI_INT, requested_modes.data(), 1,
        MPI_INT, MPI_COMM_WORLD);
    const bool same_mode = std::all_of(
        requested_modes.begin(), requested_modes.end(),
        [&](int remote) { return remote == requested_mode; });
    if (!same_mode)
    {
        invalidate_prepared_face_2form_vertex_execution_();
        face_2form_vertex_execution_mode_ =
            Face2FormVertexExecutionMode::Legacy;
        return false;
    }

    std::vector<std::string> groups;
    for (const auto &entry : face_2form_triplets_)
    {
        if (entry.second.orientation_aware &&
            static_cast<int>(entry.second.level) >=
                static_cast<int>(HaloLevel::Corner3D))
            groups.push_back(entry.first);
    }
    std::sort(groups.begin(), groups.end());

    std::ostringstream descriptor;
    descriptor << "face2-vertex-execution-v1";
    for (const std::string &group_name : groups)
    {
        const HaloTripletRequest &triplet =
            face_2form_triplets_.at(group_name);
        descriptor << "|group=" << group_name
                   << "|level=" << static_cast<int>(triplet.level)
                   << "|stage=" << static_cast<int>(HaloLevel::Corner3D)
                   << "|nghost=" << triplet.nghost;
        const std::string fields[3] = {
            triplet.xi, triplet.eta, triplet.zeta};
        for (int axis = 0; axis < 3; ++axis)
        {
            const FieldDescriptor &field = fld_->descriptor(fields[axis]);
            descriptor << "|field=" << fields[axis]
                       << "|location=" << static_cast<int>(field.location)
                       << "|kind=" << static_cast<int>(field.value_kind)
                       << "|ncomp=" << field.ncomp
                       << "|field_nghost=" << field.nghost
                       << "|halo=" << field.sync.do_halo
                       << "|oriented=" << field.sync.orientation_aware;
        }
    }
    if (!edge1_all_ranks_same_string(descriptor.str(), size))
    {
        invalidate_prepared_face_2form_vertex_execution_();
        face_2form_vertex_execution_mode_ =
            Face2FormVertexExecutionMode::Legacy;
        return false;
    }

    if (mode == Face2FormVertexExecutionMode::Legacy)
    {
        invalidate_prepared_face_2form_vertex_execution_();
        face_2form_vertex_execution_mode_ =
            Face2FormVertexExecutionMode::Legacy;
        return true;
    }

    const StaggerLocation expected_locations[3] = {
        StaggerLocation::FaceXi, StaggerLocation::FaceEt,
        StaggerLocation::FaceZe};
    bool local_valid =
        prepared_vertex_face_2form_metadata_valid_() && !groups.empty();
    for (const std::string &group_name : groups)
    {
        const auto entry =
            prepared_vertex_face_2form_plans_.find(group_name);
        if (entry == prepared_vertex_face_2form_plans_.end() ||
            entry->second.registry_generation != registry_generation_ ||
            entry->second.pattern_generation != pattern_generation_)
        {
            local_valid = false;
            continue;
        }

        const LegacyFace2FormVertexPlanSnapshot &snapshot =
            entry->second.snapshot;
        const HaloTripletRequest &triplet =
            face_2form_triplets_.at(group_name);
        const std::array<std::string, 3> fields = {
            triplet.xi, triplet.eta, triplet.zeta};
        if (snapshot.group_name != group_name ||
            snapshot.field_names != fields || snapshot.ncomp <= 0 ||
            snapshot.nghost != triplet.nghost ||
            !triplet.orientation_aware ||
            static_cast<int>(triplet.level) <
                static_cast<int>(HaloLevel::Corner3D))
        {
            local_valid = false;
            continue;
        }
        for (int axis = 0; axis < 3; ++axis)
        {
            const FieldDescriptor &field = fld_->descriptor(fields[axis]);
            if (field.location != expected_locations[axis] ||
                field.value_kind !=
                    FieldValueKind::FaceContravariant2Form ||
                field.ncomp != snapshot.ncomp ||
                field.nghost != snapshot.nghost ||
                !field.sync.do_halo || !field.sync.orientation_aware)
                local_valid = false;
        }
        for (std::size_t index = 0; index < snapshot.recv_plan.size();
             ++index)
        {
            const auto &transfer = snapshot.recv_plan[index];
            if (transfer.local_plan_order != index ||
                transfer.peer_message_order !=
                    static_cast<std::size_t>(transfer.tag) ||
                !face2_valid_vertex_transfer(
                    transfer, LegacyOrientedCornerRole::Receive, rank,
                    size, snapshot.ncomp, snapshot.nghost))
                local_valid = false;
        }
        for (std::size_t index = 0; index < snapshot.send_plan.size();
             ++index)
        {
            const auto &transfer = snapshot.send_plan[index];
            if (transfer.local_plan_order != index ||
                transfer.peer_message_order !=
                    static_cast<std::size_t>(transfer.tag) ||
                !face2_valid_vertex_transfer(
                    transfer, LegacyOrientedCornerRole::Send, rank,
                    size, snapshot.ncomp, snapshot.nghost))
                local_valid = false;
        }
    }

    if (!edge1_all_ranks_true(local_valid))
    {
        invalidate_prepared_face_2form_vertex_execution_();
        face_2form_vertex_execution_mode_ =
            Face2FormVertexExecutionMode::Legacy;
        return false;
    }

    invalidate_prepared_face_2form_vertex_execution_();
    for (const std::string &group_name : groups)
    {
        const auto entry =
            prepared_vertex_face_2form_plans_.find(group_name);
        prepared_face_2form_vertex_execution_plans_.emplace(
            group_name,
            make_prepared_face_2form_vertex_execution_view_(
                entry->second.snapshot));
    }
    reserve_prepared_exchange_views_();
    face_2form_vertex_execution_registry_generation_ =
        registry_generation_;
    face_2form_vertex_execution_pattern_generation_ =
        pattern_generation_;
    face_2form_vertex_execution_mode_ =
        Face2FormVertexExecutionMode::Prepared;
    return true;
}

bool Halo::activate_prepared_face_2form_vertex_execution()
{
    return set_face_2form_vertex_execution_mode(
        Face2FormVertexExecutionMode::Prepared);
}

bool Halo::deactivate_prepared_face_2form_vertex_execution()
{
    return set_face_2form_vertex_execution_mode(
        Face2FormVertexExecutionMode::Legacy);
}

Face2FormVertexExecutionMode
Halo::debug_face_2form_vertex_execution_mode() const noexcept
{
    return face_2form_vertex_execution_mode_;
}

bool Halo::debug_prepared_face_2form_vertex_execution_active() const noexcept
{
    return prepared_face_2form_vertex_execution_valid_();
}

Face2FormVertexExecutionCounters
Halo::debug_face_2form_vertex_execution_counters() const noexcept
{
    return face_2form_vertex_execution_counters_;
}
