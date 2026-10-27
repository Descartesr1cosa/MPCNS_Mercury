#include "4_halo/Halo.h"

#include "0_basic/Error.h"
#include "0_basic/MPI_WRAPPER.h"
#include "2_topology/TopologyView.h"
#include "4_halo/detail/HaloBuildBoxMakers.h"
#include "4_halo/detail/HaloBuildTools.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <map>
#include <optional>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace
{
struct PreparedTripletCornerMeta
{
    int recv_rank;
    int send_rank;
    int recv_block;
    int send_block;

    int recv_axis;
    int send_axis;
    int sign;

    Box3 recv_box;
    Box3 node_box_on_send;
    int dir1_send;
    int dir2_send;
    int dir3_send;

    TOPO::IndexTransform trans_recv_to_send;
    int tag;
};

struct PreparedGroupSpec
{
    std::string group_name;
    std::array<std::string, 3> fields;
    HaloLevel level = HaloLevel::Corner1D;
    int nghost = 0;
    std::array<int, 3> ncomp{};
    std::array<StaggerLocation, 3> locations{};
    std::array<FieldValueKind, 3> value_kinds{};
    bool orientation_aware = false;
};

StaggerLocation edge_location(int axis)
{
    if (axis == 0)
        return StaggerLocation::EdgeXi;
    if (axis == 1)
        return StaggerLocation::EdgeEt;
    return StaggerLocation::EdgeZe;
}

StaggerLocation face_location(int axis)
{
    if (axis == 0)
        return StaggerLocation::FaceXi;
    if (axis == 1)
        return StaggerLocation::FaceEt;
    return StaggerLocation::FaceZe;
}

int inverse_perm_axis_for_prepared(
    const TOPO::IndexTransform &transform, int recv_axis)
{
    for (int send_axis = 0; send_axis < 3; ++send_axis)
        if (transform.perm[send_axis] == recv_axis)
            return send_axis;
    ERROR::Abort("[Halo] prepared Face2Form Edge metadata found an invalid permutation");
    return -1;
}
int box_volume(const Box3 &box)
{
    return (box.hi.i - box.lo.i) *
           (box.hi.j - box.lo.j) *
           (box.hi.k - box.lo.k);
}

void append_schedule_atom(
    std::ostringstream &stream, const std::string &value)
{
    stream << value.size() << ':' << value << ';';
}

std::string schedule_descriptor(
    const char *protocol,
    bool patterns_built,
    std::uint64_t registry_generation,
    std::uint64_t pattern_generation,
    bool cache_valid,
    const std::vector<PreparedGroupSpec> &groups)
{
    std::ostringstream stream;
    stream << protocol << "|patterns="
           << (patterns_built ? 1 : 0)
           << "|registry=" << registry_generation
           << "|pattern=" << pattern_generation
           << "|cache=" << (cache_valid ? 1 : 0)
           << "|groups=" << groups.size() << ';';
    for (const PreparedGroupSpec &group : groups)
    {
        append_schedule_atom(stream, group.group_name);
        stream << static_cast<int>(group.level) << ','
               << group.nghost << ','
               << (group.orientation_aware ? 1 : 0) << ';';
        for (int axis = 0; axis < 3; ++axis)
        {
            append_schedule_atom(stream, group.fields[axis]);
            stream << group.ncomp[axis] << ','
                   << static_cast<int>(group.locations[axis]) << ','
                   << static_cast<int>(group.value_kinds[axis])
                   << ';';
        }
    }
    return stream.str();
}

void verify_global_schedule(
    const std::string &schedule, int rank, int size)
{
    const int local_length = static_cast<int>(schedule.size());
    std::vector<int> lengths(size, 0);
    MPI_Allgather(
        &local_length, 1, MPI_INT, lengths.data(), 1, MPI_INT,
        MPI_COMM_WORLD);

    std::vector<int> displacements(size, 0);
    int total_length = 0;
    for (int peer = 0; peer < size; ++peer)
    {
        displacements[peer] = total_length;
        total_length += lengths[peer];
    }

    std::vector<char> gathered(
        static_cast<std::size_t>(std::max(1, total_length)));
    MPI_Allgatherv(
        const_cast<char *>(schedule.data()), local_length, MPI_CHAR,
        gathered.data(), lengths.data(), displacements.data(), MPI_CHAR,
        MPI_COMM_WORLD);

    for (int peer = 0; peer < size; ++peer)
    {
        const std::string peer_schedule(
            gathered.data() + displacements[peer],
            static_cast<std::size_t>(lengths[peer]));
        if (peer_schedule == schedule)
            continue;

        std::ostringstream message;
        message << "[Halo][rank " << rank
                << "] prepared Edge1Form Edge-stage registry schedule "
                   "mismatch with rank "
                << peer
                << " before metadata collectives";
        ERROR::Abort(message.str());
    }
}

void exchange_prepared_triplet_corner_metadata(
    const std::map<int, std::vector<PreparedTripletCornerMeta>>
        &metadata_to_send,
    std::vector<PreparedTripletCornerMeta> &received_metadata)
{
    int size = 1;
    PARALLEL::mpi_size(&size);

    std::vector<int> send_counts(size, 0);
    std::vector<int> recv_counts(size, 0);
    for (const auto &entry : metadata_to_send)
        send_counts[entry.first] =
            static_cast<int>(entry.second.size());

    PARALLEL::mpi_alltoall(
        send_counts.data(), 1, recv_counts.data(), 1);

    std::vector<int> send_displacements(size, 0);
    std::vector<int> recv_displacements(size, 0);
    int total_send = 0;
    int total_recv = 0;
    for (int peer = 0; peer < size; ++peer)
    {
        send_displacements[peer] = total_send;
        recv_displacements[peer] = total_recv;
        total_send += send_counts[peer];
        total_recv += recv_counts[peer];
    }

    std::vector<PreparedTripletCornerMeta> send_buffer(total_send);
    received_metadata.resize(total_recv);
    for (int peer = 0; peer < size; ++peer)
    {
        const auto it = metadata_to_send.find(peer);
        if (it == metadata_to_send.end())
            continue;
        std::copy(
            it->second.begin(), it->second.end(),
            send_buffer.begin() + send_displacements[peer]);
    }

    const int bytes =
        static_cast<int>(sizeof(PreparedTripletCornerMeta));
    std::vector<int> send_counts_bytes(size);
    std::vector<int> recv_counts_bytes(size);
    std::vector<int> send_displacements_bytes(size);
    std::vector<int> recv_displacements_bytes(size);
    for (int peer = 0; peer < size; ++peer)
    {
        send_counts_bytes[peer] = send_counts[peer] * bytes;
        recv_counts_bytes[peer] = recv_counts[peer] * bytes;
        send_displacements_bytes[peer] =
            send_displacements[peer] * bytes;
        recv_displacements_bytes[peer] =
            recv_displacements[peer] * bytes;
    }

    MPI_Alltoallv(
        reinterpret_cast<const void *>(send_buffer.data()),
        send_counts_bytes.data(),
        send_displacements_bytes.data(),
        MPI_BYTE,
        reinterpret_cast<void *>(received_metadata.data()),
        recv_counts_bytes.data(),
        recv_displacements_bytes.data(),
        MPI_BYTE,
        MPI_COMM_WORLD);
}

LegacyOrientedCornerTransferSnapshot make_snapshot(
    const PreparedTripletCornerMeta &metadata,
    LegacyOrientedCornerRole role,
    const std::optional<Box3> &send_box,
    int ncomp,
    int nghost,
    std::size_t local_order)
{
    LegacyOrientedCornerTransferSnapshot snapshot;
    snapshot.role = role;
    snapshot.recv_rank = metadata.recv_rank;
    snapshot.send_rank = metadata.send_rank;
    snapshot.recv_block = metadata.recv_block;
    snapshot.send_block = metadata.send_block;
    snapshot.recv_axis = metadata.recv_axis;
    snapshot.send_axis = metadata.send_axis;
    snapshot.sign = metadata.sign;
    snapshot.recv_box = metadata.recv_box;
    snapshot.send_box = send_box;
    snapshot.node_box_on_send = metadata.node_box_on_send;
    snapshot.dir1_send = metadata.dir1_send;
    snapshot.dir2_send = metadata.dir2_send;
    snapshot.trans_recv_to_send = metadata.trans_recv_to_send;
    snapshot.peer_rank =
        role == LegacyOrientedCornerRole::Receive
            ? metadata.send_rank
            : metadata.recv_rank;
    snapshot.tag = metadata.tag;
    snapshot.ncomp = ncomp;
    snapshot.nghost = nghost;
    snapshot.message_length =
        box_volume(metadata.recv_box) * ncomp;
    snapshot.local_plan_order = local_order;
    snapshot.peer_message_order =
        static_cast<std::size_t>(metadata.tag);
    return snapshot;
}

LegacyFace2FormEdgeTransferSnapshot make_face2_edge_snapshot(
    const PreparedTripletCornerMeta &metadata,
    LegacyOrientedCornerRole role,
    const std::optional<Box3> &send_box,
    int ncomp,
    int nghost,
    std::size_t local_order)
{
    LegacyFace2FormEdgeTransferSnapshot snapshot;
    snapshot.role = role;
    snapshot.recv_rank = metadata.recv_rank;
    snapshot.send_rank = metadata.send_rank;
    snapshot.recv_block = metadata.recv_block;
    snapshot.send_block = metadata.send_block;
    snapshot.recv_axis = metadata.recv_axis;
    snapshot.send_axis = metadata.send_axis;
    snapshot.sign = metadata.sign;
    snapshot.recv_box = metadata.recv_box;
    snapshot.send_box = send_box;
    snapshot.node_box_on_send = metadata.node_box_on_send;
    snapshot.dir1_send = metadata.dir1_send;
    snapshot.dir2_send = metadata.dir2_send;
    snapshot.trans_recv_to_send = metadata.trans_recv_to_send;
    snapshot.peer_rank =
        role == LegacyOrientedCornerRole::Receive
            ? metadata.send_rank
            : metadata.recv_rank;
    snapshot.tag = metadata.tag;
    snapshot.ncomp = ncomp;
    snapshot.nghost = nghost;
    snapshot.message_length =
        box_volume(metadata.recv_box) * ncomp;
    snapshot.local_plan_order = local_order;
    snapshot.peer_message_order =
        static_cast<std::size_t>(metadata.tag);
    return snapshot;
}
LegacyOrientedVertexTransferSnapshot make_vertex_snapshot(
    const PreparedTripletCornerMeta &metadata,
    LegacyOrientedCornerRole role,
    const std::optional<Box3> &send_box,
    int ncomp,
    int nghost,
    std::size_t local_order)
{
    LegacyOrientedVertexTransferSnapshot snapshot;
    snapshot.role = role;
    snapshot.recv_rank = metadata.recv_rank;
    snapshot.send_rank = metadata.send_rank;
    snapshot.recv_block = metadata.recv_block;
    snapshot.send_block = metadata.send_block;
    snapshot.recv_axis = metadata.recv_axis;
    snapshot.send_axis = metadata.send_axis;
    snapshot.sign = metadata.sign;
    snapshot.recv_box = metadata.recv_box;
    snapshot.send_box = send_box;
    snapshot.node_box_on_send = metadata.node_box_on_send;
    snapshot.dir1_send = metadata.dir1_send;
    snapshot.dir2_send = metadata.dir2_send;
    snapshot.dir3_send = metadata.dir3_send;
    snapshot.trans_recv_to_send = metadata.trans_recv_to_send;
    snapshot.peer_rank =
        role == LegacyOrientedCornerRole::Receive
            ? metadata.send_rank
            : metadata.recv_rank;
    snapshot.tag = metadata.tag;
    snapshot.ncomp = ncomp;
    snapshot.nghost = nghost;
    snapshot.message_length =
        box_volume(metadata.recv_box) * ncomp;
    snapshot.local_plan_order = local_order;
    snapshot.peer_message_order =
        static_cast<std::size_t>(metadata.tag);
    return snapshot;
}

bool valid_edge_group(const PreparedGroupSpec &group)
{
    if (!group.orientation_aware ||
        static_cast<int>(group.level) <
            static_cast<int>(HaloLevel::Corner2D) ||
        group.nghost < 0)
        return false;

    const StaggerLocation expected_locations[3] = {
        StaggerLocation::EdgeXi,
        StaggerLocation::EdgeEt,
        StaggerLocation::EdgeZe};
    for (int axis = 0; axis < 3; ++axis)
    {
        if (group.fields[axis].empty() ||
            group.ncomp[axis] <= 0 ||
            group.locations[axis] != expected_locations[axis] ||
            group.value_kinds[axis] !=
                FieldValueKind::EdgeCovariant1Form)
            return false;
    }
    return group.ncomp[0] == group.ncomp[1] &&
           group.ncomp[0] == group.ncomp[2];
}

bool valid_face2_group(const PreparedGroupSpec &group)
{
    if (!group.orientation_aware ||
        static_cast<int>(group.level) <
            static_cast<int>(HaloLevel::Corner2D) ||
        group.nghost < 0)
        return false;

    const StaggerLocation expected_locations[3] = {
        StaggerLocation::FaceXi,
        StaggerLocation::FaceEt,
        StaggerLocation::FaceZe};
    for (int axis = 0; axis < 3; ++axis)
    {
        if (group.fields[axis].empty() ||
            group.ncomp[axis] <= 0 ||
            group.locations[axis] != expected_locations[axis] ||
            group.value_kinds[axis] !=
                FieldValueKind::FaceContravariant2Form)
            return false;
    }
    return group.ncomp[0] == group.ncomp[1] &&
           group.ncomp[0] == group.ncomp[2];
}
} // namespace
bool Halo::prepare_edge_1form_corner_metadata()
{
    int rank = 0;
    int size = 1;
    PARALLEL::mpi_rank(&rank);
    PARALLEL::mpi_size(&size);

    std::vector<PreparedGroupSpec> groups;
    groups.reserve(edge_1form_triplets_.size());
    for (const auto &entry : edge_1form_triplets_)
    {
        const HaloTripletRequest &triplet = entry.second;
        if (static_cast<int>(triplet.level) <
            static_cast<int>(HaloLevel::Corner2D))
            continue;

        PreparedGroupSpec group;
        group.group_name = triplet.group_name;
        group.fields = {triplet.xi, triplet.eta, triplet.zeta};
        group.level = triplet.level;
        group.nghost = triplet.nghost;
        group.orientation_aware = triplet.orientation_aware;
        bool descriptors_are_orientation_aware =
            triplet.orientation_aware;
        for (int axis = 0; axis < 3; ++axis)
        {
            const FieldDescriptor &descriptor =
                fld_->descriptor(group.fields[axis]);
            group.ncomp[axis] = descriptor.ncomp;
            group.locations[axis] = descriptor.location;
            group.value_kinds[axis] = descriptor.value_kind;
            descriptors_are_orientation_aware =
                descriptors_are_orientation_aware &&
                descriptor.sync.do_halo &&
                descriptor.sync.orientation_aware;
        }
        // The active Field catalog is the final eligibility check. This
        // preserves component-copy behavior when form orientation is off,
        // even if a legacy triplet cache predates a registry rebuild.
        if (!descriptors_are_orientation_aware)
            continue;
        groups.push_back(std::move(group));
    }
    std::sort(
        groups.begin(), groups.end(),
        [](const PreparedGroupSpec &left,
           const PreparedGroupSpec &right)
        {
            return left.group_name < right.group_name;
        });

    const bool cache_valid =
        prepared_edge_1form_corner_metadata_valid_();
    verify_global_schedule(
        schedule_descriptor(
            "prepared-edge1-edge-v1",
            registered_patterns_built_, registry_generation_,
            pattern_generation_, cache_valid, groups),
        rank, size);

    if (!registered_patterns_built_)
        return false;
    if (cache_valid)
        return true;

    for (const PreparedGroupSpec &group : groups)
    {
        if (valid_edge_group(group))
            continue;

        std::ostringstream message;
        message << "[Halo][rank " << rank
                << "][group " << group.group_name
                << "][stage Edge] cannot prepare an invalid "
                   "Edge1Form triplet";
        ERROR::Abort(message.str());
    }

    std::unordered_map<std::string,
                       PreparedEdge1FormCornerPlanEntry>
        prepared;
    const auto &parallel_edges =
        TOPO_VIEW::edge_patches(*topo_, TOPO::PatchKind::Parallel);
    for (const PreparedGroupSpec &group : groups)
    {
        std::map<int, std::vector<PreparedTripletCornerMeta>>
            metadata_to_send;
        std::map<int, int> next_tag;
        std::vector<PreparedTripletCornerMeta> receive_plan;

        for (const TOPO::EdgePatch &patch : parallel_edges)
        {
            if (patch.is_coupling)
                continue;

            const Direction this_direction_1 =
                HALO_TOOLS::int_to_direction(patch.dir1);
            const Direction this_direction_2 =
                HALO_TOOLS::int_to_direction(patch.dir2);
            int send_direction_1 =
                HALO_TOOLS::map_dir_to_neighbor(
                    patch.dir1, patch.trans, true);
            int send_direction_2 =
                HALO_TOOLS::map_dir_to_neighbor(
                    patch.dir2, patch.trans, false);
            const int expected_normal_axis =
                patch.trans.perm[std::abs(patch.dir1) - 1] + 1;
            if (expected_normal_axis !=
                std::abs(send_direction_1))
            {
                if (expected_normal_axis ==
                    std::abs(send_direction_2))
                {
                    std::swap(
                        send_direction_1, send_direction_2);
                }
                else
                {
                    std::ostringstream message;
                    message << "[Halo][rank " << rank
                            << "][group " << group.group_name
                            << "][stage Edge] legacy direction "
                               "mapping is invalid";
                    ERROR::Abort(message.str());
                }
            }

            for (int recv_axis = 0; recv_axis < 3; ++recv_axis)
            {
                PreparedTripletCornerMeta metadata{};
                metadata.recv_rank = patch.this_rank;
                metadata.send_rank = patch.nb_rank;
                metadata.recv_block = patch.this_block;
                metadata.send_block = patch.nb_block;
                metadata.recv_axis = recv_axis;
                metadata.send_axis = patch.trans.perm[recv_axis];
                metadata.sign = patch.trans.sign[recv_axis];
                metadata.recv_box =
                    HALO_BOX::make_2DCorner_ghost_box(
                        edge_location(recv_axis),
                        patch.this_box_node, this_direction_1,
                        this_direction_2, group.nghost);
                metadata.node_box_on_send = patch.nb_box_node;
                metadata.dir1_send = send_direction_1;
                metadata.dir2_send = send_direction_2;
                metadata.dir3_send = 0;
                metadata.trans_recv_to_send = patch.trans;
                metadata.tag = next_tag[patch.nb_rank]++;

                receive_plan.push_back(metadata);
                metadata_to_send[patch.nb_rank].push_back(
                    metadata);
            }
        }

        std::vector<PreparedTripletCornerMeta> send_plan;
        // This is the legacy sender-discovery protocol. All ranks reach it
        // for every globally verified group, even with no local patches.
        exchange_prepared_triplet_corner_metadata(
            metadata_to_send, send_plan);

        LegacyOrientedCornerPlanSnapshot snapshot;
        snapshot.group_name = group.group_name;
        snapshot.field_names = group.fields;
        snapshot.ncomp = group.ncomp[0];
        snapshot.nghost = group.nghost;
        snapshot.recv_plan.reserve(receive_plan.size());
        snapshot.send_plan.reserve(send_plan.size());

        for (std::size_t order = 0;
             order < send_plan.size(); ++order)
        {
            const PreparedTripletCornerMeta &metadata =
                send_plan[order];
            if (metadata.send_rank != rank)
            {
                std::ostringstream message;
                message << "[Halo][rank " << rank
                        << "][group " << group.group_name
                        << "][stage Edge][role send] metadata "
                           "was delivered to the wrong rank";
                ERROR::Abort(message.str());
            }

            const Direction direction_1 =
                HALO_TOOLS::int_to_direction(
                    metadata.dir1_send);
            const Direction direction_2 =
                HALO_TOOLS::int_to_direction(
                    metadata.dir2_send);
            const Box3 send_box =
                HALO_BOX::make_2DCorner_innerghost_box(
                    edge_location(metadata.send_axis),
                    metadata.node_box_on_send,
                    direction_1, direction_2, group.nghost);
            snapshot.send_plan.push_back(make_snapshot(
                metadata, LegacyOrientedCornerRole::Send,
                std::optional<Box3>{send_box},
                group.ncomp[0], group.nghost, order));
        }

        for (std::size_t order = 0;
             order < receive_plan.size(); ++order)
        {
            snapshot.recv_plan.push_back(make_snapshot(
                receive_plan[order],
                LegacyOrientedCornerRole::Receive,
                std::nullopt, group.ncomp[0],
                group.nghost, order));
        }

        PreparedEdge1FormCornerPlanEntry entry;
        entry.snapshot = std::move(snapshot);
        entry.registry_generation = registry_generation_;
        entry.pattern_generation = pattern_generation_;
        prepared.emplace(group.group_name, std::move(entry));
    }

    prepared_edge_1form_corner_plans_ = std::move(prepared);
    prepared_edge_1form_corner_metadata_ready_ = true;
    prepared_edge_1form_corner_registry_generation_ =
        registry_generation_;
    prepared_edge_1form_corner_pattern_generation_ =
        pattern_generation_;
    prepared_edge_1form_corner_metadata_exchange_count_ +=
        static_cast<std::uint64_t>(groups.size());
    return true;
}

bool Halo::prepare_vertex_edge_1form_metadata()
{
    int rank = 0;
    int size = 1;
    PARALLEL::mpi_rank(&rank);
    PARALLEL::mpi_size(&size);

    std::vector<PreparedGroupSpec> groups;
    groups.reserve(edge_1form_triplets_.size());
    for (const auto &entry : edge_1form_triplets_)
    {
        const HaloTripletRequest &triplet = entry.second;
        if (static_cast<int>(triplet.level) <
            static_cast<int>(HaloLevel::Corner3D))
            continue;

        PreparedGroupSpec group;
        group.group_name = triplet.group_name;
        group.fields = {triplet.xi, triplet.eta, triplet.zeta};
        group.level = triplet.level;
        group.nghost = triplet.nghost;
        group.orientation_aware = triplet.orientation_aware;
        bool descriptors_are_orientation_aware =
            triplet.orientation_aware;
        for (int axis = 0; axis < 3; ++axis)
        {
            const FieldDescriptor &descriptor =
                fld_->descriptor(group.fields[axis]);
            group.ncomp[axis] = descriptor.ncomp;
            group.locations[axis] = descriptor.location;
            group.value_kinds[axis] = descriptor.value_kind;
            descriptors_are_orientation_aware =
                descriptors_are_orientation_aware &&
                descriptor.sync.do_halo &&
                descriptor.sync.orientation_aware;
        }
        if (!descriptors_are_orientation_aware)
            continue;
        groups.push_back(std::move(group));
    }
    std::sort(
        groups.begin(), groups.end(),
        [](const PreparedGroupSpec &left,
           const PreparedGroupSpec &right)
        {
            return left.group_name < right.group_name;
        });

    const bool cache_valid =
        prepared_vertex_edge_1form_metadata_valid_();
    verify_global_schedule(
        schedule_descriptor(
            "prepared-edge1-vertex-v1",
            registered_patterns_built_, registry_generation_,
            pattern_generation_, cache_valid, groups),
        rank, size);

    if (!registered_patterns_built_)
        return false;
    if (cache_valid)
        return true;

    for (const PreparedGroupSpec &group : groups)
    {
        if (valid_edge_group(group))
            continue;

        std::ostringstream message;
        message << "[Halo][rank " << rank
                << "][group " << group.group_name
                << "][stage Vertex] cannot prepare an invalid "
                   "Edge1Form triplet";
        ERROR::Abort(message.str());
    }

    std::unordered_map<std::string,
                       PreparedVertexEdge1FormPlanEntry>
        prepared;
    const auto &parallel_vertices =
        TOPO_VIEW::vertex_patches(*topo_, TOPO::PatchKind::Parallel);
    for (const PreparedGroupSpec &group : groups)
    {
        std::map<int, std::vector<PreparedTripletCornerMeta>>
            metadata_to_send;
        std::map<int, int> next_tag;
        std::vector<PreparedTripletCornerMeta> receive_plan;

        for (const TOPO::VertexPatch &patch : parallel_vertices)
        {
            if (patch.is_coupling)
                continue;

            const Direction this_direction_1 =
                HALO_TOOLS::int_to_direction(patch.dir1);
            const Direction this_direction_2 =
                HALO_TOOLS::int_to_direction(patch.dir2);
            const Direction this_direction_3 =
                HALO_TOOLS::int_to_direction(patch.dir3);
            int send_direction_1 =
                HALO_TOOLS::map_dir_to_neighbor(
                    patch.dir1, patch.trans, true);
            int send_direction_2 =
                HALO_TOOLS::map_dir_to_neighbor(
                    patch.dir2, patch.trans, false);
            int send_direction_3 =
                HALO_TOOLS::map_dir_to_neighbor(
                    patch.dir3, patch.trans, false);
            const int expected_normal_axis =
                patch.trans.perm[std::abs(patch.dir1) - 1] + 1;
            if (expected_normal_axis !=
                std::abs(send_direction_1))
            {
                if (expected_normal_axis ==
                    std::abs(send_direction_2))
                {
                    std::swap(
                        send_direction_1, send_direction_2);
                }
                else if (expected_normal_axis ==
                         std::abs(send_direction_3))
                {
                    std::swap(
                        send_direction_1, send_direction_3);
                }
                else
                {
                    std::ostringstream message;
                    message << "[Halo][rank " << rank
                            << "][group " << group.group_name
                            << "][stage Vertex] legacy direction "
                               "mapping is invalid";
                    ERROR::Abort(message.str());
                }
            }

            for (int recv_axis = 0; recv_axis < 3; ++recv_axis)
            {
                PreparedTripletCornerMeta metadata{};
                metadata.recv_rank = patch.this_rank;
                metadata.send_rank = patch.nb_rank;
                metadata.recv_block = patch.this_block;
                metadata.send_block = patch.nb_block;
                metadata.recv_axis = recv_axis;
                metadata.send_axis = patch.trans.perm[recv_axis];
                metadata.sign = patch.trans.sign[recv_axis];
                metadata.recv_box =
                    HALO_BOX::make_3DCorner_ghost_box(
                        edge_location(recv_axis),
                        patch.this_box_node, this_direction_1,
                        this_direction_2, this_direction_3,
                        group.nghost);
                metadata.node_box_on_send = patch.nb_box_node;
                metadata.dir1_send = send_direction_1;
                metadata.dir2_send = send_direction_2;
                metadata.dir3_send = send_direction_3;
                metadata.trans_recv_to_send = patch.trans;
                metadata.tag = next_tag[patch.nb_rank]++;

                receive_plan.push_back(metadata);
                metadata_to_send[patch.nb_rank].push_back(
                    metadata);
            }
        }

        std::vector<PreparedTripletCornerMeta> send_plan;
        // The legacy sender-discovery protocol remains rank-grouped and all
        // ranks participate for every globally verified group.
        exchange_prepared_triplet_corner_metadata(
            metadata_to_send, send_plan);

        LegacyOrientedVertexPlanSnapshot snapshot;
        snapshot.group_name = group.group_name;
        snapshot.field_names = group.fields;
        snapshot.ncomp = group.ncomp[0];
        snapshot.nghost = group.nghost;
        snapshot.recv_plan.reserve(receive_plan.size());
        snapshot.send_plan.reserve(send_plan.size());

        for (std::size_t order = 0;
             order < send_plan.size(); ++order)
        {
            const PreparedTripletCornerMeta &metadata =
                send_plan[order];
            if (metadata.send_rank != rank)
            {
                std::ostringstream message;
                message << "[Halo][rank " << rank
                        << "][group " << group.group_name
                        << "][stage Vertex][role send] metadata "
                           "was delivered to the wrong rank";
                ERROR::Abort(message.str());
            }

            const Direction direction_1 =
                HALO_TOOLS::int_to_direction(
                    metadata.dir1_send);
            const Direction direction_2 =
                HALO_TOOLS::int_to_direction(
                    metadata.dir2_send);
            const Direction direction_3 =
                HALO_TOOLS::int_to_direction(
                    metadata.dir3_send);
            const Box3 send_box =
                HALO_BOX::make_3DCorner_innerghost_box(
                    edge_location(metadata.send_axis),
                    metadata.node_box_on_send,
                    direction_1, direction_2, direction_3,
                    group.nghost);
            snapshot.send_plan.push_back(make_vertex_snapshot(
                metadata, LegacyOrientedCornerRole::Send,
                std::optional<Box3>{send_box}, group.ncomp[0],
                group.nghost, order));
        }

        for (std::size_t order = 0;
             order < receive_plan.size(); ++order)
        {
            snapshot.recv_plan.push_back(make_vertex_snapshot(
                receive_plan[order],
                LegacyOrientedCornerRole::Receive,
                std::nullopt, group.ncomp[0], group.nghost,
                order));
        }

        PreparedVertexEdge1FormPlanEntry entry;
        entry.snapshot = std::move(snapshot);
        entry.registry_generation = registry_generation_;
        entry.pattern_generation = pattern_generation_;
        prepared.emplace(group.group_name, std::move(entry));
    }

    prepared_vertex_edge_1form_plans_ = std::move(prepared);
    prepared_vertex_edge_1form_metadata_ready_ = true;
    prepared_vertex_edge_1form_registry_generation_ =
        registry_generation_;
    prepared_vertex_edge_1form_pattern_generation_ =
        pattern_generation_;
    prepared_vertex_edge_1form_metadata_exchange_count_ +=
        static_cast<std::uint64_t>(groups.size());
    return true;
}

bool Halo::prepare_edge_face_2form_metadata()
{
    int rank = 0;
    int size = 1;
    PARALLEL::mpi_rank(&rank);
    PARALLEL::mpi_size(&size);

    std::vector<PreparedGroupSpec> groups;
    groups.reserve(face_2form_triplets_.size());
    for (const auto &entry : face_2form_triplets_)
    {
        const HaloTripletRequest &triplet = entry.second;
        if (static_cast<int>(triplet.level) <
            static_cast<int>(HaloLevel::Corner2D))
            continue;

        PreparedGroupSpec group;
        group.group_name = triplet.group_name;
        group.fields = {triplet.xi, triplet.eta, triplet.zeta};
        group.level = triplet.level;
        group.nghost = triplet.nghost;
        group.orientation_aware = triplet.orientation_aware;
        bool descriptors_are_orientation_aware =
            triplet.orientation_aware;
        for (int axis = 0; axis < 3; ++axis)
        {
            const FieldDescriptor &descriptor =
                fld_->descriptor(group.fields[axis]);
            group.ncomp[axis] = descriptor.ncomp;
            group.locations[axis] = descriptor.location;
            group.value_kinds[axis] = descriptor.value_kind;
            descriptors_are_orientation_aware =
                descriptors_are_orientation_aware &&
                descriptor.sync.do_halo &&
                descriptor.sync.orientation_aware;
        }
        if (!descriptors_are_orientation_aware)
            continue;
        groups.push_back(std::move(group));
    }
    std::sort(
        groups.begin(), groups.end(),
        [](const PreparedGroupSpec &left,
           const PreparedGroupSpec &right)
        {
            return left.group_name < right.group_name;
        });

    const bool cache_valid =
        prepared_edge_face_2form_metadata_valid_();
    verify_global_schedule(
        schedule_descriptor(
            "prepared-face2-edge-v1",
            registered_patterns_built_, registry_generation_,
            pattern_generation_, cache_valid, groups),
        rank, size);

    if (!registered_patterns_built_)
        return false;
    if (cache_valid)
        return true;

    for (const PreparedGroupSpec &group : groups)
    {
        if (valid_face2_group(group))
            continue;

        std::ostringstream message;
        message << "[Halo][rank " << rank
                << "][group " << group.group_name
                << "][stage Edge] cannot prepare an invalid "
                   "Face2Form triplet";
        ERROR::Abort(message.str());
    }

    std::unordered_map<std::string, PreparedEdgeFace2FormPlanEntry>
        prepared;
    const auto &parallel_edges =
        TOPO_VIEW::edge_patches(*topo_, TOPO::PatchKind::Parallel);
    for (const PreparedGroupSpec &group : groups)
    {
        std::map<int, std::vector<PreparedTripletCornerMeta>>
            metadata_to_send;
        std::map<int, int> next_tag;
        std::vector<PreparedTripletCornerMeta> receive_plan;

        for (const TOPO::EdgePatch &patch : parallel_edges)
        {
            if (patch.is_coupling)
                continue;

            const Direction this_direction_1 =
                HALO_TOOLS::int_to_direction(patch.dir1);
            const Direction this_direction_2 =
                HALO_TOOLS::int_to_direction(patch.dir2);
            int send_direction_1 =
                HALO_TOOLS::map_dir_to_neighbor(
                    patch.dir1, patch.trans, true);
            int send_direction_2 =
                HALO_TOOLS::map_dir_to_neighbor(
                    patch.dir2, patch.trans, false);
            const int expected_normal_axis =
                patch.trans.perm[std::abs(patch.dir1) - 1] + 1;
            if (expected_normal_axis !=
                std::abs(send_direction_1))
            {
                if (expected_normal_axis ==
                    std::abs(send_direction_2))
                {
                    std::swap(send_direction_1, send_direction_2);
                }
                else
                {
                    std::ostringstream message;
                    message << "[Halo][rank " << rank
                            << "][group " << group.group_name
                            << "][stage Edge] legacy direction "
                               "mapping is invalid";
                    ERROR::Abort(message.str());
                }
            }

            for (int recv_axis = 0; recv_axis < 3; ++recv_axis)
            {
                PreparedTripletCornerMeta metadata{};
                metadata.recv_rank = patch.this_rank;
                metadata.send_rank = patch.nb_rank;
                metadata.recv_block = patch.this_block;
                metadata.send_block = patch.nb_block;
                metadata.recv_axis = recv_axis;
                metadata.send_axis = inverse_perm_axis_for_prepared(
                    patch.trans, recv_axis);
                metadata.sign = HALO_TOOLS::face_2form_orientation_sign(
                    patch.trans, metadata.send_axis);
                metadata.recv_box =
                    HALO_BOX::make_2DCorner_ghost_box(
                        face_location(recv_axis), patch.this_box_node,
                        this_direction_1, this_direction_2,
                        group.nghost);
                metadata.node_box_on_send = patch.nb_box_node;
                metadata.dir1_send = send_direction_1;
                metadata.dir2_send = send_direction_2;
                metadata.dir3_send = 0;
                metadata.trans_recv_to_send = patch.trans;
                metadata.tag = next_tag[patch.nb_rank]++;

                receive_plan.push_back(metadata);
                metadata_to_send[patch.nb_rank].push_back(metadata);
            }
        }

        std::vector<PreparedTripletCornerMeta> send_plan;
        // All ranks retain the legacy rank-grouped sender discovery order.
        exchange_prepared_triplet_corner_metadata(
            metadata_to_send, send_plan);

        LegacyFace2FormEdgePlanSnapshot snapshot;
        snapshot.group_name = group.group_name;
        snapshot.field_names = group.fields;
        snapshot.ncomp = group.ncomp[0];
        snapshot.nghost = group.nghost;
        snapshot.recv_plan.reserve(receive_plan.size());
        snapshot.send_plan.reserve(send_plan.size());

        for (std::size_t order = 0;
             order < send_plan.size(); ++order)
        {
            const PreparedTripletCornerMeta &metadata =
                send_plan[order];
            if (metadata.send_rank != rank)
            {
                std::ostringstream message;
                message << "[Halo][rank " << rank
                        << "][group " << group.group_name
                        << "][stage Edge][role send] metadata "
                           "was delivered to the wrong rank";
                ERROR::Abort(message.str());
            }

            const Direction direction_1 =
                HALO_TOOLS::int_to_direction(metadata.dir1_send);
            const Direction direction_2 =
                HALO_TOOLS::int_to_direction(metadata.dir2_send);
            const Box3 send_box =
                HALO_BOX::make_2DCorner_innerghost_box(
                    face_location(metadata.send_axis),
                    metadata.node_box_on_send,
                    direction_1, direction_2, group.nghost);
            snapshot.send_plan.push_back(make_face2_edge_snapshot(
                metadata, LegacyOrientedCornerRole::Send,
                std::optional<Box3>{send_box}, group.ncomp[0],
                group.nghost, order));
        }

        for (std::size_t order = 0;
             order < receive_plan.size(); ++order)
        {
            snapshot.recv_plan.push_back(make_face2_edge_snapshot(
                receive_plan[order],
                LegacyOrientedCornerRole::Receive,
                std::nullopt, group.ncomp[0], group.nghost,
                order));
        }

        PreparedEdgeFace2FormPlanEntry entry;
        entry.snapshot = std::move(snapshot);
        entry.registry_generation = registry_generation_;
        entry.pattern_generation = pattern_generation_;
        prepared.emplace(group.group_name, std::move(entry));
    }

    prepared_edge_face_2form_plans_ = std::move(prepared);
    prepared_edge_face_2form_metadata_ready_ = true;
    prepared_edge_face_2form_registry_generation_ =
        registry_generation_;
    prepared_edge_face_2form_pattern_generation_ =
        pattern_generation_;
    prepared_edge_face_2form_metadata_exchange_count_ +=
        static_cast<std::uint64_t>(groups.size());
    return true;
}
namespace
{
LegacyFace2FormVertexTransferSnapshot
make_face2_vertex_snapshot(
    const PreparedTripletCornerMeta &metadata,
    LegacyOrientedCornerRole role,
    const std::optional<Box3> &send_box,
    int ncomp,
    int nghost,
    std::size_t local_order)
{
    LegacyFace2FormVertexTransferSnapshot snapshot;
    snapshot.role = role;
    snapshot.recv_rank = metadata.recv_rank;
    snapshot.send_rank = metadata.send_rank;
    snapshot.recv_block = metadata.recv_block;
    snapshot.send_block = metadata.send_block;
    snapshot.recv_axis = metadata.recv_axis;
    snapshot.send_axis = metadata.send_axis;
    snapshot.sign = metadata.sign;
    snapshot.recv_box = metadata.recv_box;
    snapshot.send_box = send_box;
    snapshot.node_box_on_send = metadata.node_box_on_send;
    snapshot.dir1_send = metadata.dir1_send;
    snapshot.dir2_send = metadata.dir2_send;
    snapshot.dir3_send = metadata.dir3_send;
    snapshot.trans_recv_to_send = metadata.trans_recv_to_send;
    snapshot.peer_rank =
        role == LegacyOrientedCornerRole::Receive
            ? metadata.send_rank
            : metadata.recv_rank;
    snapshot.tag = metadata.tag;
    snapshot.ncomp = ncomp;
    snapshot.nghost = nghost;
    snapshot.message_length =
        box_volume(metadata.recv_box) * ncomp;
    snapshot.local_plan_order = local_order;
    snapshot.peer_message_order =
        static_cast<std::size_t>(metadata.tag);
    return snapshot;
}
} // namespace

bool Halo::prepare_vertex_face_2form_metadata()
{
    int rank = 0;
    int size = 1;
    PARALLEL::mpi_rank(&rank);
    PARALLEL::mpi_size(&size);

    std::vector<PreparedGroupSpec> groups;
    groups.reserve(face_2form_triplets_.size());
    for (const auto &entry : face_2form_triplets_)
    {
        const HaloTripletRequest &triplet = entry.second;
        if (static_cast<int>(triplet.level) <
            static_cast<int>(HaloLevel::Corner3D))
            continue;

        PreparedGroupSpec group;
        group.group_name = triplet.group_name;
        group.fields = {triplet.xi, triplet.eta, triplet.zeta};
        group.level = triplet.level;
        group.nghost = triplet.nghost;
        group.orientation_aware = triplet.orientation_aware;
        bool descriptors_are_orientation_aware =
            triplet.orientation_aware;
        for (int axis = 0; axis < 3; ++axis)
        {
            const FieldDescriptor &descriptor =
                fld_->descriptor(group.fields[axis]);
            group.ncomp[axis] = descriptor.ncomp;
            group.locations[axis] = descriptor.location;
            group.value_kinds[axis] = descriptor.value_kind;
            descriptors_are_orientation_aware =
                descriptors_are_orientation_aware &&
                descriptor.sync.do_halo &&
                descriptor.sync.orientation_aware;
        }
        if (!descriptors_are_orientation_aware)
            continue;
        groups.push_back(std::move(group));
    }
    std::sort(
        groups.begin(), groups.end(),
        [](const PreparedGroupSpec &left,
           const PreparedGroupSpec &right)
        {
            return left.group_name < right.group_name;
        });

    const bool cache_valid =
        prepared_vertex_face_2form_metadata_valid_();
    verify_global_schedule(
        schedule_descriptor(
            "prepared-face2-vertex-v1",
            registered_patterns_built_, registry_generation_,
            pattern_generation_, cache_valid, groups),
        rank, size);

    if (!registered_patterns_built_)
        return false;
    if (cache_valid)
        return true;

    for (const PreparedGroupSpec &group : groups)
    {
        if (valid_face2_group(group))
            continue;

        std::ostringstream message;
        message << "[Halo][rank " << rank
                << "][group " << group.group_name
                << "][stage Vertex] cannot prepare an invalid "
                   "Face2Form triplet";
        ERROR::Abort(message.str());
    }

    std::unordered_map<std::string, PreparedVertexFace2FormPlanEntry>
        prepared;
    const auto &parallel_vertices =
        TOPO_VIEW::vertex_patches(*topo_, TOPO::PatchKind::Parallel);
    for (const PreparedGroupSpec &group : groups)
    {
        std::map<int, std::vector<PreparedTripletCornerMeta>>
            metadata_to_send;
        std::map<int, int> next_tag;
        std::vector<PreparedTripletCornerMeta> receive_plan;

        for (const TOPO::VertexPatch &patch : parallel_vertices)
        {
            if (patch.is_coupling)
                continue;

            const Direction this_direction_1 =
                HALO_TOOLS::int_to_direction(patch.dir1);
            const Direction this_direction_2 =
                HALO_TOOLS::int_to_direction(patch.dir2);
            const Direction this_direction_3 =
                HALO_TOOLS::int_to_direction(patch.dir3);
            int send_direction_1 =
                HALO_TOOLS::map_dir_to_neighbor(
                    patch.dir1, patch.trans, true);
            int send_direction_2 =
                HALO_TOOLS::map_dir_to_neighbor(
                    patch.dir2, patch.trans, false);
            int send_direction_3 =
                HALO_TOOLS::map_dir_to_neighbor(
                    patch.dir3, patch.trans, false);
            const int expected_normal_axis =
                patch.trans.perm[std::abs(patch.dir1) - 1] + 1;
            if (expected_normal_axis !=
                std::abs(send_direction_1))
            {
                if (expected_normal_axis ==
                    std::abs(send_direction_2))
                {
                    std::swap(send_direction_1, send_direction_2);
                }
                else if (expected_normal_axis ==
                         std::abs(send_direction_3))
                {
                    std::swap(send_direction_1, send_direction_3);
                }
                else
                {
                    std::ostringstream message;
                    message << "[Halo][rank " << rank
                            << "][group " << group.group_name
                            << "][stage Vertex] legacy direction "
                               "mapping is invalid";
                    ERROR::Abort(message.str());
                }
            }

            for (int recv_axis = 0; recv_axis < 3; ++recv_axis)
            {
                PreparedTripletCornerMeta metadata{};
                metadata.recv_rank = patch.this_rank;
                metadata.send_rank = patch.nb_rank;
                metadata.recv_block = patch.this_block;
                metadata.send_block = patch.nb_block;
                metadata.recv_axis = recv_axis;
                metadata.send_axis = inverse_perm_axis_for_prepared(
                    patch.trans, recv_axis);
                metadata.sign = HALO_TOOLS::face_2form_orientation_sign(
                    patch.trans, metadata.send_axis);
                metadata.recv_box =
                    HALO_BOX::make_3DCorner_ghost_box(
                        face_location(recv_axis), patch.this_box_node,
                        this_direction_1, this_direction_2,
                        this_direction_3, group.nghost);
                metadata.node_box_on_send = patch.nb_box_node;
                metadata.dir1_send = send_direction_1;
                metadata.dir2_send = send_direction_2;
                metadata.dir3_send = send_direction_3;
                metadata.trans_recv_to_send = patch.trans;
                metadata.tag = next_tag[patch.nb_rank]++;

                receive_plan.push_back(metadata);
                metadata_to_send[patch.nb_rank].push_back(metadata);
            }
        }

        std::vector<PreparedTripletCornerMeta> send_plan;
        // All ranks retain the legacy rank-grouped sender discovery order.
        exchange_prepared_triplet_corner_metadata(
            metadata_to_send, send_plan);

        LegacyFace2FormVertexPlanSnapshot snapshot;
        snapshot.group_name = group.group_name;
        snapshot.field_names = group.fields;
        snapshot.ncomp = group.ncomp[0];
        snapshot.nghost = group.nghost;
        snapshot.recv_plan.reserve(receive_plan.size());
        snapshot.send_plan.reserve(send_plan.size());

        for (std::size_t order = 0;
             order < send_plan.size(); ++order)
        {
            const PreparedTripletCornerMeta &metadata =
                send_plan[order];
            if (metadata.send_rank != rank)
            {
                std::ostringstream message;
                message << "[Halo][rank " << rank
                        << "][group " << group.group_name
                        << "][stage Vertex][role send] metadata "
                           "was delivered to the wrong rank";
                ERROR::Abort(message.str());
            }

            const Direction direction_1 =
                HALO_TOOLS::int_to_direction(metadata.dir1_send);
            const Direction direction_2 =
                HALO_TOOLS::int_to_direction(metadata.dir2_send);
            const Direction direction_3 =
                HALO_TOOLS::int_to_direction(metadata.dir3_send);
            const Box3 send_box =
                HALO_BOX::make_3DCorner_innerghost_box(
                    face_location(metadata.send_axis),
                    metadata.node_box_on_send,
                    direction_1, direction_2, direction_3,
                    group.nghost);
            snapshot.send_plan.push_back(make_face2_vertex_snapshot(
                metadata, LegacyOrientedCornerRole::Send,
                std::optional<Box3>{send_box}, group.ncomp[0],
                group.nghost, order));
        }

        for (std::size_t order = 0;
             order < receive_plan.size(); ++order)
        {
            snapshot.recv_plan.push_back(make_face2_vertex_snapshot(
                receive_plan[order],
                LegacyOrientedCornerRole::Receive,
                std::nullopt, group.ncomp[0], group.nghost,
                order));
        }

        PreparedVertexFace2FormPlanEntry entry;
        entry.snapshot = std::move(snapshot);
        entry.registry_generation = registry_generation_;
        entry.pattern_generation = pattern_generation_;
        prepared.emplace(group.group_name, std::move(entry));
    }

    prepared_vertex_face_2form_plans_ = std::move(prepared);
    prepared_vertex_face_2form_metadata_ready_ = true;
    prepared_vertex_face_2form_registry_generation_ =
        registry_generation_;
    prepared_vertex_face_2form_pattern_generation_ =
        pattern_generation_;
    prepared_vertex_face_2form_metadata_exchange_count_ +=
        static_cast<std::uint64_t>(groups.size());
    return true;
}
