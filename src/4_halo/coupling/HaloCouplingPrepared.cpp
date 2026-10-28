#include "4_halo/Halo.h"

#include "0_basic/Error.h"
#include "0_basic/LayoutTraits.h"
#include "1_grid/BlockTraits.h"
#include "2_topology/TopologyView.h"
#include "4_halo/detail/HaloBuildBoxMakers.h"
#include "4_halo/detail/HaloBuildTools.h"

#include <algorithm>
#include <cstdint>
#include <vector>

namespace
{
int coupling_box_length(const Box3 &box, int ncomp)
{
    const int ni = std::max(0, box.hi.i - box.lo.i);
    const int nj = std::max(0, box.hi.j - box.lo.j);
    const int nk = std::max(0, box.hi.k - box.lo.k);
    return ni * nj * nk * ncomp;
}

void pack_coupling_transfer(FieldBlock &field, const Box3 &box,
                            int ncomp, const TOPO::IndexTransform &transform,
                            double factor, double *buffer, int length)
{
    if (!field.is_allocated() || length == 0)
        return;

    const int local_length[3] = {
        box.hi.i - box.lo.i,
        box.hi.j - box.lo.j,
        box.hi.k - box.lo.k};
    if (local_length[0] <= 0 || local_length[1] <= 0 ||
        local_length[2] <= 0)
        return;
    if (local_length[0] * local_length[1] * local_length[2] * ncomp !=
        length)
        ERROR::Abort("[Halo] Prepared Coupling send length mismatch");

    const Int3 field_lo = field.get_lo();
    const Int3 field_hi = field.get_hi();
    const int source_lo[3] = {box.lo.i, box.lo.j, box.lo.k};
    const int source_hi[3] = {box.hi.i - 1, box.hi.j - 1,
                              box.hi.k - 1};
    const int offset[3] = {transform.offset.i, transform.offset.j,
                           transform.offset.k};
    int target_lo[3] = {0, 0, 0};
    int target_hi[3] = {0, 0, 0};
    int target_ref[3] = {0, 0, 0};
    int target_length[3] = {0, 0, 0};
    for (int axis = 0; axis < 3; ++axis)
    {
        target_lo[transform.perm[axis]] =
            transform.sign[axis] * source_lo[axis] + offset[axis];
        target_hi[transform.perm[axis]] =
            transform.sign[axis] * source_hi[axis] + offset[axis];
        target_length[transform.perm[axis]] = local_length[axis];
    }
    for (int axis = 0; axis < 3; ++axis)
        target_ref[axis] = std::min(target_lo[axis], target_hi[axis]);

    int source[3], target[3];
    for (source[0] = box.lo.i; source[0] < box.hi.i; ++source[0])
        for (source[1] = box.lo.j; source[1] < box.hi.j; ++source[1])
            for (source[2] = box.lo.k; source[2] < box.hi.k; ++source[2])
            {
                if (source[0] < field_lo.i || source[0] >= field_hi.i ||
                    source[1] < field_lo.j || source[1] >= field_hi.j ||
                    source[2] < field_lo.k || source[2] >= field_hi.k)
                    ERROR::Abort("[Halo] Prepared Coupling source box is out of bounds");
                for (int axis = 0; axis < 3; ++axis)
                    target[transform.perm[axis]] =
                        transform.sign[axis] * source[axis] + offset[axis];
                const int ri = target[0] - target_ref[0];
                const int rj = target[1] - target_ref[1];
                const int rk = target[2] - target_ref[2];
                const int base =
                    ((ri * target_length[1] + rj) * target_length[2] +
                     rk) * ncomp;
                if (base < 0 || base + ncomp > length)
                    ERROR::Abort("[Halo] Prepared Coupling pack is out of bounds");
                for (int component = 0; component < ncomp; ++component)
                    buffer[base + component] =
                        factor * field(source[0], source[1], source[2],
                                       component);
            }
}

void unpack_coupling_transfer(CouplingBufferBlock &destination,
                              const Box3 &box, int ncomp,
                              const double *buffer, int length)
{
    if (!destination.allocated || length == 0)
        return;

    const int ni = box.hi.i - box.lo.i;
    const int nj = box.hi.j - box.lo.j;
    const int nk = box.hi.k - box.lo.k;
    if (ni * nj * nk * ncomp != length)
        ERROR::Abort("[Halo] Prepared Coupling receive length mismatch");
    if (!HALO_TOOLS::box_equal(destination.box, box))
        ERROR::Abort("[Halo] Prepared Coupling destination layout changed");

    for (int ii = 0; ii < ni; ++ii)
        for (int jj = 0; jj < nj; ++jj)
            for (int kk = 0; kk < nk; ++kk)
            {
                const int base =
                    ((ii * nj + jj) * nk + kk) * ncomp;
                for (int component = 0; component < ncomp; ++component)
                    destination(box.lo.i + ii, box.lo.j + jj,
                                box.lo.k + kk, component) =
                        buffer[base + component];
            }
}

std::uint64_t fnv_append(std::uint64_t hash, const void *data,
                         std::size_t size)
{
    const unsigned char *bytes =
        static_cast<const unsigned char *>(data);
    for (std::size_t i = 0; i < size; ++i)
    {
        hash ^= bytes[i];
        hash *= 1099511628211ULL;
    }
    return hash;
}

std::uint64_t fnv_string(std::uint64_t hash, const std::string &value)
{
    hash = fnv_append(hash, value.data(), value.size());
    const unsigned char separator = 0xff;
    return fnv_append(hash, &separator, 1);
}
}

void Halo::validate_coupling_registry_consistency_()
{
    std::uint64_t hash = 1469598103934665603ULL;
    std::uint64_t count = 0;
    for (const auto &pair_entry : fld_->coupling_pairs())
    {
        hash = fnv_string(hash, pair_entry.first.first);
        hash = fnv_string(hash, pair_entry.first.second);
        for (const CouplingChannelSpec &channel :
             pair_entry.second.channels)
        {
            hash = fnv_string(hash, channel.tag);
            const int values[] = {
                static_cast<int>(channel.location),
                static_cast<int>(channel.value_kind),
                channel.ncomp, channel.nghost,
                channel.orientation_aware ? 1 : 0};
            hash = fnv_append(hash, values, sizeof(values));
            ++count;
        }
    }

    unsigned long long local[2] = {
        static_cast<unsigned long long>(count),
        static_cast<unsigned long long>(hash)};
    int size = 1;
    PARALLEL::mpi_size(&size);
    std::vector<unsigned long long> gathered(
        static_cast<std::size_t>(size) * 2, 0);
    MPI_Allgather(local, 2, MPI_UNSIGNED_LONG_LONG,
                  gathered.data(), 2, MPI_UNSIGNED_LONG_LONG,
                  MPI_COMM_WORLD);
    ++coupling_execution_counters_.initialization_collectives;

    for (int rank = 1; rank < size; ++rank)
        if (gathered[2 * rank] != gathered[0] ||
            gathered[2 * rank + 1] != gathered[1])
            ERROR::Abort(
                "[Halo] Coupling registry differs across MPI ranks");
}

int Halo::managed_coupling_tag_(int original_tag) const
{
    if (original_tag < 0 || coupling_tag_base_ < 0)
        ERROR::Abort("[Halo] invalid managed Coupling tag");
    return coupling_tag_base_ + original_tag;
}

std::vector<CouplingBufferBlock> &Halo::coupling_buffer_list_(
    CouplingBuffersForPair &buffers, HaloLevel stage,
    bool parallel, int channel) const
{
    if (channel < 0 ||
        channel >= static_cast<int>(buffers.desc.channels.size()))
        ERROR::Abort("[Halo] invalid Coupling channel id");

    if (stage == HaloLevel::Corner1D)
        return parallel ? buffers.parallel_face[channel]
                        : buffers.inner_face[channel];
    if (stage == HaloLevel::Corner2D)
        return parallel ? buffers.parallel_edge[channel]
                        : buffers.inner_edge[channel];
    if (stage == HaloLevel::Corner3D)
        return parallel ? buffers.parallel_vertex[channel]
                        : buffers.inner_vertex[channel];
    ERROR::Abort("[Halo] invalid Coupling stage");
}

const Halo::ParallelCouplingPlan &Halo::parallel_coupling_plan_(
    const std::string &src, const std::string &dst,
    HaloLevel stage, int channel) const
{
    const auto found = parallel_coupling_plans_.find(
        CouplingPlanKey{src, dst, static_cast<int>(stage), channel});
    if (found == parallel_coupling_plans_.end())
        ERROR::Abort("[Halo] Prepared parallel Coupling plan is missing");
    return found->second;
}

const std::vector<Halo::LocalCouplingTransfer> &
Halo::local_coupling_plan_(const std::string &src,
                           const std::string &dst,
                           HaloLevel stage, int channel) const
{
    const auto found = local_coupling_plans_.find(
        CouplingPlanKey{src, dst, static_cast<int>(stage), channel});
    if (found == local_coupling_plans_.end())
        ERROR::Abort("[Halo] Prepared local Coupling plan is missing");
    return found->second;
}

void Halo::build_coupling_execution_plans_(std::uint64_t generation)
{
    local_coupling_plans_.clear();
    parallel_coupling_plans_.clear();
    coupling_plans_ready_ = false;

    const int dimension = fld_->grd->dimension;
    const auto inner_faces = TOPO_VIEW::inner_faces(*topo_);
    const auto parallel_faces = TOPO_VIEW::parallel_faces(*topo_);
    const auto &inner_edges =
        TOPO_VIEW::edge_patches(*topo_, TOPO::PatchKind::Inner);
    const auto &inner_vertices =
        TOPO_VIEW::vertex_patches(*topo_, TOPO::PatchKind::Inner);

    for (const auto &pair_entry : fld_->coupling_pairs())
    {
        const std::string &src = pair_entry.first.first;
        const std::string &dst = pair_entry.first.second;
        const CouplingPairDesc &description = pair_entry.second;
        CouplingBuffersForPair &buffers =
            fld_->coupling_buffers(src, dst);

        const auto build_local_stage =
            [&](HaloLevel stage, const auto &patches,
                auto CouplingBuffersForPair::*storage_member)
        {
            auto &storage = buffers.*storage_member;
            for (int channel = 0;
                 channel < static_cast<int>(description.channels.size());
                 ++channel)
            {
                const CouplingChannelSpec &spec =
                    description.channels[channel];
                std::vector<LocalCouplingTransfer> transfers;
                for (std::size_t patch_index = 0;
                     patch_index < patches.size(); ++patch_index)
                {
                    const auto &patch = patches[patch_index];
                    if (!patch.is_coupling ||
                        patch.this_block_name != dst ||
                        patch.nb_block_name != src)
                        continue;
                    if (channel >= static_cast<int>(storage.size()) ||
                        patch_index >= storage[channel].size())
                        ERROR::Abort(
                            "[Halo] Coupling local buffer layout mismatch");
                    CouplingBufferBlock &buffer =
                        storage[channel][patch_index];
                    if (!buffer.allocated)
                        continue;

                    std::string source_name =
                        buffer.tag.empty() ? spec.tag : buffer.tag;
                    double factor = 1.0;
                    if (coupling_channel_needs_form_transfer_(spec))
                    {
                        const int destination_axis =
                            coupling_form_axis_from_location_(spec.location);
                        const int source_axis =
                            coupling_src_axis_from_dst_to_src_transform_(
                                patch.trans, destination_axis);
                        source_name = find_triplet_field_name_(
                            source_name, spec.value_kind, source_axis);
                        factor = static_cast<double>(
                            coupling_form_orientation_sign_dst_to_src_(
                                spec, patch.trans, destination_axis));
                    }
                    const int source_fid =
                        fld_->field_id(source_name);
                    if (!field_active_(source_fid, patch.nb_block))
                        continue;

                    transfers.push_back({
                        source_fid, patch.nb_block, patch_index,
                        buffer.box, patch.trans, buffer.ncomp, factor});
                }
                local_coupling_plans_.emplace(
                    CouplingPlanKey{src, dst, static_cast<int>(stage),
                                    channel},
                    std::move(transfers));
            }
        };

        build_local_stage(HaloLevel::Corner1D, inner_faces,
                          &CouplingBuffersForPair::inner_face);
        if (dimension >= 2)
            build_local_stage(HaloLevel::Corner2D, inner_edges,
                              &CouplingBuffersForPair::inner_edge);
        if (dimension >= 3)
            build_local_stage(HaloLevel::Corner3D, inner_vertices,
                              &CouplingBuffersForPair::inner_vertex);

        const auto finish_parallel_plan =
            [&](HaloLevel stage, int channel,
                ParallelCouplingPlan &&prepared)
        {
            prepared.exchange_plan.name =
                "coupling:" + src + "->" + dst + ":" +
                std::to_string(static_cast<int>(stage)) + ":" +
                std::to_string(channel);
            prepared.exchange_plan.generation = generation;
            collectively_match_exchange_plan_(
                prepared.exchange_plan);
            coupling_execution_counters_.initialization_collectives += 2;
            for (const auto &message : prepared.exchange_plan.sends)
                if (message.length == 0)
                    ++coupling_execution_counters_.zero_length_messages;
            for (const auto &message : prepared.exchange_plan.recvs)
                if (message.length == 0)
                    ++coupling_execution_counters_.zero_length_messages;
            prepared_exchange_executor_.reserve_plan(
                prepared.exchange_plan);
            parallel_coupling_plans_.emplace(
                CouplingPlanKey{src, dst, static_cast<int>(stage),
                                channel},
                std::move(prepared));
        };

        for (int channel = 0;
             channel < static_cast<int>(description.channels.size());
             ++channel)
        {
            const CouplingChannelSpec &spec =
                description.channels[channel];
            ParallelCouplingPlan prepared;
            for (std::size_t patch_index = 0;
                 patch_index < parallel_faces.size(); ++patch_index)
            {
                const auto &patch = parallel_faces[patch_index];
                if (!patch.is_coupling)
                    continue;

                if (patch.this_block_name == src &&
                    patch.nb_block_name == dst)
                {
                    std::string source_name = spec.tag;
                    StaggerLocation source_location = spec.location;
                    double factor = 1.0;
                    if (coupling_channel_needs_form_transfer_(spec))
                    {
                        const int destination_axis =
                            coupling_form_axis_from_location_(spec.location);
                        const int source_axis =
                            coupling_src_axis_from_src_to_dst_transform_(
                                patch.trans, destination_axis);
                        source_name = find_triplet_field_name_(
                            spec.tag, spec.value_kind, source_axis);
                        source_location =
                            coupling_form_location_from_axis_(
                                spec.value_kind, source_axis);
                        factor = static_cast<double>(
                            coupling_form_orientation_sign_src_to_dst_(
                                spec, patch.trans, source_axis));
                    }
                    const int source_fid =
                        fld_->field_id(source_name);
                    const Int3 block_nodes =
                        HALO_TOOLS::block_node_size(
                            fld_->grd->grids(patch.this_block));
                    const Direction direction =
                        HALO_TOOLS::detect_face_direction(
                            patch.this_box_node, block_nodes,
                            "build Prepared Coupling Face");
                    const Box3 source_box =
                        HALO_BOX::make_1DCorner_inner_box(
                            source_location, patch.this_box_node,
                            direction, spec.nghost);
                    prepared.sends.push_back({
                        source_fid, patch.this_block, source_box,
                        patch.trans, spec.ncomp, factor});
                    prepared.exchange_plan.sends.push_back({
                        patch.nb_rank,
                        managed_coupling_tag_(patch.send_flag),
                        coupling_box_length(source_box, spec.ncomp), 0,
                        field_active_(source_fid, patch.this_block)});
                }

                if (patch.this_block_name == dst &&
                    patch.nb_block_name == src)
                {
                    if (channel >=
                            static_cast<int>(buffers.parallel_face.size()) ||
                        patch_index >=
                            buffers.parallel_face[channel].size())
                        ERROR::Abort(
                            "[Halo] Coupling Face buffer layout mismatch");
                    CouplingBufferBlock &buffer =
                        buffers.parallel_face[channel][patch_index];
                    Box3 destination_box = buffer.box;
                    if (!buffer.allocated)
                    {
                        const Block &block =
                            fld_->grd->grids(patch.this_block);
                        destination_box =
                            LAYOUT::coupling_face_ghost_slab_from_cells(
                                GRID_TRAITS::cell_counts(block),
                                spec.location, patch.this_box_node,
                                patch.direction, spec.nghost);
                    }
                    prepared.recvs.push_back({
                        patch_index, destination_box, spec.ncomp});
                    prepared.exchange_plan.recvs.push_back({
                        patch.nb_rank,
                        managed_coupling_tag_(patch.recv_flag),
                        coupling_box_length(destination_box, spec.ncomp), 0,
                        buffer.allocated});
                }
            }
            finish_parallel_plan(HaloLevel::Corner1D, channel,
                                 std::move(prepared));
        }

        const auto build_parallel_corner_stage =
            [&](HaloLevel stage, const auto &send_patterns,
                const auto &recv_patterns,
                auto CouplingBuffersForPair::*storage_member)
        {
            auto &storage = buffers.*storage_member;
            for (int channel = 0;
                 channel < static_cast<int>(description.channels.size());
                 ++channel)
            {
                const CouplingChannelSpec &spec =
                    description.channels[channel];
                const CouplingPatternKey base_key{
                    src, dst, spec.location, spec.nghost};
                const auto recv_found = recv_patterns.find(base_key);
                if (recv_found == recv_patterns.end())
                    ERROR::Abort(
                        "[Halo] Coupling corner receive pattern is missing");

                std::vector<const HaloPattern *> source_patterns;
                if (coupling_channel_needs_form_transfer_(spec))
                {
                    for (int axis = 0; axis < 3; ++axis)
                    {
                        const CouplingPatternKey source_key{
                            src, dst,
                            coupling_form_location_from_axis_(
                                spec.value_kind, axis),
                            spec.nghost};
                        const auto found =
                            send_patterns.find(source_key);
                        if (found == send_patterns.end())
                            ERROR::Abort(
                                "[Halo] Coupling form send pattern is missing");
                        source_patterns.push_back(&found->second);
                    }
                }
                else
                {
                    const auto found =
                        send_patterns.find(base_key);
                    if (found == send_patterns.end())
                        ERROR::Abort(
                            "[Halo] Coupling corner send pattern is missing");
                    source_patterns.push_back(&found->second);
                }

                ParallelCouplingPlan prepared;
                for (const HaloPattern *source_pattern :
                     source_patterns)
                    for (const HaloRegion &region :
                         source_pattern->regions)
                    {
                        std::string source_name = spec.tag;
                        double factor = 1.0;
                        if (coupling_channel_needs_form_transfer_(spec))
                        {
                            const int destination_axis =
                                coupling_form_axis_from_location_(
                                    spec.location);
                            const int source_axis =
                                coupling_src_axis_from_src_to_dst_transform_(
                                    region.trans, destination_axis);
                            if (coupling_form_location_from_axis_(
                                    spec.value_kind, source_axis) !=
                                source_pattern->location)
                                continue;
                            source_name = find_triplet_field_name_(
                                spec.tag, spec.value_kind, source_axis);
                            factor = static_cast<double>(
                                coupling_form_orientation_sign_src_to_dst_(
                                    spec, region.trans, source_axis));
                        }
                        const int source_fid =
                            fld_->field_id(source_name);
                        prepared.sends.push_back({
                            source_fid, region.this_block,
                            region.send_box, region.trans,
                            spec.ncomp, factor});
                        prepared.exchange_plan.sends.push_back({
                            region.neighbor_rank,
                            managed_coupling_tag_(region.send_flag),
                            coupling_box_length(region.send_box,
                                                spec.ncomp),
                            0,
                            field_active_(source_fid,
                                          region.this_block)});
                    }

                for (const HaloRegion &region :
                     recv_found->second.regions)
                {
                    std::size_t buffer_index = 0;
                    bool allocated = false;
                    if (channel < static_cast<int>(storage.size()))
                        for (std::size_t index = 0;
                             index < storage[channel].size(); ++index)
                        {
                            CouplingBufferBlock &buffer =
                                storage[channel][index];
                            if (!buffer.allocated ||
                                buffer.this_block != region.this_block ||
                                buffer.location != spec.location ||
                                buffer.nghost != spec.nghost ||
                                !HALO_TOOLS::box_equal(
                                    buffer.box, region.recv_box))
                                continue;
                            buffer_index = index;
                            allocated = true;
                            break;
                        }
                    prepared.recvs.push_back({
                        buffer_index, region.recv_box, spec.ncomp});
                    prepared.exchange_plan.recvs.push_back({
                        region.neighbor_rank,
                        managed_coupling_tag_(region.recv_flag),
                        coupling_box_length(region.recv_box,
                                            spec.ncomp),
                        0, allocated});
                }
                finish_parallel_plan(stage, channel,
                                     std::move(prepared));
            }
        };

        if (dimension >= 2)
            build_parallel_corner_stage(
                HaloLevel::Corner2D,
                coupling_parallel_edge_patterns_send,
                coupling_parallel_edge_patterns_recv,
                &CouplingBuffersForPair::parallel_edge);
        if (dimension >= 3)
            build_parallel_corner_stage(
                HaloLevel::Corner3D,
                coupling_parallel_vertex_patterns_send,
                coupling_parallel_vertex_patterns_recv,
                &CouplingBuffersForPair::parallel_vertex);
    }

    coupling_field_generation_ = fld_->coupling_generation();
    coupling_pattern_generation_ = generation;
    coupling_plans_ready_ = true;
}

void Halo::validate_coupling_execution_ready_() const
{
    if (!registered_patterns_built_ || !coupling_plans_ready_ ||
        coupling_pattern_generation_ != pattern_generation_ ||
        coupling_field_generation_ != fld_->coupling_generation() ||
        !prepared_exchange_executor_.ready(pattern_generation_))
        ERROR::Abort(
            "[Halo] Coupling plans are stale; explicitly rebuild Halo patterns");
}

void Halo::execute_coupling_inner_(
    const std::string &src, const std::string &dst,
    HaloLevel stage, const std::vector<int32_t> *channels)
{
    if (static_cast<int>(stage) > fld_->grd->dimension)
        return;
    if (!fld_->has_coupling_pair(src, dst))
        return;
    validate_coupling_execution_ready_();

    CouplingBuffersForPair &buffers =
        fld_->coupling_buffers(src, dst);
    const int channel_count =
        static_cast<int>(buffers.desc.channels.size());
    const auto execute_channel = [&](int channel)
    {
        if (channel < 0 || channel >= channel_count)
            ERROR::Abort("[Halo] invalid Coupling channel subset");
        std::vector<CouplingBufferBlock> &destination_buffers =
            coupling_buffer_list_(buffers, stage, false, channel);
        const auto &transfers =
            local_coupling_plan_(src, dst, stage, channel);
        for (const LocalCouplingTransfer &transfer : transfers)
        {
            if (transfer.buffer_index >= destination_buffers.size())
                ERROR::Abort(
                    "[Halo] stale local Coupling buffer index");
            CouplingBufferBlock &destination =
                destination_buffers[transfer.buffer_index];
            if (!destination.allocated ||
                !HALO_TOOLS::box_equal(
                    destination.box, transfer.destination_box))
                ERROR::Abort(
                    "[Halo] local Coupling buffer layout changed");
            FieldBlock &source = fld_->field(
                transfer.source_fid, transfer.source_block);
            const Box3 &box = transfer.destination_box;
            for (int i = box.lo.i; i < box.hi.i; ++i)
                for (int j = box.lo.j; j < box.hi.j; ++j)
                    for (int k = box.lo.k; k < box.hi.k; ++k)
                    {
                        int si = 0, sj = 0, sk = 0;
                        HALO_TOOLS::apply_transform(
                            transfer.destination_to_source,
                            i, j, k, si, sj, sk);
                        for (int component = 0;
                             component < transfer.ncomp; ++component)
                            destination(i, j, k, component) =
                                transfer.factor *
                                source(si, sj, sk, component);
                    }
        }
        ++coupling_execution_counters_.
            prepared_inner_channel_executions;
    };
    if (channels)
        for (const int32_t channel : *channels)
            execute_channel(static_cast<int>(channel));
    else
        for (int channel = 0; channel < channel_count; ++channel)
            execute_channel(channel);
}

void Halo::execute_coupling_parallel_(
    const std::string &src, const std::string &dst,
    HaloLevel stage, const std::vector<int32_t> *channels)
{
    if (static_cast<int>(stage) > fld_->grd->dimension)
        return;
    if (!fld_->has_coupling_pair(src, dst))
        return;
    validate_coupling_execution_ready_();

    CouplingBuffersForPair &buffers =
        fld_->coupling_buffers(src, dst);
    const int channel_count =
        static_cast<int>(buffers.desc.channels.size());
    const auto execute_channel = [&](int channel)
    {
        if (channel < 0 || channel >= channel_count)
            ERROR::Abort("[Halo] invalid Coupling channel subset");
        std::vector<CouplingBufferBlock> &destination_buffers =
            coupling_buffer_list_(buffers, stage, true, channel);
        const ParallelCouplingPlan &prepared =
            parallel_coupling_plan_(src, dst, stage, channel);

        prepared_exchange_executor_.execute(
            prepared.exchange_plan,
            [&](std::size_t index, double *buffer, int length)
            {
                const ParallelCouplingSendTransfer &transfer =
                    prepared.sends[index];
                FieldBlock &source = fld_->field(
                    transfer.source_fid, transfer.source_block);
                pack_coupling_transfer(
                    source, transfer.source_box, transfer.ncomp,
                    transfer.source_to_destination, transfer.factor,
                    buffer, length);
            },
            [&](std::size_t index, const double *buffer, int length)
            {
                const ParallelCouplingRecvTransfer &transfer =
                    prepared.recvs[index];
                if (transfer.buffer_index >=
                    destination_buffers.size())
                    ERROR::Abort(
                        "[Halo] stale parallel Coupling buffer index");
                unpack_coupling_transfer(
                    destination_buffers[transfer.buffer_index],
                    transfer.destination_box, transfer.ncomp,
                    buffer, length);
            });

        for (const auto &message : prepared.exchange_plan.sends)
            if (message.enabled)
            {
                ++coupling_execution_counters_.sent_messages;
                coupling_execution_counters_.sent_doubles +=
                    message.length;
            }
        for (const auto &message : prepared.exchange_plan.recvs)
            if (message.enabled)
            {
                ++coupling_execution_counters_.received_messages;
                coupling_execution_counters_.received_doubles +=
                    message.length;
            }
        ++coupling_execution_counters_.
            prepared_parallel_channel_executions;
    };
    if (channels)
        for (const int32_t channel : *channels)
            execute_channel(static_cast<int>(channel));
    else
        for (int channel = 0; channel < channel_count; ++channel)
            execute_channel(channel);
}

void Halo::set_coupling_execution_mode(
    CouplingExecutionMode mode) noexcept
{
    coupling_execution_mode_ = mode;
}

CouplingExecutionMode
Halo::debug_coupling_execution_mode() const noexcept
{
    return coupling_execution_mode_;
}

CouplingExecutionCounters
Halo::debug_coupling_execution_counters() const noexcept
{
    return coupling_execution_counters_;
}

bool Halo::debug_coupling_execution_ready() const noexcept
{
    return registered_patterns_built_ && coupling_plans_ready_ &&
           coupling_pattern_generation_ == pattern_generation_ &&
           coupling_field_generation_ == fld_->coupling_generation() &&
           prepared_exchange_executor_.ready(pattern_generation_);
}

#define HALO_COUPLING_WRAPPERS(STAGE_NAME, LEVEL)                         \
void Halo::coupling_inner_##STAGE_NAME(                                  \
    const std::string &src, const std::string &dst)                       \
{                                                                         \
    if (coupling_execution_mode_ == CouplingExecutionMode::LegacyValidation) \
    {                                                                     \
        ++coupling_execution_counters_.legacy_entry_executions;          \
        coupling_inner_##STAGE_NAME##_legacy_(src, dst);                  \
        return;                                                           \
    }                                                                     \
    execute_coupling_inner_(src, dst, LEVEL, nullptr);                    \
}                                                                         \
void Halo::coupling_inner_##STAGE_NAME(                                  \
    const std::string &src, const std::string &dst,                       \
    std::vector<int32_t> &channels)                                      \
{                                                                         \
    if (coupling_execution_mode_ == CouplingExecutionMode::LegacyValidation) \
    {                                                                     \
        ++coupling_execution_counters_.legacy_entry_executions;          \
        coupling_inner_##STAGE_NAME##_legacy_(src, dst, channels);        \
        return;                                                           \
    }                                                                     \
    execute_coupling_inner_(src, dst, LEVEL, &channels);                  \
}                                                                         \
void Halo::coupling_parallel_##STAGE_NAME(                               \
    const std::string &src, const std::string &dst)                       \
{                                                                         \
    if (coupling_execution_mode_ == CouplingExecutionMode::LegacyValidation) \
    {                                                                     \
        ++coupling_execution_counters_.legacy_entry_executions;          \
        coupling_parallel_##STAGE_NAME##_legacy_(src, dst);               \
        return;                                                           \
    }                                                                     \
    execute_coupling_parallel_(src, dst, LEVEL, nullptr);                 \
}                                                                         \
void Halo::coupling_parallel_##STAGE_NAME(                               \
    const std::string &src, const std::string &dst,                       \
    std::vector<int32_t> &channels)                                      \
{                                                                         \
    if (coupling_execution_mode_ == CouplingExecutionMode::LegacyValidation) \
    {                                                                     \
        ++coupling_execution_counters_.legacy_entry_executions;          \
        coupling_parallel_##STAGE_NAME##_legacy_(src, dst, channels);     \
        return;                                                           \
    }                                                                     \
    execute_coupling_parallel_(src, dst, LEVEL, &channels);               \
}

HALO_COUPLING_WRAPPERS(face, HaloLevel::Corner1D)
HALO_COUPLING_WRAPPERS(edge, HaloLevel::Corner2D)
HALO_COUPLING_WRAPPERS(vertex, HaloLevel::Corner3D)

#undef HALO_COUPLING_WRAPPERS
