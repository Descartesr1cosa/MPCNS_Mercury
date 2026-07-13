#pragma once
#include "2_topology/Topology.h"
#include "4_halo/HaloTypes.h"
#include "4_halo/HaloSyncTypes.h"
#include "3_field/Field.h"
#include "0_basic/MPI_WRAPPER.h"

#include <iosfwd>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>


class Halo
{
public:
    Halo(Field *field, TOPO::Topology *topo);

    // Registration and pattern lifecycle.
    void register_halo_field(const FieldHaloRequest &request);
    void register_halo_field(const std::string &field_name,
                             HaloLevel level = HaloLevel::Corner3D);
    void register_halo_fields(const std::vector<FieldHaloRequest> &requests);

    void build_registered_patterns();

    // High-level sync API.
    void sync_registered();
    void sync_registered(HaloLevel stage);
    void sync_field(const std::string &field_name);
    void sync_field(const std::string &field_name, HaloLevel stage);
    void sync_group(const std::string &group_name);
    void sync_group(const std::string &group_name, HaloLevel stage);
    void sync_owner_alias();
    void sync_owner_alias_field(const std::string &field_name);
    void sync_owner_alias_group(const std::string &group_name);
    void dump_sync_registry(std::ostream &os) const;

    // Test/debug helper: expose the recv/send boxes used by the registered
    // halo patterns without changing the production sync path.
    std::vector<HaloRegion> debug_halo_regions(StaggerLocation location,
                                                int nghost,
                                                HaloLevel stage) const;
    std::vector<HaloRegion> debug_halo_send_regions(StaggerLocation location,
                                                     int nghost,
                                                     HaloLevel stage) const;

    // Topology equivalence is required by owner-alias synchronization.
    void set_topology_equiv(const TOPO::Topology *equiv);
    const TOPO::Topology *topology_equiv() const;

    // Legacy halo transfer API. Prefer sync_registered/sync_field for new code.
    void data_trans_1DCorner(std::string &field_name);
    void data_trans_2DCorner(std::string &field_name);
    void data_trans_3DCorner(std::string &field_name);
    void data_trans_edge_1form_triplet(const std::vector<std::string> &fields,
                                       HaloLevel stage);
    void data_trans_face_axis_triplet(const std::vector<std::string> &fields,
                                      HaloLevel stage);

    // Coupling transfer API. These fill coupling buffers and do not write ghosts.
    void coupling_trans_1DCorner(std::string &src, std::string &dst);
    void coupling_trans_1DCorner(std::string &src, std::string &dst, std::vector<int32_t> &field_cids);
    void coupling_trans_2DCorner(std::string &src, std::string &dst);
    void coupling_trans_2DCorner(std::string &src, std::string &dst, std::vector<int32_t> &field_cids);
    void coupling_trans_3DCorner(std::string &src, std::string &dst);
    void coupling_trans_3DCorner(std::string &src, std::string &dst, std::vector<int32_t> &field_cids);

private:
    // Same-field component copy exchange.
    void exchange_inner(const std::string &field_name);
    void exchange_parallel(const std::string &field_name);
    void exchange_inner_edge(const std::string &field_name);
    void exchange_parallel_edge(const std::string &field_name);
    void exchange_inner_vertex(const std::string &field_name);
    void exchange_parallel_vertex(const std::string &field_name);

    // Orientation-aware triplet exchange.
    void exchange_inner_face_edge_1form_triplet_(const std::vector<std::string> &fields);
    void exchange_parallel_face_edge_1form_triplet_(const std::vector<std::string> &fields);
    void exchange_inner_edge_edge_1form_triplet_(const std::vector<std::string> &fields);
    void exchange_parallel_edge_edge_1form_triplet_(const std::vector<std::string> &fields);
    void exchange_inner_vertex_edge_1form_triplet_(const std::vector<std::string> &fields);
    void exchange_parallel_vertex_edge_1form_triplet_(const std::vector<std::string> &fields);

    void exchange_inner_face_face_2form_triplet_(const std::vector<std::string> &fields, bool oriented = true);
    void exchange_parallel_face_face_2form_triplet_(const std::vector<std::string> &fields, bool oriented = true);
    void exchange_inner_edge_face_2form_triplet_(const std::vector<std::string> &fields, bool oriented = true);
    void exchange_parallel_edge_face_2form_triplet_(const std::vector<std::string> &fields, bool oriented = true);
    void exchange_inner_vertex_face_2form_triplet_(const std::vector<std::string> &fields, bool oriented = true);
    void exchange_parallel_vertex_face_2form_triplet_(const std::vector<std::string> &fields, bool oriented = true);

    // Coupling buffer exchange.
    void coupling_inner_face(const std::string &src, const std::string &dst);
    void coupling_parallel_face(const std::string &src, const std::string &dst);
    void coupling_inner_face(const std::string &src, const std::string &dst, std::vector<int32_t> &field_cids);
    void coupling_parallel_face(const std::string &src, const std::string &dst, std::vector<int32_t> &field_cids);
    void coupling_inner_edge(const std::string &src, const std::string &dst);
    void coupling_parallel_edge(const std::string &src, const std::string &dst);
    void coupling_inner_edge(const std::string &src, const std::string &dst, std::vector<int32_t> &field_cids);
    void coupling_parallel_edge(const std::string &src, const std::string &dst, std::vector<int32_t> &field_cids);
    void coupling_inner_vertex(const std::string &src, const std::string &dst);
    void coupling_parallel_vertex(const std::string &src, const std::string &dst);
    void coupling_inner_vertex(const std::string &src, const std::string &dst, std::vector<int32_t> &field_cids);
    void coupling_parallel_vertex(const std::string &src, const std::string &dst, std::vector<int32_t> &field_cids);

    // Orientation-aware coupling helpers.
    bool coupling_channel_needs_form_transfer_(const CouplingChannelSpec &ch) const;
    int coupling_form_axis_from_location_(StaggerLocation loc) const;
    StaggerLocation coupling_form_location_from_axis_(FieldValueKind value_kind, int axis) const;
    int coupling_src_axis_from_dst_to_src_transform_(const TOPO::IndexTransform &tr, int dst_axis) const;
    int coupling_src_axis_from_src_to_dst_transform_(const TOPO::IndexTransform &tr, int dst_axis) const;
    int coupling_edge_1form_sign_dst_to_src_(const TOPO::IndexTransform &tr, int dst_axis) const;
    int coupling_edge_1form_sign_src_to_dst_(const TOPO::IndexTransform &tr, int src_axis) const;
    int coupling_face_2form_sign_dst_to_src_(const TOPO::IndexTransform &tr, int dst_axis) const;
    int coupling_face_2form_sign_src_to_dst_(const TOPO::IndexTransform &tr, int src_axis) const;
    int coupling_form_orientation_sign_dst_to_src_(const CouplingChannelSpec &ch,
                                                   const TOPO::IndexTransform &tr,
                                                   int dst_axis) const;
    int coupling_form_orientation_sign_src_to_dst_(const CouplingChannelSpec &ch,
                                                   const TOPO::IndexTransform &tr,
                                                   int src_axis) const;
    std::string find_triplet_field_name_(const std::string &dst_field_name,
                                         FieldValueKind value_kind,
                                         int wanted_src_axis) const;

private:
    Field *fld_;
    TOPO::Topology *topo_;
    const TOPO::Topology *equiv_ = nullptr;

    using HaloSyncSemantics = HALO_SYNC::Semantics;
    using HaloTripletRequest = HALO_SYNC::TripletRequest;
    using HaloOwnerRequest = HALO_SYNC::OwnerRequest;
    using OwnerSyncLocalOp = HALO_SYNC::OwnerLocalOp;
    using OwnerSyncSendOp = HALO_SYNC::OwnerSendOp;
    using OwnerSyncRecvOp = HALO_SYNC::OwnerRecvOp;
    using OwnerSyncPattern = HALO_SYNC::OwnerPattern;
    using CouplingPatternKey = HALO_SYNC::CouplingPatternKey;

    // 同一个 field_name 多次注册时取更高等级（Corner1D < Corner2D < Corner3D）
    std::unordered_map<std::string, FieldHaloRequest> halo_registry_;

    std::vector<std::string> component_copy_fields_;

    std::unordered_map<std::string, HaloTripletRequest> edge_1form_triplets_;
    std::unordered_map<std::string, HaloTripletRequest> face_2form_triplets_;

    std::vector<HaloOwnerRequest> owner_sync_requests_;

    std::unordered_map<std::string, std::string> field_to_group_;

    std::unordered_map<std::string, OwnerSyncPattern> owner_sync_patterns_;

    // Reused by owner synchronization; MPI requests are overwritten on every
    // call, while vector capacity remains available.
    std::vector<MPI_Request> owner_recv_requests_;
    std::vector<MPI_Request> owner_send_requests_;
    std::vector<MPI_Status> owner_recv_statuses_;
    std::vector<MPI_Status> owner_send_statuses_;

    int owner_sync_tag_base_ = 2100;

    // Sync registry helpers.
    HaloSyncSemantics sync_semantics_(const FieldHaloRequest &req) const;

    void upsert_halo_request_(const FieldHaloRequest &request);

    void rebuild_sync_registry_();

    void classify_registered_request_(const FieldHaloRequest &req);

    void validate_triplet_registry_() const;

    void validate_sync_registry_consistency_() const;

    const FieldHaloRequest &halo_request_(const std::string &field_name) const;

    void sync_component_copy_field_(const std::string &field_name);
    void sync_component_copy_field_stage_(const std::string &field_name, HaloLevel stage);

    void sync_component_copy_registered_();
    void sync_component_copy_registered_stage_(HaloLevel stage);

    void sync_edge_1form_triplet_(const HaloTripletRequest &tri);
    void sync_edge_1form_triplet_stage_(const HaloTripletRequest &tri, HaloLevel stage);

    void sync_edge_1form_triplets_registered_();
    void sync_edge_1form_triplets_registered_stage_(HaloLevel stage);

    void sync_face_2form_triplet_(const HaloTripletRequest &tri);
    void sync_face_2form_triplet_stage_(const HaloTripletRequest &tri, HaloLevel stage);

    void sync_face_2form_triplet_face_level_(const HaloTripletRequest &tri);

    void sync_face_2form_triplet_edge_level_(const HaloTripletRequest &tri);

    void sync_face_2form_triplet_vertex_level_(const HaloTripletRequest &tri);

    void sync_face_2form_triplets_registered_();
    void sync_face_2form_triplets_registered_stage_(HaloLevel stage);

    void sync_owner_alias_request_(const HaloOwnerRequest &req);

    void sync_owner_alias_registered_();
    bool sync_owner_alias_field_(const std::string &field_name);
    bool sync_owner_alias_group_(const std::string &group_name);

    bool field_is_component_copy_(const std::string &field_name) const;

    TOPO::EntityDim owner_policy_to_entity_dim_(OwnerSyncPolicy policy) const;

    void require_owner_equiv_available_(OwnerSyncPolicy policy,
                                        const std::string &field_name) const;

    int owner_alias_sign_(const HaloOwnerRequest &req,
                          const TOPO::EquivMember &owner,
                          const TOPO::EquivMember &alias) const;

    bool owner_member_matches_field_(const HaloOwnerRequest &req,
                                     const TOPO::EquivMember &m) const;

    int owner_sync_member_field_id_(const HaloOwnerRequest &req,
                                    const TOPO::EquivMember &member) const;

    void copy_owner_to_alias_local_(const HaloOwnerRequest &req,
                                    int fid,
                                    const TOPO::EquivMember &owner,
                                    const TOPO::EquivMember &alias,
                                    int sign);

    void build_owner_sync_patterns_();

    OwnerSyncPattern build_owner_sync_pattern_for_request_(const HaloOwnerRequest &req) const;

    int owner_sync_tag_(const HaloOwnerRequest &req) const;

    void resize_owner_sync_buffers_(OwnerSyncPattern &pat) const;

    void execute_owner_sync_local_ops_(const OwnerSyncPattern &pat);

    void pack_owner_sync_send_buffer_(OwnerSyncPattern &pat);

    void unpack_owner_sync_recv_buffer_(OwnerSyncPattern &pat);

    void execute_owner_sync_mpi_ops_(OwnerSyncPattern &pat);

    bool halo_level_includes_edge_(HaloLevel level) const;

    bool halo_level_includes_vertex_(HaloLevel level) const;

    // key = (StaggerLocation, nghost)，value = HaloPattern
    using PatternKey = std::pair<StaggerLocation, int>;
    std::map<PatternKey, HaloPattern> inner_patterns_;
    std::map<PatternKey, HaloPattern> parallel_patterns_;

    // For 2D Corner
    std::map<PatternKey, HaloPattern> inner_edge_patterns_;
    std::map<PatternKey, HaloPattern> parallel_edge_patterns_send;
    std::map<PatternKey, HaloPattern> parallel_edge_patterns_recv;

    // For 3D Corner
    std::map<PatternKey, HaloPattern> inner_vertex_patterns_;
    std::map<PatternKey, HaloPattern> parallel_vertex_patterns_send;
    std::map<PatternKey, HaloPattern> parallel_vertex_patterns_recv;

    // For Coupling (Parallel Corner): patterns are directed (src -> dst)
    std::map<CouplingPatternKey, HaloPattern> coupling_parallel_edge_patterns_send;
    std::map<CouplingPatternKey, HaloPattern> coupling_parallel_edge_patterns_recv;
    std::map<CouplingPatternKey, HaloPattern> coupling_parallel_vertex_patterns_send;
    std::map<CouplingPatternKey, HaloPattern> coupling_parallel_vertex_patterns_recv;

    // 复用的 send / recv 缓冲区（MPI 并行用）
    std::vector<std::vector<double>> send_buf;
    std::vector<std::vector<double>> recv_buf;
    std::vector<MPI_Request> req_send, req_recv;
    std::vector<MPI_Status> stat_send, stat_recv;
    std::vector<int32_t> length, length_corner_recv;

    //=========================================================================

    void build_inner_1DCorner_pattern(StaggerLocation loc, int nghost);
    void build_parallel_1DCorner_pattern(StaggerLocation loc, int nghost);
    void build_inner_2DCorner_pattern(StaggerLocation loc, int nghost);
    void build_parallel_2DCorner_pattern(StaggerLocation loc, int nghost);
    void build_inner_3DCorner_pattern(StaggerLocation loc, int nghost);
    void build_parallel_3DCorner_pattern(StaggerLocation loc, int nghost);

    // Coupling (Parallel) corner patterns (directed src -> dst)
    void build_coupling_parallel_2DCorner_pattern(const std::string &src,
                                                  const std::string &dst,
                                                  StaggerLocation loc,
                                                  int nghost);
    void build_coupling_parallel_3DCorner_pattern(const std::string &src,
                                                  const std::string &dst,
                                                  StaggerLocation loc,
                                                  int nghost);

    void mpi_exchange_edge_meta(
        const std::map<int, std::vector<EdgeMeta>> &meta_to_send,
        std::vector<EdgeMeta> &recv_metas);
    void mpi_exchange_vertex_meta(
        const std::map<int, std::vector<VertexMeta>> &meta_to_send,
        std::vector<VertexMeta> &recv_metas);

    bool field_active_(int fid, int iblock) const
    {
        return fld_->field(fid, iblock).is_allocated();
    }
};
