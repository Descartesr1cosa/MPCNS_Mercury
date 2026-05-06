#pragma once
#include "4_halo/Halo_Type.h"
#include "3_field/2_MPCNS_Field.h"
#include "0_basic/MPI_WRAPPER.h"

class Halo
{
public:
    Halo(Field *field, TOPO::Topology *topo)
    {
        fld_ = field;
        topo_ = topo;
    };

    // 记录哪个物理场需要 halo，做到哪个 corner 等级
    void register_halo_field(const FieldHaloRequest &request);
    void register_halo_field(const std::string &field_name,
                             HaloLevel level = HaloLevel::Vertex);
    void register_halo_fields(const std::vector<FieldHaloRequest> &requests)
    {
        for (const auto &request : requests)
            register_halo_field(request);
    }

    // 统一构建 pattern（只构建 registry 中需要的）
    void build_registered_patterns();

    void sync_registered();

    void sync_field(const std::string &field_name);

    void sync_group(const std::string &group_name);

    //=========================================================================
    // 普通面虚网格halo通信
    void data_trans_1DCorner(std::string &field_name)
    {
        exchange_inner(field_name);
        exchange_parallel(field_name);
    }

    // 2D棱虚网格halo通信
    void data_trans_2DCorner(std::string &field_name)
    {
        exchange_inner_edge(field_name);
        exchange_parallel_edge(field_name);
    }

    // 3D角虚网格halo通信
    void data_trans_3DCorner(std::string &field_name)
    {
        exchange_inner_vertex(field_name);
        exchange_parallel_vertex(field_name);
    }

    // Orientation-aware triplet halo for edge 1-forms such as E/J/dE.
    // fields must be ordered {xi, eta, zeta}. Across block transforms, the
    // source component is selected by IndexTransform::perm and multiplied by
    // IndexTransform::sign.
    void data_trans_edge_1form_triplet(const std::vector<std::string> &fields,
                                       HaloLevel stage);
    //=========================================================================

    //=========================================================================
    // 普通面虚网格耦合通信（fill coupling buffer, no ghost write）
    void coupling_trans_1DCorner(std::string &src, std::string &dst)
    {
        coupling_inner_face(src, dst);
        coupling_parallel_face(src, dst);
    }

    // fill coupling buffer, no ghost write, only for field: field_name -> field_cids
    // field_cids: corresponding channel ids for field_name(s)
    void coupling_trans_1DCorner(std::string &src, std::string &dst, std::vector<int32_t> &field_cids)
    {
        coupling_inner_face(src, dst, field_cids);
        coupling_parallel_face(src, dst, field_cids);
    }

    // 2D棱虚网格耦合通信（fill coupling buffer, no ghost write）
    void coupling_trans_2DCorner(std::string &src, std::string &dst)
    {
        coupling_inner_edge(src, dst);
        coupling_parallel_edge(src, dst);
    }

    // fill coupling buffer, no ghost write, only for field: field_name -> field_cids
    // field_cids: corresponding channel ids for field_name(s)
    void coupling_trans_2DCorner(std::string &src, std::string &dst, std::vector<int32_t> &field_cids)
    {
        coupling_inner_edge(src, dst, field_cids);
        coupling_parallel_edge(src, dst, field_cids);
    }

    // 3D角虚网格耦合通信（fill coupling buffer, no ghost write）
    void coupling_trans_3DCorner(std::string &src, std::string &dst)
    {
        coupling_inner_vertex(src, dst);
        coupling_parallel_vertex(src, dst);
    }

    // fill coupling buffer, no ghost write, only for field: field_name -> field_cids
    // field_cids: corresponding channel ids for field_name(s)
    void coupling_trans_3DCorner(std::string &src, std::string &dst, std::vector<int32_t> &field_cids)
    {
        coupling_inner_vertex(src, dst, field_cids);
        coupling_parallel_vertex(src, dst, field_cids);
    }
    //=========================================================================

private:
    //=========================================================================
    // face
    void exchange_inner(std::string field_name);
    void exchange_parallel(std::string field_name);

    // edge
    void exchange_inner_edge(std::string field_name);
    void exchange_parallel_edge(std::string field_name);

    // vertex
    void exchange_inner_vertex(std::string field_name);
    void exchange_parallel_vertex(std::string field_name);

    void exchange_inner_face_edge_1form_triplet_(const std::vector<std::string> &fields);
    void exchange_parallel_face_edge_1form_triplet_(const std::vector<std::string> &fields);
    void exchange_inner_edge_edge_1form_triplet_(const std::vector<std::string> &fields);
    void exchange_parallel_edge_edge_1form_triplet_(const std::vector<std::string> &fields);
    void exchange_inner_vertex_edge_1form_triplet_(const std::vector<std::string> &fields);
    void exchange_parallel_vertex_edge_1form_triplet_(const std::vector<std::string> &fields);
    //=========================================================================

    //=========================================================================
    // coupling face
    void coupling_inner_face(const std::string &src, const std::string &dst);
    void coupling_parallel_face(const std::string &src, const std::string &dst);
    void coupling_inner_face(const std::string &src, const std::string &dst, std::vector<int32_t> &field_cids);
    void coupling_parallel_face(const std::string &src, const std::string &dst, std::vector<int32_t> &field_cids);

    // coupling edge
    void coupling_inner_edge(const std::string &src, const std::string &dst);
    void coupling_parallel_edge(const std::string &src, const std::string &dst);
    void coupling_inner_edge(const std::string &src, const std::string &dst, std::vector<int32_t> &field_cids);
    void coupling_parallel_edge(const std::string &src, const std::string &dst, std::vector<int32_t> &field_cids);

    // coupling vertex
    void coupling_inner_vertex(const std::string &src, const std::string &dst);
    void coupling_parallel_vertex(const std::string &src, const std::string &dst);
    void coupling_inner_vertex(const std::string &src, const std::string &dst, std::vector<int32_t> &field_cids);
    void coupling_parallel_vertex(const std::string &src, const std::string &dst, std::vector<int32_t> &field_cids);
    //=========================================================================

private:
    Field *fld_;
    TOPO::Topology *topo_;

    enum class HaloSyncSemantics
    {
        ComponentCopy,
        Edge1FormTriplet,
        Face2FormTriplet
    };

    struct HaloTripletRequest
    {
        std::string group_name;

        std::string xi;
        std::string eta;
        std::string zeta;

        FieldValueKind value_kind = FieldValueKind::Scalar;

        HaloLevel level = HaloLevel::Vertex;
        int nghost = 0;

        bool orientation_aware = true;
    };

    struct HaloOwnerRequest
    {
        std::string field_name;
        std::string sync_group;

        OwnerSyncPolicy policy = OwnerSyncPolicy::None;

        FieldValueKind value_kind = FieldValueKind::Scalar;
        StaggerLocation location = StaggerLocation::Cell;

        bool orientation_aware = false;
    };

    // 同一个 field_name 多次注册时取更高等级（FaceOnly < Edge < Vertex）
    std::unordered_map<std::string, FieldHaloRequest> halo_registry_;

    std::vector<std::string> component_copy_fields_;

    std::unordered_map<std::string, HaloTripletRequest> edge_1form_triplets_;
    std::unordered_map<std::string, HaloTripletRequest> face_2form_triplets_;

    std::vector<HaloOwnerRequest> owner_sync_requests_;

    std::unordered_map<std::string, std::string> field_to_group_;

    HaloSyncSemantics sync_semantics_(const FieldHaloRequest &req) const;

    void rebuild_sync_registry_();

    void classify_registered_request_(const FieldHaloRequest &req);

    void validate_triplet_registry_() const;

    const FieldHaloRequest &halo_request_(const std::string &field_name) const;

    void sync_component_copy_field_(const std::string &field_name);

    void sync_component_copy_registered_();

    void sync_edge_1form_triplet_(const HaloTripletRequest &tri);

    void sync_edge_1form_triplets_registered_();

    void sync_face_2form_triplet_(const HaloTripletRequest &tri);

    void sync_face_2form_triplets_registered_();

    void sync_owner_alias_request_(const HaloOwnerRequest &req);

    void sync_owner_alias_registered_();

    bool field_is_component_copy_(const std::string &field_name) const;

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

    // Coupling (src->dst) parallel corner patterns need their own cache,
    // because the pattern depends on the directed physical pair in addition
    // to (loc, nghost).
    struct CouplingPatternKey
    {
        std::string src;
        std::string dst;
        StaggerLocation loc;
        int nghost;

        bool operator<(const CouplingPatternKey &o) const
        {
            if (src != o.src)
                return src < o.src;
            if (dst != o.dst)
                return dst < o.dst;
            if ((int)loc != (int)o.loc)
                return (int)loc < (int)o.loc;
            return nghost < o.nghost;
        }
    };

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
