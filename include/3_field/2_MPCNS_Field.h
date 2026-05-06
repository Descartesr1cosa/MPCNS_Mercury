#pragma once
#include <vector>
#include <string>

#include "3_field/FieldCatalog.h"
#include "3_field/FieldStorage.h"
#include "3_field/Coupling_Type.h"

class Field
{
public:
    // physic pair唯一表
    using PairKey = std::pair<std::string, std::string>; // (src,dst)

    Field() = default;
    ~Field() = default;

    Field(Grid *grd_, Param *par_, int geomtry_ghost_)
    {
        grd = grd_;
        par = par_;
        set_blocks(grd);
        build_geometry(geomtry_ghost_);
    };

    // 注册一个物理场（记录 desc，立刻分配）
    void register_field(const FieldDescriptor &desc);

    // 分配所有 fieldid 下 block 的数据
    void allocate(int32_t fieldID);

    int field_id(const std::string &field_name) const { return catalog_.field_id(field_name); }
    bool has_field(const std::string &name) const
    {
        return catalog_.has_field(name);
    }
    int num_fields() const { return catalog_.size(); }
    int num_blocks() const { return storage_.num_blocks(); }

    const FieldDescriptor &descriptor(int32_t fid) const { return catalog_.descriptor(fid); }
    const FieldDescriptor &descriptor(const std::string &field_name) const { return catalog_.descriptor(field_name); }
    const std::vector<FieldDescriptor> &descriptors() const { return catalog_.descriptors(); }

    std::vector<std::string> boundary_field_names() const;
    std::vector<std::string> coupled_field_names() const;
    std::vector<FieldHaloRequest> halo_requests() const;

    // 按 ID 访问所有block
    std::vector<FieldBlock> &field(int32_t fid)
    {
        return storage_.field(fid);
    }
    // 按 ID 访问
    FieldBlock &field(int32_t fid, int iblock)
    {
        return storage_.field(fid, iblock);
    }
    // 按名字访问所有block
    std::vector<FieldBlock> &field(const std::string &name)
    {
        return field(field_id(name));
    }
    // 按名字访问
    FieldBlock &field(const std::string &name, int iblock)
    {
        return field(field_id(name), iblock);
    }

    //===================================================================================
    void build_geometry(int geomtry_ghost_);
    //===================================================================================

    //===================================================================================
    // 为多物理场面耦合开辟缓冲区域
    //-----------------------------------------------------------------------------------
    void register_coupling_channel(const std::string &src,   // 源物理块
                                   const std::string &dst,   // 目标物理块
                                   const std::string &tag,   // 所需传输物理场的名称
                                   StaggerLocation location, // 物理场挂载的几何位置
                                   int ncomp,
                                   int nghost);
    void register_coupling_channel(const std::string &src,
                                   const std::string &dst,
                                   const std::string &field_name);
    void register_declared_coupling_channels(const std::vector<PairKey> &directed_pairs);
    //-----------------------------------------------------------------------------------
    // 为所注册的耦合方式（Channel）开辟缓冲空间，调用一次即可
    void build_coupling_buffers(const TOPO::Topology &topo, int dimension);
    //-----------------------------------------------------------------------------------
    // 查询（后面分配缓冲/Halo 会用）
    bool has_coupling_pair(const std::string &src, const std::string &dst) const;
    const CouplingPairDesc &coupling_pair(const std::string &src, const std::string &dst) const;
    CouplingBuffersForPair &coupling_buffers(const std::string &src, const std::string &dst);
    const std::map<PairKey, CouplingPairDesc> &coupling_pairs() const;
    //===================================================================================

private:
    // 存储网格指针
    void set_blocks(Grid *grd);

    FieldCatalog catalog_;
    FieldStorage storage_;

    std::map<PairKey, CouplingPairDesc> coupling_pairs_;

    std::map<PairKey, CouplingBuffersForPair> coupling_buffers_;

public:
    Grid *grd;
    Param *par;
};
