#pragma once
#include <vector>
#include <string>
#include <unordered_map>

#include "1_grid/1_MPCNS_Grid.h"   // Block
#include "3_field/1_Field_Block.h" // FieldBlock
#include "3_field/Field_Type.h"    // FieldDescriptor
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

    int field_id(std::string field_name) { return name_to_id_[field_name]; }
    bool has_field(const std::string &name) const
    {
        return name_to_id_.find(name) != name_to_id_.end();
    }
    int num_fields() const { return static_cast<int>(field_descs_.size()); }
    int num_blocks() const { return static_cast<int>(blocks_.size()); }

    const FieldDescriptor &descriptor(int32_t fid) const { return field_descs_[fid]; }

    // 按 ID 访问所有block
    std::vector<FieldBlock> &field(int32_t fid)
    {
        return field_blocks_[fid];
    }
    // 按 ID 访问
    FieldBlock &field(int32_t fid, int iblock)
    {
        return field_blocks_[fid][iblock];
    }
    // 按名字访问所有block
    std::vector<FieldBlock> &field(std::string name)
    {
        return field(name_to_id_.at(name));
    }
    // 按名字访问
    FieldBlock &field(std::string name, int iblock)
    {
        return field(name_to_id_.at(name), iblock);
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

    // 这个 rank 上的所有 Block（只存指针，不拥有）
    std::vector<Block *> blocks_;

    // 所有场的描述
    std::vector<FieldDescriptor> field_descs_;
    std::unordered_map<std::string, int32_t> name_to_id_;

    std::unordered_map<std::string, std::vector<int>> blocks_by_name_;

    // 真正的数据：field_blocks_[fid][iblock]
    std::vector<std::vector<FieldBlock>> field_blocks_;

    std::map<PairKey, CouplingPairDesc> coupling_pairs_;

    std::map<PairKey, CouplingBuffersForPair> coupling_buffers_;

public:
    Grid *grd;
    Param *par;
};