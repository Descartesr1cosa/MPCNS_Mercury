#pragma once
#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include "3_field/FieldCatalog.h"
#include "3_field/FieldStorage.h"
#include "3_field/CouplingTypes.h"

namespace TOPO
{
    struct Topology;
}

class Field
{
public:
    // Directed physics pair: source block name -> destination block name.
    using PairKey = std::pair<std::string, std::string>;

    Field() = default;
    ~Field() = default;
    Field(Grid *grid, Param *param, int geometry_ghost);

    // Metadata registration. Storage is allocated immediately after a grid is bound.
    void register_field(const FieldDescriptor &desc);
    void allocate(int32_t fieldID);

    // Catalog queries.
    int field_id(const std::string &field_name) const;
    bool has_field(const std::string &name) const;
    int num_fields() const { return catalog_.size(); }
    int num_blocks() const { return storage_.num_blocks(); }

    const FieldDescriptor &descriptor(int32_t fid) const;
    const FieldDescriptor &descriptor(const std::string &field_name) const;
    const std::vector<FieldDescriptor> &descriptors() const;

    std::vector<std::string> boundary_field_names() const;
    std::vector<std::string> coupled_field_names() const;
    std::vector<FieldHaloRequest> halo_requests() const;

    // Field data access.
    std::vector<FieldBlock> &field(int32_t fid);
    const std::vector<FieldBlock> &field(int32_t fid) const;
    FieldBlock &field(int32_t fid, int iblock);
    const FieldBlock &field(int32_t fid, int iblock) const;
    std::vector<FieldBlock> &field(const std::string &name);
    const std::vector<FieldBlock> &field(const std::string &name) const;
    FieldBlock &field(const std::string &name, int iblock);
    const FieldBlock &field(const std::string &name, int iblock) const;

    // Geometry fields.
    void build_geometry(int geomtry_ghost_);

    // Coupling channel registration.
    void register_coupling_channel(const std::string &src,
                                   const std::string &dst,
                                   const std::string &tag,
                                   StaggerLocation location,
                                   FieldValueKind value_kind,
                                   int ncomp,
                                   int nghost,
                                   bool orientation_aware);
    void register_coupling_channel(const std::string &src,
                                   const std::string &dst,
                                   const std::string &tag,
                                   StaggerLocation location,
                                   int ncomp,
                                   int nghost);
    void register_coupling_channel(const std::string &src,
                                   const std::string &dst,
                                   const std::string &field_name);
    void register_declared_coupling_channels(const std::vector<PairKey> &directed_pairs);

    // Coupling buffers. Call after topology and coupling channels are ready.
    void build_coupling_buffers(const TOPO::Topology &topo, int dimension);
    bool has_coupling_pair(const std::string &src, const std::string &dst) const;
    const CouplingPairDesc &coupling_pair(const std::string &src, const std::string &dst) const;
    CouplingBuffersForPair &coupling_buffers(const std::string &src, const std::string &dst);
    const std::map<PairKey, CouplingPairDesc> &coupling_pairs() const;

private:
    void bind_grid(Grid *grid);

    FieldCatalog catalog_;
    FieldStorage storage_;
    std::map<PairKey, CouplingPairDesc> coupling_pairs_;
    std::map<PairKey, CouplingBuffersForPair> coupling_buffers_;

public:
    // Legacy public handles kept for existing solver code.
    Grid *grd = nullptr;
    Param *par = nullptr;
};
