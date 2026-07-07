#pragma once

#include <cstdint>
#include <unordered_map>
#include <vector>

#include "2_topology/Topology.h"

class Grid;

namespace TOPO::detail
{
    struct EdgeBuildScratch
    {
        std::unordered_map<EdgeKey, std::vector<EntityKey>, EdgeKey::Hash> qkey_to_local_members;
    };

    struct FaceBuildScratch
    {
        std::unordered_map<FaceKey, std::vector<EntityKey>, FaceKey::Hash> qkey_to_local_members;
    };

    void pack_node(std::vector<int> &buf, const EntityKey &x);
    EntityKey unpack_node(const int *p);
    std::vector<int> allgather_packed_records(std::vector<int> send_buf,
                                              int record_size,
                                              const char *context);

    void pack_edge_local(std::vector<int> &buf, const EntityKey &e);
    EntityKey unpack_edge_local(const int *p);
    void pack_edge_key(std::vector<int> &buf, const EdgeKey &key);
    EdgeKey unpack_edge_key(const int *p);

    void pack_face_local(std::vector<int> &buf, const EntityKey &f);
    EntityKey unpack_face_local(const int *p);
    void pack_face_key(std::vector<int> &buf, const FaceKey &key);
    FaceKey unpack_face_key(const int *p);

    std::vector<EntityKey> collect_all_local_nodes(Grid &grid, int my_rank);
    std::vector<EntityKey> collect_all_local_edges(Grid &grid, int my_rank, int dimension);
    std::vector<EntityKey> collect_all_local_faces(Grid &grid, int my_rank, int dimension);
    std::vector<EntityKey> collect_all_local_cells(Grid &grid, int my_rank, int dimension);

    void reconcile_node_equivalence_parallel(
        const Topology &topo,
        int my_rank,
        const std::vector<EntityKey> &all_local_nodes,
        Topology &equiv,
        std::unordered_map<EntityKey, int, EntityKey::Hash> &rep_count);
    void build_node_entity_ids(Topology &equiv);

    void build_edge_equivalence_from_nodes(
        const std::vector<EntityKey> &all_local_edges,
        Topology &equiv,
        EdgeBuildScratch &scratch);
    void build_edge_entity_ids(Topology &equiv);
    void select_edge_owner_parallel(
        Topology &equiv,
        const EdgeBuildScratch &scratch,
        const std::unordered_map<EntityKey, int, EntityKey::Hash> &rep_count);
    void build_edge_owner_gid(int my_rank, Topology &equiv);

    void build_face_equivalence_from_nodes(
        const std::vector<EntityKey> &all_local_faces,
        Topology &equiv,
        FaceBuildScratch &scratch);
    void build_face_entity_ids(Topology &equiv);
    void select_face_owner_parallel(
        Topology &equiv,
        const FaceBuildScratch &scratch,
        const std::unordered_map<EntityKey, int, EntityKey::Hash> &rep_count);
    void build_face_owner_gid(int my_rank, Topology &equiv);

    void build_cell_entity_ids(
        const std::vector<EntityKey> &all_local_cells,
        Topology &equiv);
} // namespace TOPO::detail
