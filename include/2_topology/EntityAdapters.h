#pragma once

#include <stdexcept>

#include "2_topology/Entity.h"
#include "2_topology/TopologyEquiv.h"

namespace TOPO
{
    inline EntityAxis to_entity_axis(int dir)
    {
        switch (dir)
        {
        case 1:
            return EntityAxis::Xi;
        case 2:
            return EntityAxis::Eta;
        case 3:
            return EntityAxis::Zeta;
        default:
            throw std::invalid_argument("TOPO::to_entity_axis: legacy direction must be 1, 2, or 3.");
        }
    }

    inline int to_legacy_dir(EntityAxis axis)
    {
        switch (axis)
        {
        case EntityAxis::Xi:
            return 1;
        case EntityAxis::Eta:
            return 2;
        case EntityAxis::Zeta:
            return 3;
        default:
            throw std::invalid_argument("TOPO::to_legacy_dir: entity axis must be Xi, Eta, or Zeta.");
        }
    }

    inline EntityKey to_entity_key(const LocalNodeID &node)
    {
        return make_node(node.rank, node.gblock, node.i, node.j, node.k);
    }

    inline EntityKey to_entity_key(const EdgeLocalID &edge)
    {
        return make_edge(edge.rank, edge.gblock, edge.i, edge.j, edge.k,
                         to_entity_axis(edge.dir));
    }

    inline EntityKey to_entity_key(const FaceLocalID &face)
    {
        return make_face(face.rank, face.gblock, face.i, face.j, face.k,
                         to_entity_axis(face.dir));
    }

    inline LocalNodeID to_local_node_id(const EntityKey &entity)
    {
        if (entity.dim != EntityDim::Node || entity.axis != EntityAxis::None)
        {
            throw std::invalid_argument("TOPO::to_local_node_id: EntityKey is not a node.");
        }
        return LocalNodeID{entity.rank, entity.block, entity.i, entity.j, entity.k};
    }

    inline EdgeLocalID to_edge_local_id(const EntityKey &entity)
    {
        if (entity.dim != EntityDim::Edge)
        {
            throw std::invalid_argument("TOPO::to_edge_local_id: EntityKey is not an edge.");
        }
        return EdgeLocalID{entity.rank, entity.block, entity.i, entity.j, entity.k,
                           to_legacy_dir(entity.axis)};
    }

    inline FaceLocalID to_face_local_id(const EntityKey &entity)
    {
        if (entity.dim != EntityDim::Face)
        {
            throw std::invalid_argument("TOPO::to_face_local_id: EntityKey is not a face.");
        }
        return FaceLocalID{entity.rank, entity.block, entity.i, entity.j, entity.k,
                           to_legacy_dir(entity.axis)};
    }
} // namespace TOPO
