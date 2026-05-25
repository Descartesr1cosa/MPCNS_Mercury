#pragma once

#include <map>
#include <stdexcept>
#include <vector>

#include "2_topology/Entity.h"

namespace TOPO
{
    struct IncidenceEntry
    {
        EntityKey entity;
        int sign;
    };

    inline std::vector<IncidenceEntry> boundary_of_edge(const EntityKey &edge)
    {
        if (edge.dim != EntityDim::Edge || !is_valid_entity_key_basic(edge))
        {
            throw std::invalid_argument("TOPO::boundary_of_edge: entity is not an edge.");
        }

        const EntityKey tail = make_node(edge.rank, edge.block, edge.i, edge.j, edge.k);
        EntityKey head = tail;
        switch (edge.axis)
        {
        case EntityAxis::Xi:
            ++head.i;
            break;
        case EntityAxis::Eta:
            ++head.j;
            break;
        case EntityAxis::Zeta:
            ++head.k;
            break;
        default:
            throw std::invalid_argument("TOPO::boundary_of_edge: invalid edge axis.");
        }

        return {{head, +1}, {tail, -1}};
    }

    inline std::vector<IncidenceEntry> boundary_of_face(const EntityKey &face)
    {
        if (face.dim != EntityDim::Face || !is_valid_entity_key_basic(face))
        {
            throw std::invalid_argument("TOPO::boundary_of_face: entity is not a face.");
        }

        const int r = face.rank;
        const int b = face.block;
        const int i = face.i;
        const int j = face.j;
        const int k = face.k;

        switch (face.axis)
        {
        case EntityAxis::Xi:
            return {
                {make_edge(r, b, i, j, k, EntityAxis::Eta), +1},
                {make_edge(r, b, i, j + 1, k, EntityAxis::Zeta), +1},
                {make_edge(r, b, i, j, k + 1, EntityAxis::Eta), -1},
                {make_edge(r, b, i, j, k, EntityAxis::Zeta), -1}};
        case EntityAxis::Eta:
            return {
                {make_edge(r, b, i, j, k, EntityAxis::Zeta), +1},
                {make_edge(r, b, i, j, k + 1, EntityAxis::Xi), +1},
                {make_edge(r, b, i + 1, j, k, EntityAxis::Zeta), -1},
                {make_edge(r, b, i, j, k, EntityAxis::Xi), -1}};
        case EntityAxis::Zeta:
            return {
                {make_edge(r, b, i, j, k, EntityAxis::Xi), +1},
                {make_edge(r, b, i + 1, j, k, EntityAxis::Eta), +1},
                {make_edge(r, b, i, j + 1, k, EntityAxis::Xi), -1},
                {make_edge(r, b, i, j, k, EntityAxis::Eta), -1}};
        default:
            throw std::invalid_argument("TOPO::boundary_of_face: invalid face axis.");
        }
    }

    inline std::vector<IncidenceEntry> boundary_of_cell(const EntityKey &cell)
    {
        if (cell.dim != EntityDim::Cell || !is_valid_entity_key_basic(cell))
        {
            throw std::invalid_argument("TOPO::boundary_of_cell: entity is not a cell.");
        }

        const int r = cell.rank;
        const int b = cell.block;
        const int i = cell.i;
        const int j = cell.j;
        const int k = cell.k;
        return {
            {make_face(r, b, i + 1, j, k, EntityAxis::Xi), +1},
            {make_face(r, b, i, j, k, EntityAxis::Xi), -1},
            {make_face(r, b, i, j + 1, k, EntityAxis::Eta), +1},
            {make_face(r, b, i, j, k, EntityAxis::Eta), -1},
            {make_face(r, b, i, j, k + 1, EntityAxis::Zeta), +1},
            {make_face(r, b, i, j, k, EntityAxis::Zeta), -1}};
    }

    inline bool check_boundary_boundary_face_zero(const EntityKey &face)
    {
        std::map<EntityKey, int> node_coefficients;
        for (const IncidenceEntry &edge_entry : boundary_of_face(face))
        {
            for (const IncidenceEntry &node_entry : boundary_of_edge(edge_entry.entity))
            {
                node_coefficients[node_entry.entity] += edge_entry.sign * node_entry.sign;
            }
        }

        for (const auto &entry : node_coefficients)
        {
            if (entry.second != 0)
            {
                return false;
            }
        }
        return true;
    }

    inline bool check_boundary_boundary_cell_zero(const EntityKey &cell)
    {
        std::map<EntityKey, int> edge_coefficients;
        for (const IncidenceEntry &face_entry : boundary_of_cell(cell))
        {
            for (const IncidenceEntry &edge_entry : boundary_of_face(face_entry.entity))
            {
                edge_coefficients[edge_entry.entity] += face_entry.sign * edge_entry.sign;
            }
        }

        for (const auto &entry : edge_coefficients)
        {
            if (entry.second != 0)
            {
                return false;
            }
        }
        return true;
    }
} // namespace TOPO
