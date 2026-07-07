#pragma once

#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "2_topology/Topology.h"
#include "2_topology/LocalIncidence.h"

namespace TOPO
{
    struct GlobalIncidenceEntry
    {
        EntityId entity;
        int sign;
    };

    // Metric-free incidence on quotient entities.  This view performs no
    // communication: it builds rows for the canonical entities represented
    // in the supplied, already-built Topology instance.
    class GlobalIncidence
    {
    public:
        explicit GlobalIncidence(const Topology &equiv)
            : equiv_(equiv)
        {
            build_edge_rows();
            build_face_rows();
            build_cell_rows();
        }

        GlobalIncidence(const Topology &equiv,
                        const std::vector<EntityKey> &local_cells)
            : GlobalIncidence(equiv)
        {
            for (const EntityKey &cell : local_cells)
                add_local_cell(cell);
        }

        std::vector<GlobalIncidenceEntry>
        boundary_of_global_edge(EntityId edge_id) const
        {
            require_dim(edge_id, EntityDim::Edge, "boundary_of_global_edge");
            return boundary_row(edge_rows_, edge_id, "boundary_of_global_edge");
        }

        std::vector<GlobalIncidenceEntry>
        boundary_of_global_face(EntityId face_id) const
        {
            require_dim(face_id, EntityDim::Face, "boundary_of_global_face");
            return boundary_row(face_rows_, face_id, "boundary_of_global_face");
        }

        std::vector<GlobalIncidenceEntry>
        boundary_of_global_cell(EntityId cell_id) const
        {
            require_dim(cell_id, EntityDim::Cell, "boundary_of_global_cell");
            return boundary_row(cell_rows_, cell_id, "boundary_of_global_cell");
        }

        void add_local_cell(const EntityKey &cell)
        {
            if (cell.dim != EntityDim::Cell)
                throw std::invalid_argument("GlobalIncidence::add_local_cell: entity is not a cell.");
            insert_row(cell_rows_, cell, boundary_of_cell(cell), "cell");
        }

        bool check_d1_d0_zero() const
        {
            for (const auto &[face_id, face_row] : face_rows_)
            {
                (void)face_id;
                Row node_coefficients;
                for (const auto &[edge_id, face_edge_sign] : face_row)
                {
                    const auto edge_it = edge_rows_.find(edge_id);
                    if (edge_it == edge_rows_.end())
                        return false;
                    for (const auto &[node_id, edge_node_sign] : edge_it->second)
                        add_term(node_coefficients, node_id, face_edge_sign * edge_node_sign);
                }
                if (!node_coefficients.empty())
                    return false;
            }
            return true;
        }

        bool check_d2_d1_zero() const
        {
            for (const auto &[cell_id, cell_row] : cell_rows_)
            {
                (void)cell_id;
                Row edge_coefficients;
                for (const auto &[face_id, cell_face_sign] : cell_row)
                {
                    const auto face_it = face_rows_.find(face_id);
                    if (face_it == face_rows_.end())
                        return false;
                    for (const auto &[edge_id, face_edge_sign] : face_it->second)
                        add_term(edge_coefficients, edge_id, cell_face_sign * face_edge_sign);
                }
                if (!edge_coefficients.empty())
                    return false;
            }
            return true;
        }

    private:
        using Row = std::map<EntityId, int>;
        using Rows = std::map<EntityId, Row>;

        const Topology &equiv_;
        Rows edge_rows_;
        Rows face_rows_;
        Rows cell_rows_;

        static void require_dim(EntityId id, EntityDim dim, const char *operation)
        {
            if (id.dim != dim)
                throw std::invalid_argument(std::string("GlobalIncidence::") + operation +
                                            ": unexpected entity dimension.");
        }

        static void add_term(Row &row, EntityId lower, int sign)
        {
            const int coefficient = (row[lower] += sign);
            if (coefficient == 0)
                row.erase(lower);
        }

        static std::vector<GlobalIncidenceEntry>
        boundary_row(const Rows &rows, EntityId upper, const char *operation)
        {
            const auto row_it = rows.find(upper);
            if (row_it == rows.end())
                throw std::runtime_error(std::string("GlobalIncidence::") + operation +
                                         ": quotient entity row is not available.");

            std::vector<GlobalIncidenceEntry> result;
            result.reserve(row_it->second.size());
            for (const auto &[entity, sign] : row_it->second)
                result.push_back(GlobalIncidenceEntry{entity, sign});
            return result;
        }

        void insert_row(Rows &rows,
                        const EntityKey &local_upper,
                        const std::vector<IncidenceEntry> &local_boundary,
                        const char *kind)
        {
            const EntityId upper_id = equiv_.id_of(local_upper);
            const int upper_sign = equiv_.sign_to_owner(local_upper);

            Row row;
            for (const IncidenceEntry &entry : local_boundary)
            {
                const EntityId lower_id = equiv_.id_of(entry.entity);
                const int lower_sign = equiv_.sign_to_owner(entry.entity);
                add_term(row, lower_id, upper_sign * entry.sign * lower_sign);
            }

            const auto existing = rows.find(upper_id);
            if (existing == rows.end())
            {
                rows.emplace(upper_id, std::move(row));
                return;
            }
            if (existing->second != row)
            {
                throw std::runtime_error(std::string("GlobalIncidence: inconsistent quotient ") +
                                         kind + " boundary across local representatives.");
            }
        }

        void build_edge_rows()
        {
            for (const auto &[edge, key] : equiv_.edges.local_to_qkey)
            {
                (void)key;
                const EntityKey local_edge = edge;
                if (!has_local_base_node(local_edge))
                    continue;
                insert_row(edge_rows_, local_edge, boundary_of_edge(local_edge), "edge");
            }
        }

        void build_face_rows()
        {
            for (const auto &[face, key] : equiv_.faces.local_to_qkey)
            {
                (void)key;
                const EntityKey local_face = face;
                if (!has_local_base_node(local_face))
                    continue;
                insert_row(face_rows_, local_face, boundary_of_face(local_face), "face");
            }
        }

        void build_cell_rows()
        {
            for (const auto &[cell, id] : equiv_.cells.local_to_qid)
            {
                (void)id;
                const EntityKey local_cell = cell;
                if (!has_local_base_node(local_cell))
                    continue;
                insert_row(cell_rows_, local_cell, boundary_of_cell(local_cell), "cell");
            }
        }

        bool has_local_base_node(const EntityKey &entity) const
        {
            const EntityKey base = make_node(entity.rank, entity.block, entity.i, entity.j, entity.k);
            return equiv_.nodes.local_to_rep.find(base) != equiv_.nodes.local_to_rep.end();
        }
    };
} // namespace TOPO
