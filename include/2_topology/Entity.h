#pragma once

#include <compare>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <stdexcept>

namespace TOPO
{
    // Topological dimension of an entity in the discrete complex.
    enum class EntityDim
    {
        Node = 0,
        Edge = 1,
        Face = 2,
        Cell = 3
    };

    // Reference coordinate axis for oriented edge and face entities.
    enum class EntityAxis
    {
        None = -1,
        Xi = 0,
        Eta = 1,
        Zeta = 2
    };

    // Identity of one real entity in a local rank/block topology.
    // Indices are in node space; ordinary ghost storage is not represented.
    struct EntityKey
    {
        EntityDim dim;
        int rank;
        int block;
        int i;
        int j;
        int k;
        EntityAxis axis;

        auto operator<=>(const EntityKey &) const = default;

        struct Hash
        {
            std::size_t operator()(const EntityKey &x) const
            {
                std::size_t h = 0;
                const auto combine = [&h](std::size_t value)
                {
                    h ^= value + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
                };
                combine(std::hash<int>{}(static_cast<int>(x.dim)));
                combine(std::hash<int>{}(x.rank));
                combine(std::hash<int>{}(x.block));
                combine(std::hash<int>{}(x.i));
                combine(std::hash<int>{}(x.j));
                combine(std::hash<int>{}(x.k));
                combine(std::hash<int>{}(static_cast<int>(x.axis)));
                return h;
            }
        };
    };

    // Per-block node-space bounds.  Each value is the maximum node index.
    struct BlockTopoSize
    {
        int imax;
        int jmax;
        int kmax;
    };

    // Identity of one entity after quotienting interfaces and duplicates.
    // The integer id space is separate for each entity dimension.
    struct EntityId
    {
        EntityDim dim;
        int id;

        auto operator<=>(const EntityId &) const = default;

        struct Hash
        {
            std::size_t operator()(const EntityId &x) const
            {
                std::size_t h = 0;
                const auto combine = [&h](std::size_t value)
                {
                    h ^= value + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
                };
                combine(std::hash<int>{}(static_cast<int>(x.dim)));
                combine(std::hash<int>{}(x.id));
                return h;
            }
        };
    };

    enum class EntityActivity
    {
        Active,
        Inactive
    };

    enum class EntityFlag : std::uint32_t
    {
        PhysicalBoundary = 1u << 0,
        BlockInterface = 1u << 1,
        PeriodicInterface = 1u << 2,
        Axis = 1u << 3,
        Collapsed = 1u << 4,
        Singular = 1u << 5
    };

    using EntityFlags = std::uint32_t;

    constexpr EntityFlags entity_flag(EntityFlag value)
    {
        return static_cast<EntityFlags>(value);
    }

    constexpr EntityFlags operator|(EntityFlag lhs, EntityFlag rhs)
    {
        return entity_flag(lhs) | entity_flag(rhs);
    }

    constexpr EntityFlags operator|(EntityFlags lhs, EntityFlag rhs)
    {
        return lhs | entity_flag(rhs);
    }

    constexpr bool has_flag(EntityFlags flags, EntityFlag flag)
    {
        return (flags & entity_flag(flag)) != 0u;
    }

    struct EntityMetadata
    {
        EntityFlags flags = 0u;
        EntityActivity activity = EntityActivity::Active;
        int boundary_marker = -1;
    };

    inline bool is_oriented_axis(EntityAxis axis)
    {
        return axis == EntityAxis::Xi ||
               axis == EntityAxis::Eta ||
               axis == EntityAxis::Zeta;
    }

    inline bool valid_axis_for_dim(EntityDim dim, EntityAxis axis)
    {
        switch (dim)
        {
        case EntityDim::Node:
        case EntityDim::Cell:
            return axis == EntityAxis::None;
        case EntityDim::Edge:
        case EntityDim::Face:
            return is_oriented_axis(axis);
        default:
            return false;
        }
    }

    inline bool is_valid_entity_key_basic(const EntityKey &key)
    {
        return key.rank >= 0 &&
               key.block >= 0 &&
               key.i >= 0 &&
               key.j >= 0 &&
               key.k >= 0 &&
               valid_axis_for_dim(key.dim, key.axis);
    }

    inline bool is_valid_entity_key(const EntityKey &key, const BlockTopoSize &size)
    {
        if (!is_valid_entity_key_basic(key) ||
            size.imax < 0 || size.jmax < 0 || size.kmax < 0)
        {
            return false;
        }

        const auto at_most = [](int value, int maximum)
        {
            return value <= maximum;
        };

        switch (key.dim)
        {
        case EntityDim::Node:
            return at_most(key.i, size.imax) &&
                   at_most(key.j, size.jmax) &&
                   at_most(key.k, size.kmax);
        case EntityDim::Edge:
            switch (key.axis)
            {
            case EntityAxis::Xi:
                return at_most(key.i, size.imax - 1) &&
                       at_most(key.j, size.jmax) &&
                       at_most(key.k, size.kmax);
            case EntityAxis::Eta:
                return at_most(key.i, size.imax) &&
                       at_most(key.j, size.jmax - 1) &&
                       at_most(key.k, size.kmax);
            case EntityAxis::Zeta:
                return at_most(key.i, size.imax) &&
                       at_most(key.j, size.jmax) &&
                       at_most(key.k, size.kmax - 1);
            default:
                return false;
            }
        case EntityDim::Face:
            switch (key.axis)
            {
            case EntityAxis::Xi:
                return at_most(key.i, size.imax) &&
                       at_most(key.j, size.jmax - 1) &&
                       at_most(key.k, size.kmax - 1);
            case EntityAxis::Eta:
                return at_most(key.i, size.imax - 1) &&
                       at_most(key.j, size.jmax) &&
                       at_most(key.k, size.kmax - 1);
            case EntityAxis::Zeta:
                return at_most(key.i, size.imax - 1) &&
                       at_most(key.j, size.jmax - 1) &&
                       at_most(key.k, size.kmax);
            default:
                return false;
            }
        case EntityDim::Cell:
            return at_most(key.i, size.imax - 1) &&
                   at_most(key.j, size.jmax - 1) &&
                   at_most(key.k, size.kmax - 1);
        default:
            return false;
        }
    }

    inline EntityKey make_node(int rank, int block, int i, int j, int k)
    {
        return EntityKey{EntityDim::Node, rank, block, i, j, k, EntityAxis::None};
    }

    inline EntityKey make_edge(int rank, int block, int i, int j, int k, EntityAxis axis)
    {
        if (!is_oriented_axis(axis))
        {
            throw std::invalid_argument("TOPO::make_edge: edge axis must be Xi, Eta, or Zeta.");
        }
        return EntityKey{EntityDim::Edge, rank, block, i, j, k, axis};
    }

    inline EntityKey make_face(int rank, int block, int i, int j, int k, EntityAxis axis)
    {
        if (!is_oriented_axis(axis))
        {
            throw std::invalid_argument("TOPO::make_face: face axis must be Xi, Eta, or Zeta.");
        }
        return EntityKey{EntityDim::Face, rank, block, i, j, k, axis};
    }

    inline EntityKey make_cell(int rank, int block, int i, int j, int k)
    {
        return EntityKey{EntityDim::Cell, rank, block, i, j, k, EntityAxis::None};
    }

} // namespace TOPO
