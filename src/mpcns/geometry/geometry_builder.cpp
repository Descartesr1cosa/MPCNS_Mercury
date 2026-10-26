#include "mpcns/geometry/geometry_builder.hpp"

#include <stdexcept>
#include <utility>

namespace mpcns::geometry
{

GeometryBuildResult GeometryBuilder::finalize(GeometryBuildInput input) const
{
    if (input.blocks.empty())
        throw std::invalid_argument("GeometryBuilder requires at least one geometry block");

    for (std::size_t block_index = 0; block_index < input.blocks.size(); ++block_index)
    {
        const GeometryBlock &block = input.blocks[block_index];
        if (block.name.empty())
            throw std::invalid_argument("GeometryBuilder received an unnamed geometry block");
        if (!block.is_complete())
            throw std::invalid_argument(
                "GeometryBuilder received an incomplete geometry block: " +
                block.name);
    }

    const std::size_t completed_block_count = input.blocks.size();
    return GeometryBuildResult{
        GeometrySnapshot{std::move(input.blocks)},
        completed_block_count};
}

} // namespace mpcns::geometry
