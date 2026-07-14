#pragma once

#include "6_io/PostBinaryWriter.h"

#include <cstdint>
#include <filesystem>
#include <map>
#include <string>
#include <vector>

class Field;
class Grid;
class RunData;
namespace TOPO { struct Topology; }
namespace METRIC { class SingularEdgeRegistry; }

namespace POST
{
struct WriteOptions
{
    Uuid case_uuid{};
    Uuid mesh_uuid{};
    std::vector<std::string> constant_fields;
    // Actual case normalization/scalar constants supplied by the caller.
    // Keeping these explicit avoids guessing solver-specific parameter names.
    std::map<std::string, double> normalization;
    std::map<std::string, double> physical_constants;
    std::vector<std::string> species;
    // Metadata only: describe the solver's existing checkpoint files so the
    // Python reader can reuse them. PostDataWriter never writes these files.
    std::string existing_flow_path_pattern{"./DATA/flow_field{rank:04d}.bin"};
    std::vector<std::string> existing_flow_fields;
    bool validate = true;
};

// Read-only adapter over the solver's existing Grid/Topology/Field objects.
// It never rebuilds topology or reconstruction weights. MPI output is chunked
// by rank and merged by global IDs by the serial reader.
class PostDataWriter
{
public:
    PostDataWriter(const Grid &grid, const TOPO::Topology &topology,
                   const Field &fields, const RunData &run_data, int rank,
                   const METRIC::SingularEdgeRegistry *singular_edges = nullptr);

    void WriteStaticData(const std::filesystem::path &output_directory,
                         WriteOptions options = {}) const;
private:
    const Grid &grid_;
    const TOPO::Topology &topology_;
    const Field &fields_;
    const RunData &run_data_;
    const METRIC::SingularEdgeRegistry *singular_edges_;
    int rank_ = 0;
};
} // namespace POST
