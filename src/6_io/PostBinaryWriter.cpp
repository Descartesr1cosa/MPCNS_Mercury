#include "6_io/PostBinaryWriter.h"

#include <bit>
#include <cstring>
#include <limits>
#include <stdexcept>

namespace POST
{
namespace
{
template <class T> void write_raw(std::ofstream &out, const T &value)
{
    out.write(reinterpret_cast<const char *>(&value), sizeof(value));
}

void write_header(std::ofstream &out, const FileHeader &h)
{
    out.write(h.magic.data(), h.magic.size());
    write_raw(out, h.format_version); write_raw(out, h.file_type);
    write_raw(out, h.header_bytes); write_raw(out, h.payload_bytes);
    write_raw(out, h.case_uuid_hi); write_raw(out, h.case_uuid_lo);
    write_raw(out, h.mesh_uuid_hi); write_raw(out, h.mesh_uuid_lo);
    write_raw(out, h.endian_flag); write_raw(out, h.float_bytes);
    write_raw(out, h.index_bytes); write_raw(out, h.reserved);
}
}

PostBinaryWriter::PostBinaryWriter(const std::filesystem::path &path, FileType type,
                                   Uuid case_uuid, Uuid mesh_uuid)
    : path_(path), out_(path, std::ios::binary | std::ios::trunc)
{
    static_assert(sizeof(double) == 8 && sizeof(std::int64_t) == 8);
    if constexpr (std::endian::native != std::endian::little)
        throw std::runtime_error("PostBinaryWriter: only little-endian hosts are supported by format v1");
    if (!out_) throw std::runtime_error("PostBinaryWriter: cannot open " + path.string());
    header_.file_type = static_cast<std::uint32_t>(type);
    header_.case_uuid_hi = case_uuid.hi; header_.case_uuid_lo = case_uuid.lo;
    header_.mesh_uuid_hi = mesh_uuid.hi; header_.mesh_uuid_lo = mesh_uuid.lo;
    write_header(out_, header_);
    if (!out_) throw std::runtime_error("PostBinaryWriter: cannot write header: " + path.string());
}

PostBinaryWriter::~PostBinaryWriter() noexcept(false)
{
    if (!closed_) Close();
}

void PostBinaryWriter::WriteBytes(const void *data, std::size_t bytes)
{
    if (closed_) throw std::logic_error("PostBinaryWriter: write after close");
    if (bytes > static_cast<std::size_t>(std::numeric_limits<std::streamsize>::max()))
        throw std::overflow_error("PostBinaryWriter: write exceeds streamsize");
    out_.write(static_cast<const char *>(data), static_cast<std::streamsize>(bytes));
    if (!out_) throw std::runtime_error("PostBinaryWriter: write failed: " + path_.string());
}

void PostBinaryWriter::WriteSectionHeader(std::string_view name, ScalarType type,
                                          std::uint32_t components, std::uint64_t count,
                                          std::uint64_t bytes)
{
    if (name.empty() || name.size() >= 32)
        throw std::invalid_argument("PostBinaryWriter: section name must contain 1..31 bytes");
    if (components == 0) throw std::invalid_argument("PostBinaryWriter: zero section components");
    std::array<char, 32> section_name{};
    std::memcpy(section_name.data(), name.data(), name.size());
    WriteBytes(section_name.data(), section_name.size());
    const std::uint32_t scalar = static_cast<std::uint32_t>(type);
    WriteBytes(&scalar, sizeof(scalar)); WriteBytes(&components, sizeof(components));
    WriteBytes(&count, sizeof(count)); WriteBytes(&bytes, sizeof(bytes));
}

std::uint64_t PostBinaryWriter::CurrentOffset()
{
    const auto p = out_.tellp();
    if (p < 0) throw std::runtime_error("PostBinaryWriter: tellp failed");
    return static_cast<std::uint64_t>(p);
}

void PostBinaryWriter::RewriteHeader()
{
    out_.seekp(0);
    write_header(out_, header_);
    if (!out_) throw std::runtime_error("PostBinaryWriter: header finalization failed: " + path_.string());
}

void PostBinaryWriter::Close()
{
    if (closed_) return;
    header_.payload_bytes = CurrentOffset() - header_.header_bytes;
    RewriteHeader();
    out_.flush();
    if (!out_) throw std::runtime_error("PostBinaryWriter: flush failed: " + path_.string());
    out_.close();
    if (out_.fail()) throw std::runtime_error("PostBinaryWriter: close failed: " + path_.string());
    closed_ = true;
}
} // namespace POST
