#pragma once

#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <span>
#include <string_view>
#include <type_traits>

namespace POST
{
enum class FileType : std::uint32_t
{
    Geometry = 1,
    Topology = 2,
    Reconstruction = 3,
    ConstantField = 4
};

enum class ScalarType : std::uint32_t
{
    Int8 = 1,
    UInt8 = 2,
    Int32 = 3,
    UInt32 = 4,
    Int64 = 5,
    UInt64 = 6,
    Float64 = 7
};

struct Uuid
{
    std::uint64_t hi = 0;
    std::uint64_t lo = 0;
};

struct FileHeader
{
    std::array<char, 8> magic{{'M','P','C','N','S','B','I','N'}};
    std::uint32_t format_version = 1;
    std::uint32_t file_type = 0;
    std::uint64_t header_bytes = 80;
    std::uint64_t payload_bytes = 0;
    std::uint64_t case_uuid_hi = 0;
    std::uint64_t case_uuid_lo = 0;
    std::uint64_t mesh_uuid_hi = 0;
    std::uint64_t mesh_uuid_lo = 0;
    std::uint32_t endian_flag = 1;
    std::uint32_t float_bytes = 8;
    std::uint32_t index_bytes = 8;
    std::uint32_t reserved = 0;
};

class PostBinaryWriter
{
public:
    PostBinaryWriter(const std::filesystem::path &path, FileType type,
                     Uuid case_uuid, Uuid mesh_uuid);
    ~PostBinaryWriter() noexcept(false);

    PostBinaryWriter(const PostBinaryWriter &) = delete;
    PostBinaryWriter &operator=(const PostBinaryWriter &) = delete;

    template <class T>
    void WriteArray(std::span<const T> values)
    {
        static_assert(std::is_trivially_copyable_v<T>);
        if (!values.empty())
            WriteBytes(values.data(), values.size_bytes());
    }

    template <class T>
    void WriteSection(std::string_view name, std::uint32_t components,
                      std::span<const T> values)
    {
        WriteSectionHeader(name, ScalarCode<T>(), components,
                           components == 0 ? 0 : values.size() / components,
                           values.size_bytes());
        WriteArray(values);
    }

    void Close();
    std::uint64_t CurrentOffset();

private:
    template <class T> static constexpr ScalarType ScalarCode()
    {
        if constexpr (std::is_same_v<T, std::int8_t>) return ScalarType::Int8;
        else if constexpr (std::is_same_v<T, std::uint8_t>) return ScalarType::UInt8;
        else if constexpr (std::is_same_v<T, std::int32_t>) return ScalarType::Int32;
        else if constexpr (std::is_same_v<T, std::uint32_t>) return ScalarType::UInt32;
        else if constexpr (std::is_same_v<T, std::int64_t>) return ScalarType::Int64;
        else if constexpr (std::is_same_v<T, std::uint64_t>) return ScalarType::UInt64;
        else if constexpr (std::is_same_v<T, double>) return ScalarType::Float64;
        else static_assert(!sizeof(T), "unsupported post-data scalar type");
    }

    void WriteBytes(const void *data, std::size_t bytes);
    void WriteSectionHeader(std::string_view name, ScalarType type,
                            std::uint32_t components, std::uint64_t count,
                            std::uint64_t bytes);
    void RewriteHeader();

    std::filesystem::path path_;
    std::ofstream out_;
    FileHeader header_;
    bool closed_ = false;
};
} // namespace POST
