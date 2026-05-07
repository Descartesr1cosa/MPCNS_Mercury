#pragma once

#include <cstdint>
#include <cstring>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace VTKXML
{
    enum class DataType
    {
        Float64,
        Int32,
        Int64
    };

    struct AppendedArray
    {
        std::string name;
        DataType type = DataType::Float64;
        int number_of_components = 1;
        std::uint64_t offset = 0;
        std::uint64_t bytes = 0;
    };

    const char *DataTypeName(DataType type);

    class AppendedRawWriter
    {
    public:
        AppendedArray AddFloat64(const std::string &name,
                                 int number_of_components,
                                 const std::vector<double> &values)
        {
            return Add(name, DataType::Float64, number_of_components, values);
        }

        AppendedArray AddInt32(const std::string &name,
                               int number_of_components,
                               const std::vector<std::int32_t> &values)
        {
            return Add(name, DataType::Int32, number_of_components, values);
        }

        AppendedArray AddInt64(const std::string &name,
                               int number_of_components,
                               const std::vector<std::int64_t> &values)
        {
            return Add(name, DataType::Int64, number_of_components, values);
        }

        void WriteDataArrayTag(std::ostream &out,
                               const AppendedArray &array,
                               int indent_spaces) const;
        void WriteAppendedData(std::ostream &out) const;

    private:
        struct Block
        {
            AppendedArray array;
            std::vector<unsigned char> payload;
        };

        template <class T>
        AppendedArray Add(const std::string &name,
                          DataType type,
                          int number_of_components,
                          const std::vector<T> &values)
        {
            if (number_of_components <= 0)
                throw std::invalid_argument("[VTKXML] NumberOfComponents must be positive");

            Block block;
            block.array.name = name;
            block.array.type = type;
            block.array.number_of_components = number_of_components;
            block.array.offset = next_offset_;
            block.array.bytes = static_cast<std::uint64_t>(values.size() * sizeof(T));
            block.payload.resize(values.size() * sizeof(T));
            if (!values.empty())
                std::memcpy(block.payload.data(), values.data(), block.payload.size());

            next_offset_ += sizeof(std::uint64_t) + block.array.bytes;
            blocks_.push_back(block);
            return blocks_.back().array;
        }

        std::uint64_t next_offset_ = 0;
        std::vector<Block> blocks_;
    };
}
