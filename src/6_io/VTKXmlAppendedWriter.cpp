#include "6_io/VTKXmlAppendedWriter.h"

#include <iomanip>

namespace VTKXML
{
    const char *DataTypeName(DataType type)
    {
        switch (type)
        {
        case DataType::Float64:
            return "Float64";
        case DataType::Int32:
            return "Int32";
        case DataType::Int64:
            return "Int64";
        }
        return "Float64";
    }

    void AppendedRawWriter::WriteDataArrayTag(std::ostream &out,
                                              const AppendedArray &array,
                                              int indent_spaces) const
    {
        out << std::string(static_cast<std::size_t>(indent_spaces), ' ')
            << "<DataArray type=\"" << DataTypeName(array.type) << "\"";
        if (!array.name.empty())
            out << " Name=\"" << array.name << "\"";
        out << " NumberOfComponents=\"" << array.number_of_components << "\""
            << " format=\"appended\" offset=\"" << array.offset << "\"/>\n";
    }

    void AppendedRawWriter::WriteAppendedData(std::ostream &out) const
    {
        out << "  <AppendedData encoding=\"raw\">\n_";
        for (const Block &block : blocks_)
        {
            const std::uint64_t bytes = block.array.bytes;
            out.write(reinterpret_cast<const char *>(&bytes), sizeof(bytes));
            if (!block.payload.empty())
            {
                out.write(reinterpret_cast<const char *>(block.payload.data()),
                          static_cast<std::streamsize>(block.payload.size()));
            }
        }
        out << "\n  </AppendedData>\n";
    }
}
