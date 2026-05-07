#include "5_io/IOModule.h"
#include "5_io/VTKXmlAppendedWriter.h"

#include <cstdio>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <vector>

namespace fs = std::filesystem;

namespace
{
    void WriteParaViewAppendedSmokeVTP(const fs::path &path)
    {
        VTKXML::AppendedRawWriter appended;

        const VTKXML::AppendedArray temperature =
            appended.AddFloat64("temperature", 1, std::vector<double>{300.0});
        const VTKXML::AppendedArray velocity =
            appended.AddFloat64("velocity", 3, std::vector<double>{1.0, 2.0, 3.0});
        const VTKXML::AppendedArray id32 =
            appended.AddInt32("id32", 1, std::vector<std::int32_t>{7});
        const VTKXML::AppendedArray id64 =
            appended.AddInt64("id64", 1, std::vector<std::int64_t>{7000000000LL});
        const VTKXML::AppendedArray cell_marker =
            appended.AddInt32("cell_marker", 1, std::vector<std::int32_t>{1});
        const VTKXML::AppendedArray points =
            appended.AddFloat64("", 3, std::vector<double>{0.0, 0.0, 0.0});
        const VTKXML::AppendedArray connectivity =
            appended.AddInt64("connectivity", 1, std::vector<std::int64_t>{0});
        const VTKXML::AppendedArray offsets =
            appended.AddInt64("offsets", 1, std::vector<std::int64_t>{1});

        std::ofstream out(path, std::ios::binary);
        if (!out)
            throw std::runtime_error("[IOModule][VTK] cannot open: " + path.string());

        out << "<?xml version=\"1.0\"?>\n"
            << "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n"
            << "  <PolyData>\n"
            << "    <Piece NumberOfPoints=\"1\" NumberOfVerts=\"1\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n"
            << "      <PointData Scalars=\"temperature\" Vectors=\"velocity\">\n";
        appended.WriteDataArrayTag(out, temperature, 8);
        appended.WriteDataArrayTag(out, velocity, 8);
        appended.WriteDataArrayTag(out, id32, 8);
        appended.WriteDataArrayTag(out, id64, 8);
        out << "      </PointData>\n"
            << "      <CellData Scalars=\"cell_marker\">\n";
        appended.WriteDataArrayTag(out, cell_marker, 8);
        out << "      </CellData>\n"
            << "      <Points>\n";
        appended.WriteDataArrayTag(out, points, 8);
        out << "      </Points>\n"
            << "      <Verts>\n";
        appended.WriteDataArrayTag(out, connectivity, 8);
        appended.WriteDataArrayTag(out, offsets, 8);
        out << "      </Verts>\n"
            << "      <Lines>\n"
            << "      </Lines>\n"
            << "      <Strips>\n"
            << "      </Strips>\n"
            << "      <Polys>\n"
            << "      </Polys>\n"
            << "    </Piece>\n"
            << "  </PolyData>\n";
        appended.WriteAppendedData(out);
        out << "</VTKFile>\n";
    }
}

void IOModule::SetParaViewFields(const std::vector<std::string> &names)
{
    paraview_fields_ = names;
}

void IOModule::ClearParaViewFields()
{
    paraview_fields_.clear();
}

void IOModule::SetParaViewPath(const std::string &path)
{
    paraview_path_ = path;
}

void IOModule::SetParaViewIncludeGhost(bool flag)
{
    paraview_include_ghost_ = flag;
}

void IOModule::WriteParaViewFile()
{
    const int rank = par_ ? par_->GetInt("myid") : 0;

    if (rank == 0)
    {
        fs::create_directories(paraview_path_);

        const fs::path output_path = fs::path(paraview_path_) / "paradata.vtm";
        const fs::path smoke_path = fs::path(paraview_path_) / "appended_smoke.vtp";

        try
        {
            WriteParaViewAppendedSmokeVTP(smoke_path);
        }
        catch (const std::exception &e)
        {
            Fail_(e.what());
        }

        std::ofstream out(output_path);
        if (!out)
            Fail_("[IOModule][VTK] cannot open: " + output_path.string());

        out << "<?xml version=\"1.0\"?>\n"
            << "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\">\n"
            << "  <vtkMultiBlockDataSet>\n"
            << "    <DataSet index=\"0\" name=\"appended_smoke\" file=\"appended_smoke.vtp\"/>\n"
            << "  </vtkMultiBlockDataSet>\n"
            << "</VTKFile>\n";

        std::printf("[IOModule][VTK] ParaView output is called: %s\n",
                    output_path.string().c_str());
        std::fflush(stdout);
    }

    if (par_)
        PARALLEL::mpi_barrier();
}
