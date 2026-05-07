#include "5_io/IOModule.h"
#include "5_io/VTKXmlAppendedWriter.h"

#include <array>
#include <cstdio>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace
{
    struct ParaViewDatasetSpec
    {
        const char *name;
        const char *suffix;
        const char *extension;
    };

    const std::array<ParaViewDatasetSpec, 7> kParaViewDatasets = {
        {{"volume", "volume", ".vts"},
         {"edge_xi", "edge_xi", ".vtp"},
         {"edge_eta", "edge_eta", ".vtp"},
         {"edge_zeta", "edge_zeta", ".vtp"},
         {"face_xi", "face_xi", ".vtp"},
         {"face_eta", "face_eta", ".vtp"},
         {"face_zeta", "face_zeta", ".vtp"}}};

    std::string RankBlockPrefix(int rank, int block)
    {
        std::ostringstream os;
        os << "rank" << std::setw(4) << std::setfill('0') << rank
           << "_block" << std::setw(4) << std::setfill('0') << block;
        return os.str();
    }

    std::string DatasetFilename(const std::string &prefix, const ParaViewDatasetSpec &spec)
    {
        return prefix + "_" + spec.suffix + spec.extension;
    }

    void WritePlaceholderVTS(const fs::path &path)
    {
        VTKXML::AppendedRawWriter appended;
        const VTKXML::AppendedArray points =
            appended.AddFloat64("", 3, std::vector<double>{0.0, 0.0, 0.0});

        std::ofstream out(path, std::ios::binary);
        if (!out)
            throw std::runtime_error("[IOModule][VTK] cannot open: " + path.string());

        out << "<?xml version=\"1.0\"?>\n"
            << "<VTKFile type=\"StructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n"
            << "  <StructuredGrid WholeExtent=\"0 0 0 0 0 0\">\n"
            << "    <Piece Extent=\"0 0 0 0 0 0\">\n"
            << "      <PointData>\n"
            << "      </PointData>\n"
            << "      <CellData>\n"
            << "      </CellData>\n"
            << "      <Points>\n";
        appended.WriteDataArrayTag(out, points, 8);
        out << "      </Points>\n"
            << "    </Piece>\n"
            << "  </StructuredGrid>\n";
        appended.WriteAppendedData(out);
        out << "</VTKFile>\n";
    }

    void WritePlaceholderVTP(const fs::path &path)
    {
        VTKXML::AppendedRawWriter appended;
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
            << "      <PointData>\n"
            << "      </PointData>\n"
            << "      <CellData>\n"
            << "      </CellData>\n"
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

    void WritePlaceholderDataset(const fs::path &path, const ParaViewDatasetSpec &spec)
    {
        const std::string extension = spec.extension;
        if (extension == ".vts")
            WritePlaceholderVTS(path);
        else if (extension == ".vtp")
            WritePlaceholderVTP(path);
        else
            throw std::runtime_error("[IOModule][VTK] unsupported placeholder extension: " + extension);
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
    if (!par_ || !grd_)
        Fail_("[IOModule][VTK] Setup() must be called before WriteParaViewFile()");

    const int rank = par_->GetInt("myid");

    try
    {
        fs::create_directories(paraview_path_);

        const int nblock = grd_->nblock;
        for (int ib = 0; ib < nblock; ++ib)
        {
            const std::string prefix = RankBlockPrefix(rank, ib);
            for (const ParaViewDatasetSpec &spec : kParaViewDatasets)
                WritePlaceholderDataset(fs::path(paraview_path_) / DatasetFilename(prefix, spec), spec);
        }
    }
    catch (const std::exception &e)
    {
        Fail_(e.what());
    }

    if (rank == 0)
    {
        const fs::path output_path = fs::path(paraview_path_) / "paradata.vtm";
        std::ofstream out(output_path);
        if (!out)
            Fail_("[IOModule][VTK] cannot open: " + output_path.string());

        out << "<?xml version=\"1.0\"?>\n"
            << "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\">\n"
            << "  <vtkMultiBlockDataSet>\n";

        const int nblock = grd_->nblock;
        for (int ib = 0; ib < nblock; ++ib)
        {
            const std::string prefix = RankBlockPrefix(rank, ib);
            out << "    <Block index=\"" << ib << "\" name=\"" << prefix << "\">\n";
            for (std::size_t idata = 0; idata < kParaViewDatasets.size(); ++idata)
            {
                const ParaViewDatasetSpec &spec = kParaViewDatasets[idata];
                out << "      <DataSet index=\"" << idata << "\" name=\"" << spec.name
                    << "\" file=\"" << DatasetFilename(prefix, spec) << "\"/>\n";
            }
            out << "    </Block>\n";
        }

        out << "  </vtkMultiBlockDataSet>\n"
            << "</VTKFile>\n";

        std::printf("[IOModule][VTK] ParaView output is called: %s\n",
                    output_path.string().c_str());
        std::fflush(stdout);
    }

    PARALLEL::mpi_barrier();
}
