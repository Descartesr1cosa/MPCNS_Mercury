#include "5_io/IOModule.h"
#include "5_io/VTKXmlAppendedWriter.h"

#include <array>
#include <algorithm>
#include <cstdio>
#include <cstdint>
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

    bool FieldSelected(const std::vector<std::string> &white_list, const std::string &name)
    {
        return white_list.empty() ||
               std::find(white_list.begin(), white_list.end(), name) != white_list.end();
    }

    void AddFieldDataArrays(VTKXML::AppendedRawWriter &appended,
                            std::vector<VTKXML::AppendedArray> &arrays,
                            FieldBlock &field_block,
                            const std::string &field_name,
                            int ncomp,
                            int ni,
                            int nj,
                            int nk)
    {
        if (ncomp == 1 || ncomp == 3)
        {
            std::vector<double> values;
            values.reserve(static_cast<std::size_t>(ni) * nj * nk * ncomp);
            for (int k = 0; k < nk; ++k)
                for (int j = 0; j < nj; ++j)
                    for (int i = 0; i < ni; ++i)
                        for (int m = 0; m < ncomp; ++m)
                            values.push_back(field_block(i, j, k, m));
            arrays.push_back(appended.AddFloat64(field_name, ncomp, values));
            return;
        }

        for (int m = 0; m < ncomp; ++m)
        {
            std::vector<double> values;
            values.reserve(static_cast<std::size_t>(ni) * nj * nk);
            for (int k = 0; k < nk; ++k)
                for (int j = 0; j < nj; ++j)
                    for (int i = 0; i < ni; ++i)
                        values.push_back(field_block(i, j, k, m));
            arrays.push_back(appended.AddFloat64(field_name + "_" + std::to_string(m), 1, values));
        }
    }

    void WriteVolumeVTS(const fs::path &path,
                        Block &block,
                        Field &field,
                        int iblock,
                        const std::vector<std::string> &white_list)
    {
        VTKXML::AppendedRawWriter appended;

        const int dim = block.dimension;
        const int mx = block.mx;
        const int my = block.my;
        const int mz = (dim == 2) ? 1 : block.mz;
        const int point_ni = mx + 1;
        const int point_nj = my + 1;
        const int point_nk = (dim == 2) ? 1 : (block.mz + 1);
        const int cell_ni = mx;
        const int cell_nj = my;
        const int cell_nk = mz;

        std::vector<double> point_values;
        point_values.reserve(static_cast<std::size_t>(point_ni) * point_nj * point_nk * 3);
        for (int k = 0; k < point_nk; ++k)
        {
            for (int j = 0; j < point_nj; ++j)
            {
                for (int i = 0; i < point_ni; ++i)
                {
                    point_values.push_back(block.x(i, j, k));
                    point_values.push_back(block.y(i, j, k));
                    point_values.push_back(block.z(i, j, k));
                }
            }
        }
        const VTKXML::AppendedArray points =
            appended.AddFloat64("", 3, point_values);

        std::vector<VTKXML::AppendedArray> point_arrays;
        std::vector<VTKXML::AppendedArray> cell_arrays;

        const int nfield = field.num_fields();
        for (int fid = 0; fid < nfield; ++fid)
        {
            const FieldDescriptor &desc = field.descriptor(fid);
            if (!FieldSelected(white_list, desc.name))
                continue;
            if (desc.location != StaggerLocation::Node && desc.location != StaggerLocation::Cell)
                continue;

            FieldBlock &field_block = field.field(fid, iblock);
            if (!field_block.is_allocated())
                continue;

            if (desc.location == StaggerLocation::Node)
            {
                AddFieldDataArrays(appended, point_arrays, field_block, desc.name, desc.ncomp,
                                   point_ni, point_nj, point_nk);
            }
            else
            {
                AddFieldDataArrays(appended, cell_arrays, field_block, desc.name, desc.ncomp,
                                   cell_ni, cell_nj, cell_nk);
            }
        }

        std::ofstream out(path, std::ios::binary);
        if (!out)
            throw std::runtime_error("[IOModule][VTK] cannot open: " + path.string());

        const int extent_z_hi = (dim == 2) ? 0 : block.mz;
        out << "<?xml version=\"1.0\"?>\n"
            << "<VTKFile type=\"StructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n"
            << "  <StructuredGrid WholeExtent=\"0 " << mx << " 0 " << my << " 0 " << extent_z_hi << "\">\n"
            << "    <Piece Extent=\"0 " << mx << " 0 " << my << " 0 " << extent_z_hi << "\">\n"
            << "      <PointData>\n";
        for (const VTKXML::AppendedArray &array : point_arrays)
            appended.WriteDataArrayTag(out, array, 8);
        out << "      </PointData>\n"
            << "      <CellData>\n";
        for (const VTKXML::AppendedArray &array : cell_arrays)
            appended.WriteDataArrayTag(out, array, 8);
        out << "      </CellData>\n"
            << "      <Points>\n";
        appended.WriteDataArrayTag(out, points, 8);
        out << "      </Points>\n"
            << "    </Piece>\n"
            << "  </StructuredGrid>\n";
        appended.WriteAppendedData(out);
        out << "</VTKFile>\n";
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
            {
                const fs::path path = fs::path(paraview_path_) / DatasetFilename(prefix, spec);
                if (std::string(spec.extension) == ".vts")
                    WriteVolumeVTS(path, grd_->grids(ib), *fld_, ib, paraview_fields_);
                else
                    WritePlaceholderDataset(path, spec);
            }
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
