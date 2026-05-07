#include "5_io/IOModule.h"
#include "5_io/VTKXmlAppendedWriter.h"

#include <array>
#include <algorithm>
#include <cmath>
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

    int LocationCode(StaggerLocation loc)
    {
        switch (loc)
        {
        case StaggerLocation::EdgeXi:
            return 1;
        case StaggerLocation::EdgeEt:
            return 2;
        case StaggerLocation::EdgeZe:
            return 3;
        case StaggerLocation::FaceXi:
            return 4;
        case StaggerLocation::FaceEt:
            return 5;
        case StaggerLocation::FaceZe:
            return 6;
        default:
            return 0;
        }
    }

    void EdgeDofSize(const Block &block, StaggerLocation loc, int &ni, int &nj, int &nk)
    {
        const int mx = block.mx;
        const int my = block.my;
        const int mz = (block.dimension == 2) ? 0 : block.mz;

        if (loc == StaggerLocation::EdgeXi)
        {
            ni = mx;
            nj = my + 1;
            nk = mz + 1;
        }
        else if (loc == StaggerLocation::EdgeEt)
        {
            ni = mx + 1;
            nj = my;
            nk = mz + 1;
        }
        else if (loc == StaggerLocation::EdgeZe)
        {
            ni = mx + 1;
            nj = my + 1;
            nk = mz;
        }
        else
        {
            ni = 0;
            nj = 0;
            nk = 0;
        }
    }

    void EdgeEndpoint1(StaggerLocation loc, int i, int j, int k, int &i1, int &j1, int &k1)
    {
        i1 = i;
        j1 = j;
        k1 = k;

        if (loc == StaggerLocation::EdgeXi)
            ++i1;
        else if (loc == StaggerLocation::EdgeEt)
            ++j1;
        else if (loc == StaggerLocation::EdgeZe)
            ++k1;
    }

    void AddEdgeDebugArrays(VTKXML::AppendedRawWriter &appended,
                            std::vector<VTKXML::AppendedArray> &cell_arrays,
                            const std::vector<std::int32_t> &i_values,
                            const std::vector<std::int32_t> &j_values,
                            const std::vector<std::int32_t> &k_values,
                            const std::vector<std::int32_t> &block_values,
                            const std::vector<std::int32_t> &rank_values,
                            const std::vector<std::int32_t> &location_values,
                            const std::vector<double> &edge_dx,
                            const std::vector<double> &edge_dy,
                            const std::vector<double> &edge_dz,
                            const std::vector<double> &edge_length)
    {
        cell_arrays.push_back(appended.AddInt32("i", 1, i_values));
        cell_arrays.push_back(appended.AddInt32("j", 1, j_values));
        cell_arrays.push_back(appended.AddInt32("k", 1, k_values));
        cell_arrays.push_back(appended.AddInt32("block_id", 1, block_values));
        cell_arrays.push_back(appended.AddInt32("rank", 1, rank_values));
        cell_arrays.push_back(appended.AddInt32("location_code", 1, location_values));
        cell_arrays.push_back(appended.AddFloat64("edge_dx", 1, edge_dx));
        cell_arrays.push_back(appended.AddFloat64("edge_dy", 1, edge_dy));
        cell_arrays.push_back(appended.AddFloat64("edge_dz", 1, edge_dz));
        cell_arrays.push_back(appended.AddFloat64("edge_length", 1, edge_length));
    }

    void FaceDofSize(const Block &block, StaggerLocation loc, int &ni, int &nj, int &nk)
    {
        const int mx = block.mx;
        const int my = block.my;
        const int mz = (block.dimension == 2) ? 0 : block.mz;

        if (loc == StaggerLocation::FaceXi)
        {
            ni = mx + 1;
            nj = my;
            nk = mz;
        }
        else if (loc == StaggerLocation::FaceEt)
        {
            ni = mx;
            nj = my + 1;
            nk = mz;
        }
        else if (loc == StaggerLocation::FaceZe)
        {
            ni = mx;
            nj = my;
            nk = mz + 1;
        }
        else
        {
            ni = 0;
            nj = 0;
            nk = 0;
        }
    }

    void FaceCorner(StaggerLocation loc, int i, int j, int k, int corner, int &ic, int &jc, int &kc)
    {
        ic = i;
        jc = j;
        kc = k;

        if (loc == StaggerLocation::FaceXi)
        {
            if (corner == 1 || corner == 2)
                ++jc;
            if (corner == 2 || corner == 3)
                ++kc;
        }
        else if (loc == StaggerLocation::FaceEt)
        {
            if (corner == 2 || corner == 3)
                ++ic;
            if (corner == 1 || corner == 2)
                ++kc;
        }
        else if (loc == StaggerLocation::FaceZe)
        {
            if (corner == 1 || corner == 2)
                ++ic;
            if (corner == 2 || corner == 3)
                ++jc;
        }
    }

    void AddFaceDebugArrays(VTKXML::AppendedRawWriter &appended,
                            std::vector<VTKXML::AppendedArray> &cell_arrays,
                            const std::vector<std::int32_t> &i_values,
                            const std::vector<std::int32_t> &j_values,
                            const std::vector<std::int32_t> &k_values,
                            const std::vector<std::int32_t> &block_values,
                            const std::vector<std::int32_t> &rank_values,
                            const std::vector<std::int32_t> &location_values,
                            const std::vector<double> &normal_x,
                            const std::vector<double> &normal_y,
                            const std::vector<double> &normal_z,
                            const std::vector<double> &face_area)
    {
        cell_arrays.push_back(appended.AddInt32("i", 1, i_values));
        cell_arrays.push_back(appended.AddInt32("j", 1, j_values));
        cell_arrays.push_back(appended.AddInt32("k", 1, k_values));
        cell_arrays.push_back(appended.AddInt32("block_id", 1, block_values));
        cell_arrays.push_back(appended.AddInt32("rank", 1, rank_values));
        cell_arrays.push_back(appended.AddInt32("location_code", 1, location_values));
        cell_arrays.push_back(appended.AddFloat64("normal_x", 1, normal_x));
        cell_arrays.push_back(appended.AddFloat64("normal_y", 1, normal_y));
        cell_arrays.push_back(appended.AddFloat64("normal_z", 1, normal_z));
        cell_arrays.push_back(appended.AddFloat64("face_area", 1, face_area));
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

    void WriteEdgeVTP(const fs::path &path,
                      Block &block,
                      Field &field,
                      int iblock,
                      int rank,
                      StaggerLocation edge_location,
                      const std::vector<std::string> &white_list)
    {
        VTKXML::AppendedRawWriter appended;

        int ni = 0;
        int nj = 0;
        int nk = 0;
        EdgeDofSize(block, edge_location, ni, nj, nk);

        const std::size_t nedge = static_cast<std::size_t>(ni) * nj * nk;

        std::vector<double> point_values;
        std::vector<std::int64_t> connectivity;
        std::vector<std::int64_t> offsets;
        std::vector<std::int32_t> i_values;
        std::vector<std::int32_t> j_values;
        std::vector<std::int32_t> k_values;
        std::vector<std::int32_t> block_values;
        std::vector<std::int32_t> rank_values;
        std::vector<std::int32_t> location_values;
        std::vector<double> edge_dx;
        std::vector<double> edge_dy;
        std::vector<double> edge_dz;
        std::vector<double> edge_length;

        point_values.reserve(nedge * 2 * 3);
        connectivity.reserve(nedge * 2);
        offsets.reserve(nedge);
        i_values.reserve(nedge);
        j_values.reserve(nedge);
        k_values.reserve(nedge);
        block_values.reserve(nedge);
        rank_values.reserve(nedge);
        location_values.reserve(nedge);
        edge_dx.reserve(nedge);
        edge_dy.reserve(nedge);
        edge_dz.reserve(nedge);
        edge_length.reserve(nedge);

        std::int64_t point_index = 0;
        for (int k = 0; k < nk; ++k)
        {
            for (int j = 0; j < nj; ++j)
            {
                for (int i = 0; i < ni; ++i)
                {
                    int i1 = i;
                    int j1 = j;
                    int k1 = k;
                    EdgeEndpoint1(edge_location, i, j, k, i1, j1, k1);

                    const double x0 = block.x(i, j, k);
                    const double y0 = block.y(i, j, k);
                    const double z0 = block.z(i, j, k);
                    const double x1 = block.x(i1, j1, k1);
                    const double y1 = block.y(i1, j1, k1);
                    const double z1 = block.z(i1, j1, k1);
                    const double dx = x1 - x0;
                    const double dy = y1 - y0;
                    const double dz = z1 - z0;

                    point_values.push_back(x0);
                    point_values.push_back(y0);
                    point_values.push_back(z0);
                    point_values.push_back(x1);
                    point_values.push_back(y1);
                    point_values.push_back(z1);

                    connectivity.push_back(point_index);
                    connectivity.push_back(point_index + 1);
                    point_index += 2;
                    offsets.push_back(point_index);

                    i_values.push_back(i);
                    j_values.push_back(j);
                    k_values.push_back(k);
                    block_values.push_back(iblock);
                    rank_values.push_back(rank);
                    location_values.push_back(LocationCode(edge_location));
                    edge_dx.push_back(dx);
                    edge_dy.push_back(dy);
                    edge_dz.push_back(dz);
                    edge_length.push_back(std::sqrt(dx * dx + dy * dy + dz * dz));
                }
            }
        }

        const VTKXML::AppendedArray points =
            appended.AddFloat64("", 3, point_values);
        const VTKXML::AppendedArray line_connectivity =
            appended.AddInt64("connectivity", 1, connectivity);
        const VTKXML::AppendedArray line_offsets =
            appended.AddInt64("offsets", 1, offsets);

        std::vector<VTKXML::AppendedArray> cell_arrays;
        AddEdgeDebugArrays(appended, cell_arrays,
                           i_values, j_values, k_values,
                           block_values, rank_values, location_values,
                           edge_dx, edge_dy, edge_dz, edge_length);

        const int nfield = field.num_fields();
        for (int fid = 0; fid < nfield; ++fid)
        {
            const FieldDescriptor &desc = field.descriptor(fid);
            if (!FieldSelected(white_list, desc.name))
                continue;
            if (desc.location != edge_location)
                continue;

            FieldBlock &field_block = field.field(fid, iblock);
            if (!field_block.is_allocated())
                continue;

            AddFieldDataArrays(appended, cell_arrays, field_block, desc.name, desc.ncomp, ni, nj, nk);
        }

        std::ofstream out(path, std::ios::binary);
        if (!out)
            throw std::runtime_error("[IOModule][VTK] cannot open: " + path.string());

        out << "<?xml version=\"1.0\"?>\n"
            << "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n"
            << "  <PolyData>\n"
            << "    <Piece NumberOfPoints=\"" << nedge * 2
            << "\" NumberOfVerts=\"0\" NumberOfLines=\"" << nedge
            << "\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n"
            << "      <PointData>\n"
            << "      </PointData>\n"
            << "      <CellData>\n";
        for (const VTKXML::AppendedArray &array : cell_arrays)
            appended.WriteDataArrayTag(out, array, 8);
        out << "      </CellData>\n"
            << "      <Points>\n";
        appended.WriteDataArrayTag(out, points, 8);
        out << "      </Points>\n"
            << "      <Lines>\n";
        appended.WriteDataArrayTag(out, line_connectivity, 8);
        appended.WriteDataArrayTag(out, line_offsets, 8);
        out << "      </Lines>\n"
            << "      <Verts>\n"
            << "      </Verts>\n"
            << "      <Strips>\n"
            << "      </Strips>\n"
            << "      <Polys>\n"
            << "      </Polys>\n"
            << "    </Piece>\n"
            << "  </PolyData>\n";
        appended.WriteAppendedData(out);
        out << "</VTKFile>\n";
    }

    void WriteFaceVTP(const fs::path &path,
                      Block &block,
                      Field &field,
                      int iblock,
                      int rank,
                      StaggerLocation face_location,
                      const std::vector<std::string> &white_list)
    {
        VTKXML::AppendedRawWriter appended;

        int ni = 0;
        int nj = 0;
        int nk = 0;
        FaceDofSize(block, face_location, ni, nj, nk);

        const std::size_t nface = static_cast<std::size_t>(ni) * nj * nk;

        std::vector<double> point_values;
        std::vector<std::int64_t> connectivity;
        std::vector<std::int64_t> offsets;
        std::vector<std::int32_t> i_values;
        std::vector<std::int32_t> j_values;
        std::vector<std::int32_t> k_values;
        std::vector<std::int32_t> block_values;
        std::vector<std::int32_t> rank_values;
        std::vector<std::int32_t> location_values;
        std::vector<double> normal_x;
        std::vector<double> normal_y;
        std::vector<double> normal_z;
        std::vector<double> face_area;

        point_values.reserve(nface * 4 * 3);
        connectivity.reserve(nface * 4);
        offsets.reserve(nface);
        i_values.reserve(nface);
        j_values.reserve(nface);
        k_values.reserve(nface);
        block_values.reserve(nface);
        rank_values.reserve(nface);
        location_values.reserve(nface);
        normal_x.reserve(nface);
        normal_y.reserve(nface);
        normal_z.reserve(nface);
        face_area.reserve(nface);

        std::int64_t point_index = 0;
        for (int k = 0; k < nk; ++k)
        {
            for (int j = 0; j < nj; ++j)
            {
                for (int i = 0; i < ni; ++i)
                {
                    double x[4] = {};
                    double y[4] = {};
                    double z[4] = {};
                    for (int c = 0; c < 4; ++c)
                    {
                        int ic = i;
                        int jc = j;
                        int kc = k;
                        FaceCorner(face_location, i, j, k, c, ic, jc, kc);
                        x[c] = block.x(ic, jc, kc);
                        y[c] = block.y(ic, jc, kc);
                        z[c] = block.z(ic, jc, kc);

                        point_values.push_back(x[c]);
                        point_values.push_back(y[c]);
                        point_values.push_back(z[c]);
                        connectivity.push_back(point_index++);
                    }
                    offsets.push_back(point_index);

                    const double ax = x[1] - x[0];
                    const double ay = y[1] - y[0];
                    const double az = z[1] - z[0];
                    const double bx = x[3] - x[0];
                    const double by = y[3] - y[0];
                    const double bz = z[3] - z[0];
                    double nx = ay * bz - az * by;
                    double ny = az * bx - ax * bz;
                    double nz = ax * by - ay * bx;
                    const double area1 = std::sqrt(nx * nx + ny * ny + nz * nz);

                    const double cx = x[2] - x[1];
                    const double cy = y[2] - y[1];
                    const double cz = z[2] - z[1];
                    const double dx = x[0] - x[1];
                    const double dy = y[0] - y[1];
                    const double dz = z[0] - z[1];
                    const double nx2 = cy * dz - cz * dy;
                    const double ny2 = cz * dx - cx * dz;
                    const double nz2 = cx * dy - cy * dx;
                    const double area2 = std::sqrt(nx2 * nx2 + ny2 * ny2 + nz2 * nz2);
                    const double area = 0.5 * (area1 + area2);

                    const double nmag = std::sqrt(nx * nx + ny * ny + nz * nz);
                    if (nmag > 0.0)
                    {
                        nx /= nmag;
                        ny /= nmag;
                        nz /= nmag;
                    }

                    i_values.push_back(i);
                    j_values.push_back(j);
                    k_values.push_back(k);
                    block_values.push_back(iblock);
                    rank_values.push_back(rank);
                    location_values.push_back(LocationCode(face_location));
                    normal_x.push_back(nx);
                    normal_y.push_back(ny);
                    normal_z.push_back(nz);
                    face_area.push_back(area);
                }
            }
        }

        const VTKXML::AppendedArray points =
            appended.AddFloat64("", 3, point_values);
        const VTKXML::AppendedArray poly_connectivity =
            appended.AddInt64("connectivity", 1, connectivity);
        const VTKXML::AppendedArray poly_offsets =
            appended.AddInt64("offsets", 1, offsets);

        std::vector<VTKXML::AppendedArray> cell_arrays;
        AddFaceDebugArrays(appended, cell_arrays,
                           i_values, j_values, k_values,
                           block_values, rank_values, location_values,
                           normal_x, normal_y, normal_z, face_area);

        const int nfield = field.num_fields();
        for (int fid = 0; fid < nfield; ++fid)
        {
            const FieldDescriptor &desc = field.descriptor(fid);
            if (!FieldSelected(white_list, desc.name))
                continue;
            if (desc.location != face_location)
                continue;

            FieldBlock &field_block = field.field(fid, iblock);
            if (!field_block.is_allocated())
                continue;

            AddFieldDataArrays(appended, cell_arrays, field_block, desc.name, desc.ncomp, ni, nj, nk);
        }

        std::ofstream out(path, std::ios::binary);
        if (!out)
            throw std::runtime_error("[IOModule][VTK] cannot open: " + path.string());

        out << "<?xml version=\"1.0\"?>\n"
            << "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n"
            << "  <PolyData>\n"
            << "    <Piece NumberOfPoints=\"" << nface * 4
            << "\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"" << nface << "\">\n"
            << "      <PointData>\n"
            << "      </PointData>\n"
            << "      <CellData>\n";
        for (const VTKXML::AppendedArray &array : cell_arrays)
            appended.WriteDataArrayTag(out, array, 8);
        out << "      </CellData>\n"
            << "      <Points>\n";
        appended.WriteDataArrayTag(out, points, 8);
        out << "      </Points>\n"
            << "      <Polys>\n";
        appended.WriteDataArrayTag(out, poly_connectivity, 8);
        appended.WriteDataArrayTag(out, poly_offsets, 8);
        out << "      </Polys>\n"
            << "      <Verts>\n"
            << "      </Verts>\n"
            << "      <Lines>\n"
            << "      </Lines>\n"
            << "      <Strips>\n"
            << "      </Strips>\n"
            << "    </Piece>\n"
            << "  </PolyData>\n";
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

    bool IsEdgeDataset(const ParaViewDatasetSpec &spec, StaggerLocation &loc)
    {
        const std::string name = spec.name;
        if (name == "edge_xi")
        {
            loc = StaggerLocation::EdgeXi;
            return true;
        }
        if (name == "edge_eta")
        {
            loc = StaggerLocation::EdgeEt;
            return true;
        }
        if (name == "edge_zeta")
        {
            loc = StaggerLocation::EdgeZe;
            return true;
        }
        return false;
    }

    bool IsFaceDataset(const ParaViewDatasetSpec &spec, StaggerLocation &loc)
    {
        const std::string name = spec.name;
        if (name == "face_xi")
        {
            loc = StaggerLocation::FaceXi;
            return true;
        }
        if (name == "face_eta")
        {
            loc = StaggerLocation::FaceEt;
            return true;
        }
        if (name == "face_zeta")
        {
            loc = StaggerLocation::FaceZe;
            return true;
        }
        return false;
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
                StaggerLocation edge_location = StaggerLocation::EdgeXi;
                StaggerLocation face_location = StaggerLocation::FaceXi;
                if (std::string(spec.extension) == ".vts")
                    WriteVolumeVTS(path, grd_->grids(ib), *fld_, ib, paraview_fields_);
                else if (IsEdgeDataset(spec, edge_location))
                    WriteEdgeVTP(path, grd_->grids(ib), *fld_, ib, rank, edge_location, paraview_fields_);
                else if (IsFaceDataset(spec, face_location))
                    WriteFaceVTP(path, grd_->grids(ib), *fld_, ib, rank, face_location, paraview_fields_);
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
