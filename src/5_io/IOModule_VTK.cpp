#include "5_io/IOModule.h"

#include <cstdio>
#include <filesystem>
#include <fstream>

namespace fs = std::filesystem;

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
        std::ofstream out(output_path);
        if (!out)
            Fail_("[IOModule][VTK] cannot open: " + output_path.string());

        out << "<?xml version=\"1.0\"?>\n"
            << "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\">\n"
            << "  <vtkMultiBlockDataSet>\n"
            << "  </vtkMultiBlockDataSet>\n"
            << "</VTKFile>\n";

        std::printf("[IOModule][VTK] ParaView output is called: %s\n",
                    output_path.string().c_str());
        std::fflush(stdout);
    }

    if (par_)
        PARALLEL::mpi_barrier();
}
