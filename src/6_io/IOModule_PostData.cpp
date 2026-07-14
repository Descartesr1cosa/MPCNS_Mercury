#include "6_io/IOModule.h"

#include <stdexcept>

void IOModule::WritePostStaticData(const std::string &output_directory,
                                   POST::WriteOptions options) const
{
    if (!post_topology_)
        throw std::logic_error("IOModule::WritePostStaticData: call SetPostDataContext first");
    POST::PostDataWriter writer(*grd_, *post_topology_, *fld_, run_, myid_,
                                post_singular_edges_);
    writer.WriteStaticData(output_directory, std::move(options));
}
