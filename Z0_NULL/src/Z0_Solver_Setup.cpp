#include "Z0_Solver.h"

#include "Z0_Control.h"

#include "3_field/Field.h"

#include <set>
#include <string>
#include <vector>

void Z0_Solver::Setup_()
{
    if (setup_done_)
        return;
    control_ = Z0::LoadControl(*param_);

    param_->AddParam("Archive_Output_Time", 0.0);
    io_.Setup(param_, grid_, field_, 1);
    io_.SetTecplotOutputMode(IOModule::TecplotMode::AllNode);
    io_.SetTecplotFormReconstruction(IOModule::TecplotFormReconstruction::ToNode);

    std::vector<std::string> output_fields;
    std::set<std::string> seen;
    for (const FieldDescriptor &desc : field_->descriptors())
    {
        if (desc.sync.group.empty())
            continue;
        if (seen.insert(desc.name).second)
            output_fields.push_back(desc.name);
    }
    io_.SetTecplotPhysicalOutputs(output_fields);

    setup_done_ = true;
}
