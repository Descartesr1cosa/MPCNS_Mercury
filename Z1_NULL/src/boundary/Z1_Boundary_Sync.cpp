#include "Z1_Boundary.h"

#include "4_halo/Halo.h"

void Z1_Boundary::SyncGroup(const std::string &group)
{
    if (halo_)
        halo_->sync_group(group);
}

void Z1_Boundary::SyncAllRegistered()
{
    if (halo_)
        halo_->sync_registered();
}
