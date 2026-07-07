#include "Z0_Boundary.h"

#include "4_halo/Halo.h"

void Z0_Boundary::SyncGroup(const std::string &group)
{
    if (halo_)
        halo_->sync_group(group);
}

void Z0_Boundary::SyncGroup(const std::string &group, HaloLevel stage)
{
    if (halo_)
        halo_->sync_group(group, stage);
}

void Z0_Boundary::SyncAllRegistered()
{
    if (halo_)
        halo_->sync_registered();
}

void Z0_Boundary::SyncAllRegistered(HaloLevel stage)
{
    if (halo_)
        halo_->sync_registered(stage);
}
