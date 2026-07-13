#pragma once

#include "0_basic/StaggerLocation.h"
#include "3_field/FieldDescriptor.h"
#include "3_field/FieldValueKind.h"

#include <string>
#include <vector>

namespace HALO_SYNC
{
    // Registered fields are classified into one of these execution families:
    // plain component-copy fields, orientation-aware edge 1-form triplets,
    // or orientation-aware face 2-form triplets.
    enum class Semantics
    {
        ComponentCopy,
        Edge1FormTriplet,
        Face2FormTriplet
    };

    struct TripletRequest
    {
        std::string group_name;

        std::string xi;
        std::string eta;
        std::string zeta;

        FieldValueKind value_kind = FieldValueKind::Scalar;

        HaloLevel level = HaloLevel::Corner3D;
        int nghost = 0;

        bool orientation_aware = true;
    };

    struct OwnerRequest
    {
        std::string field_name;
        std::string sync_group;

        OwnerSyncPolicy policy = OwnerSyncPolicy::None;

        FieldValueKind value_kind = FieldValueKind::Scalar;
        StaggerLocation location = StaggerLocation::Cell;
        int ncomp = 1;

        bool orientation_aware = false;
    };

    struct OwnerLocalOp
    {
        int owner_fid = -1;
        int alias_fid = -1;

        int owner_block = -1;
        int owner_i = 0;
        int owner_j = 0;
        int owner_k = 0;

        int alias_block = -1;
        int alias_i = 0;
        int alias_j = 0;
        int alias_k = 0;

        int ncomp = 1;
        int sign = +1;
    };

    struct OwnerSendOp
    {
        int owner_fid = -1;
        int class_gid = -1;

        int owner_block = -1;
        int owner_i = 0;
        int owner_j = 0;
        int owner_k = 0;

        int alias_rank = -1;
        int alias_block = -1;
        int alias_i = 0;
        int alias_j = 0;
        int alias_k = 0;

        int ncomp = 1;
        int sign_for_alias = +1;

        int dst_rank = -1;
        int tag = -1;
        int buffer_offset = 0;
    };

    struct OwnerRecvOp
    {
        int alias_fid = -1;
        int class_gid = -1;

        int alias_block = -1;
        int alias_i = 0;
        int alias_j = 0;
        int alias_k = 0;

        int ncomp = 1;

        int src_rank = -1;
        int tag = -1;
        int buffer_offset = 0;
    };

    struct OwnerPattern
    {
        std::string field_name;
        std::string sync_group;

        OwnerSyncPolicy policy = OwnerSyncPolicy::None;

        FieldValueKind value_kind = FieldValueKind::Scalar;
        StaggerLocation location = StaggerLocation::Cell;

        bool orientation_aware = false;

        std::vector<OwnerLocalOp> local_ops;
        std::vector<OwnerSendOp> send_ops;
        std::vector<OwnerRecvOp> recv_ops;

        std::vector<double> send_buffer;
        std::vector<double> recv_buffer;

        // Pattern sizes are immutable after build.  The cross-rank length
        // handshake is a construction-time invariant and only needs to pass
        // once, not on every field synchronization.
        bool mpi_lengths_validated = false;
    };

    struct CouplingPatternKey
    {
        std::string src;
        std::string dst;
        StaggerLocation loc = StaggerLocation::Cell;
        int nghost = 0;

        bool operator<(const CouplingPatternKey &o) const
        {
            if (src != o.src)
                return src < o.src;
            if (dst != o.dst)
                return dst < o.dst;
            if (static_cast<int>(loc) != static_cast<int>(o.loc))
                return static_cast<int>(loc) < static_cast<int>(o.loc);
            return nghost < o.nghost;
        }
    };
}
