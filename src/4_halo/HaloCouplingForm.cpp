#include "4_halo/1_MPCNS_Halo.h"
#include "0_basic/Error.h"
#include "4_halo/detail/halo_build_tools.h"

namespace
{
std::string descriptor_group(const FieldDescriptor &desc)
{
    return desc.sync.group.empty() ? desc.name : desc.sync.group;
}
} // namespace

bool Halo::coupling_channel_needs_form_transfer_(const CouplingChannelSpec &ch) const
{
    return ch.orientation_aware &&
           (ch.value_kind == FieldValueKind::EdgeCovariant1Form ||
            ch.value_kind == FieldValueKind::FaceContravariant2Form);
}

int Halo::coupling_form_axis_from_location_(StaggerLocation loc) const
{
    switch (loc)
    {
    case StaggerLocation::EdgeXi:
    case StaggerLocation::FaceXi:
        return 0;
    case StaggerLocation::EdgeEt:
    case StaggerLocation::FaceEt:
        return 1;
    case StaggerLocation::EdgeZe:
    case StaggerLocation::FaceZe:
        return 2;
    default:
        ERROR::Abort("[Halo] coupling form transfer requires edge/face location");
    }
}

StaggerLocation Halo::coupling_form_location_from_axis_(FieldValueKind value_kind, int axis) const
{
    if (axis < 0 || axis > 2)
        ERROR::Abort("[Halo] coupling form transfer got invalid axis");

    if (value_kind == FieldValueKind::EdgeCovariant1Form)
    {
        if (axis == 0)
            return StaggerLocation::EdgeXi;
        if (axis == 1)
            return StaggerLocation::EdgeEt;
        return StaggerLocation::EdgeZe;
    }

    if (value_kind == FieldValueKind::FaceContravariant2Form)
    {
        if (axis == 0)
            return StaggerLocation::FaceXi;
        if (axis == 1)
            return StaggerLocation::FaceEt;
        return StaggerLocation::FaceZe;
    }

    ERROR::Abort("[Halo] coupling form transfer requires edge 1-form or face 2-form");
}

int Halo::coupling_src_axis_from_dst_to_src_transform_(const TOPO::IndexTransform &tr,
                                                       int dst_axis) const
{
    if (dst_axis < 0 || dst_axis > 2)
        ERROR::Abort("[Halo] coupling form transfer got invalid destination axis");
    return tr.perm[dst_axis];
}

int Halo::coupling_src_axis_from_src_to_dst_transform_(const TOPO::IndexTransform &tr,
                                                       int dst_axis) const
{
    if (dst_axis < 0 || dst_axis > 2)
        ERROR::Abort("[Halo] coupling form transfer got invalid destination axis");

    for (int a = 0; a < 3; ++a)
        if (tr.perm[a] == dst_axis)
            return a;

    ERROR::Abort("[Halo] coupling form transfer got invalid IndexTransform permutation");
}

int Halo::coupling_edge_1form_sign_dst_to_src_(const TOPO::IndexTransform &tr,
                                               int dst_axis) const
{
    if (dst_axis < 0 || dst_axis > 2)
        ERROR::Abort("[Halo] coupling edge 1-form sign got invalid destination axis");
    return tr.sign[dst_axis] >= 0 ? +1 : -1;
}

int Halo::coupling_edge_1form_sign_src_to_dst_(const TOPO::IndexTransform &tr,
                                               int src_axis) const
{
    if (src_axis < 0 || src_axis > 2)
        ERROR::Abort("[Halo] coupling edge 1-form sign got invalid source axis");
    return tr.sign[src_axis] >= 0 ? +1 : -1;
}

int Halo::coupling_face_2form_sign_dst_to_src_(const TOPO::IndexTransform &tr,
                                               int dst_axis) const
{
    const TOPO::IndexTransform src_to_dst = HALO_TOOLS::inverse_transform(tr);
    const int src_axis = coupling_src_axis_from_dst_to_src_transform_(tr, dst_axis);
    return HALO_TOOLS::face_2form_orientation_sign(src_to_dst, src_axis);
}

int Halo::coupling_face_2form_sign_src_to_dst_(const TOPO::IndexTransform &tr,
                                               int src_axis) const
{
    return HALO_TOOLS::face_2form_orientation_sign(tr, src_axis);
}

int Halo::coupling_form_orientation_sign_dst_to_src_(const CouplingChannelSpec &ch,
                                                     const TOPO::IndexTransform &tr,
                                                     int dst_axis) const
{
    if (ch.value_kind == FieldValueKind::EdgeCovariant1Form)
        return coupling_edge_1form_sign_dst_to_src_(tr, dst_axis);
    if (ch.value_kind == FieldValueKind::FaceContravariant2Form)
        return coupling_face_2form_sign_dst_to_src_(tr, dst_axis);

    return +1;
}

int Halo::coupling_form_orientation_sign_src_to_dst_(const CouplingChannelSpec &ch,
                                                     const TOPO::IndexTransform &tr,
                                                     int src_axis) const
{
    if (ch.value_kind == FieldValueKind::EdgeCovariant1Form)
        return coupling_edge_1form_sign_src_to_dst_(tr, src_axis);
    if (ch.value_kind == FieldValueKind::FaceContravariant2Form)
        return coupling_face_2form_sign_src_to_dst_(tr, src_axis);

    return +1;
}

std::string Halo::find_triplet_field_name_(const std::string &dst_field_name,
                                           FieldValueKind value_kind,
                                           int wanted_src_axis) const
{
    if (!fld_->has_field(dst_field_name))
        ERROR::Abort("[Halo] coupling form transfer cannot find field: " + dst_field_name);

    const FieldDescriptor &dst_desc = fld_->descriptor(dst_field_name);
    const std::string group = descriptor_group(dst_desc);
    if (group.empty())
        ERROR::Abort("[Halo] coupling form transfer field has no sync group: " + dst_field_name);

    const StaggerLocation wanted_loc =
        coupling_form_location_from_axis_(value_kind, wanted_src_axis);

    for (const FieldDescriptor &desc : fld_->descriptors())
    {
        if (descriptor_group(desc) != group)
            continue;
        if (desc.value_kind != value_kind)
            continue;
        if (desc.location != wanted_loc)
            continue;
        return desc.name;
    }

    ERROR::Abort("[Halo] coupling form transfer cannot find triplet field in group " +
                 group + " for field " + dst_field_name);
}
