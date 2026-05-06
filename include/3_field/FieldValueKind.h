#pragma once

enum class FieldValueKind
{
    Scalar,
    CartesianVector,
    ConservedVector,
    EdgeCovariant1Form,
    FaceContravariant2Form,
    GeometryMetric,
    CouplingBufferOnly,
    Temporary
};

inline const char *field_value_kind_name(FieldValueKind kind)
{
    switch (kind)
    {
    case FieldValueKind::Scalar:
        return "Scalar";
    case FieldValueKind::CartesianVector:
        return "CartesianVector";
    case FieldValueKind::ConservedVector:
        return "ConservedVector";
    case FieldValueKind::EdgeCovariant1Form:
        return "EdgeCovariant1Form";
    case FieldValueKind::FaceContravariant2Form:
        return "FaceContravariant2Form";
    case FieldValueKind::GeometryMetric:
        return "GeometryMetric";
    case FieldValueKind::CouplingBufferOnly:
        return "CouplingBufferOnly";
    case FieldValueKind::Temporary:
        return "Temporary";
    }

    return "Unknown";
}
