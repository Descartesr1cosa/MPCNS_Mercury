#include "3_field/Field.h"

#include "7_metric/Metric.h"

void Field::build_geometry(int geometry_ghost)
{
    METRIC::build_field_geometry(*this, *grd, geometry_ghost);
}
