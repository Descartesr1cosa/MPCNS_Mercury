#include "2_topology/TopologyBuilder.h"
#include "2_topology/TopologyOps.h"

void TOPO::node_box_from_subsup(const int sub[3], const int sup[3], Box3 &box)
{
    TOPO::fill_node_box_from_subsup(sub, sup, box);
}
