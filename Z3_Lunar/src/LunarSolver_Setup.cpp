#include "LunarSolver.h"

#include "0_LunarFieldCatalog.h"
#include "2_topology/TopologyBuilder.h"
#include "3_field/Field.h"
#include "4_halo/Halo.h"

void LunarSolver::RegisterFields(Field *fld, int ngg)
{
    LUNAR_FIELD::RegisterFields(fld, ngg);
}

void LunarSolver::RegisterCouplingChannels(Field *fld, const TOPO::Topology &topology, int dimension, int ngg)
{
    (void)ngg;
    LUNAR_FIELD::RegisterCouplingChannels(fld, topology, dimension);
}

void LunarSolver::RegisterHaloFields(Field *fld, Halo *halo)
{
    LUNAR_FIELD::RegisterHaloFields(fld, halo);
}
