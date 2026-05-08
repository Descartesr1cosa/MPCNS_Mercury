#include "MercurySolver.h"

#include "0_MercuryFieldCatalog.h"
#include "2_topology/TopologyBuilder.h"
#include "3_field/Field.h"
#include "4_halo/Halo.h"

void MercurySolver::RegisterFields(Field *fld, int ngg)
{
    MERCURY_FIELD::RegisterFields(fld, ngg);
}

void MercurySolver::RegisterCouplingChannels(Field *fld, const TOPO::Topology &topology, int dimension, int ngg)
{
    (void)ngg;
    MERCURY_FIELD::RegisterCouplingChannels(fld, topology, dimension);
}

void MercurySolver::RegisterHaloFields(Field *fld, Halo *halo)
{
    MERCURY_FIELD::RegisterHaloFields(fld, halo);
}
