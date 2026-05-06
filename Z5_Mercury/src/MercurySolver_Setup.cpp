#include "MercurySolver.h"

#include "0_MercuryFieldCatalog.h"
#include "2_topology/2_MPCNS_Topology.h"
#include "3_field/2_MPCNS_Field.h"
#include "4_halo/1_MPCNS_Halo.h"

void MercurySolver::RegisterFields(Field *fld, int ngg)
{
    MERCURY_FIELD::RegisterFields(fld, ngg);
}

void MercurySolver::RegisterCouplingChannels(Field *fld, const TOPO::Topology &topology, int dimension, int ngg)
{
    (void)ngg;
    MERCURY_FIELD::RegisterCouplingChannels(fld, topology, dimension);
}

void MercurySolver::RegisterHaloFields(Halo *halo)
{
    MERCURY_FIELD::RegisterHaloFields(halo);
}
