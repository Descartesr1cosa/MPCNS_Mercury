#pragma once
#include "mpcns/mesh/mesh_partition.hpp"
class Grid; namespace TOPO { struct Topology; }
namespace mpcns::compat { class LegacyMeshPartitionAdapter { public: static mesh::MeshPartition build(const Grid&,const TOPO::Topology&,int,int); }; }