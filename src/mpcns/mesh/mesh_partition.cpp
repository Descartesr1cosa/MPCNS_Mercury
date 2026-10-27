#include "mpcns/mesh/mesh_partition.hpp"
#include <stdexcept>
namespace mpcns::mesh { MeshPartition::MeshPartition(int r,int s,std::vector<LocalBlock>b,std::vector<Patch>f,std::vector<Patch>e,std::vector<Patch>v):rank_(r),size_(s),blocks_(std::move(b)),faces_(std::move(f)),edges_(std::move(e)),vertices_(std::move(v)){if(r<0||s<1||r>=s)throw std::invalid_argument("invalid MPI rank snapshot");} }