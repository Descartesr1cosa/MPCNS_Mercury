#pragma once
#include "0_basic/TYPES.h"
#include <array>
#include <string>
#include <vector>
namespace mpcns::mesh {
struct BlockAddress { int rank=-1; int local_block=-1; };
struct IndexTransform { std::array<int,3> permutation{}; std::array<int,3> sign{}; Int3 offset{}; };
enum class PatchKind { Inner, Parallel, Physical };
struct LocalBlock { BlockAddress address; std::string name; Box3 node_extent; };
struct Patch { PatchKind kind{}; BlockAddress source; BlockAddress neighbor; std::string source_name; std::string neighbor_name; Box3 source_node_box{}; Box3 neighbor_node_box{}; IndexTransform transform{}; bool has_neighbor=false; bool coupling=false; int direction=0; int neighbor_direction=0; int directions[3]={0,0,0}; int boundary_id=0; std::string boundary_name; };
class MeshPartition { public: MeshPartition(int,int,std::vector<LocalBlock>,std::vector<Patch>,std::vector<Patch>,std::vector<Patch>); int rank()const noexcept{return rank_;} int size()const noexcept{return size_;} const std::vector<LocalBlock>& local_blocks()const noexcept{return blocks_;} const std::vector<Patch>& face_patches()const noexcept{return faces_;} const std::vector<Patch>& edge_patches()const noexcept{return edges_;} const std::vector<Patch>& vertex_patches()const noexcept{return vertices_;} private:int rank_,size_;std::vector<LocalBlock> blocks_;std::vector<Patch> faces_,edges_,vertices_;};
}