#include "Z0_Tests.h"
#include "Z0_TestCommon.h"
#include "0_basic/MPI_WRAPPER.h"
#include "1_grid/1_MPCNS_Grid.h"
#include "2_topology/Topology.h"
#include "3_field/Field.h"
#include "4_halo/Halo.h"
#include "mpcns/compat/legacy_halo_plan_adapter.hpp"
#include "mpcns/compat/legacy_mesh_partition_adapter.hpp"
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {
using namespace mpcns;
using parallel::ExchangeKind;
struct Summary { long long checked=0, failed=0; };
int diagnostic_stage=-1;
int diagnostic_location=-1;
int diagnostic_nghost=-1;
const char* role(ExchangeKind k) {
 switch(k) {
 case ExchangeKind::InnerCopy:return "inner-copy";
 case ExchangeKind::ParallelFace:return "send+recv";
 case ExchangeKind::ParallelEdgeSend:
 case ExchangeKind::ParallelVertexSend:return "send";
 case ExchangeKind::ParallelEdgeRecv:
 case ExchangeKind::ParallelVertexRecv:return "recv";
 }
 return "unknown";
}
const char* name(ExchangeKind k) {
 switch(k) {
 case ExchangeKind::InnerCopy:return "InnerCopy";
 case ExchangeKind::ParallelFace:return "ParallelFace";
 case ExchangeKind::ParallelEdgeSend:return "ParallelEdgeSend";
 case ExchangeKind::ParallelEdgeRecv:return "ParallelEdgeRecv";
 case ExchangeKind::ParallelVertexSend:return "ParallelVertexSend";
 case ExchangeKind::ParallelVertexRecv:return "ParallelVertexRecv";
 } return "Unknown";
}
parallel::FieldLocation loc(StaggerLocation x) {
 switch(x) {
 case StaggerLocation::Cell:return parallel::FieldLocation::Cell;
 case StaggerLocation::Node:return parallel::FieldLocation::Node;
 case StaggerLocation::FaceXi:return parallel::FieldLocation::FaceXi;
 case StaggerLocation::FaceEt:return parallel::FieldLocation::FaceEta;
 case StaggerLocation::FaceZe:return parallel::FieldLocation::FaceZeta;
 case StaggerLocation::EdgeXi:return parallel::FieldLocation::EdgeXi;
 case StaggerLocation::EdgeEt:return parallel::FieldLocation::EdgeEta;
 case StaggerLocation::EdgeZe:return parallel::FieldLocation::EdgeZeta;
 } throw std::invalid_argument("location");
}
ExchangeKind kind(LegacyHaloPatternKind x) {
 switch(x) {
 case LegacyHaloPatternKind::InnerCopy:return ExchangeKind::InnerCopy;
 case LegacyHaloPatternKind::ParallelFace:return ExchangeKind::ParallelFace;
 case LegacyHaloPatternKind::ParallelEdgeSend:return ExchangeKind::ParallelEdgeSend;
 case LegacyHaloPatternKind::ParallelEdgeRecv:return ExchangeKind::ParallelEdgeRecv;
 case LegacyHaloPatternKind::ParallelVertexSend:return ExchangeKind::ParallelVertexSend;
 case LegacyHaloPatternKind::ParallelVertexRecv:return ExchangeKind::ParallelVertexRecv;
 } throw std::invalid_argument("kind");
}
parallel::HaloStage stage(HaloLevel x) {
 return x==HaloLevel::Corner1D?parallel::HaloStage::Face:
        x==HaloLevel::Corner2D?parallel::HaloStage::Edge:parallel::HaloStage::Vertex;
}
bool has_send(LegacyHaloPatternKind k) {
 return k==LegacyHaloPatternKind::InnerCopy||k==LegacyHaloPatternKind::ParallelFace||
        k==LegacyHaloPatternKind::ParallelEdgeSend||k==LegacyHaloPatternKind::ParallelVertexSend;
}
bool has_recv(LegacyHaloPatternKind k) {
 return k==LegacyHaloPatternKind::InnerCopy||k==LegacyHaloPatternKind::ParallelFace||
        k==LegacyHaloPatternKind::ParallelEdgeRecv||k==LegacyHaloPatternKind::ParallelVertexRecv;
}
template<class A,class B> void eq(Summary& s,int rank,const std::string& request,
 ExchangeKind k,size_t pi,size_t ri,const std::string& field,const A& e,const B& a) {
 ++s.checked; if(e==a)return; ++s.failed;
 std::cerr<<"[ExchangePlan][rank "<<rank<<"][request "<<request<<"][stage "
          <<diagnostic_stage<<"][location "<<diagnostic_location<<"][nghost "
          <<diagnostic_nghost<<"]["<<name(k)<<"][role "<<role(k)
          <<"][pattern "<<pi<<"][region "<<ri<<"] field="<<field
          <<" expected="<<e<<" actual="<<a<<"\n";
}
void box(Summary& s,int rank,const std::string& req,ExchangeKind k,size_t pi,size_t ri,
 const char* p,const Box3& e,const core::LogicalBox& a) {
 const int ev[6]={e.lo.i,e.lo.j,e.lo.k,e.hi.i,e.hi.j,e.hi.k};
 const int av[6]={a.lo.i,a.lo.j,a.lo.k,a.hi.i,a.hi.j,a.hi.k};
 const char* n[6]={".lo.i",".lo.j",".lo.k",".hi.i",".hi.j",".hi.k"};
 for(int i=0;i<6;++i)eq(s,rank,req,k,pi,ri,std::string(p)+n[i],ev[i],av[i]);
}
void region(Summary& s,int rank,const FieldHaloRequest& q,size_t pi,size_t ri,
 const LegacyHaloPatternSnapshot& p,const HaloRegion& e,const parallel::ExchangeRegion& a) {
 const ExchangeKind k=kind(p.kind);
 diagnostic_stage=(int)stage(p.stage);
 diagnostic_location=(int)loc(p.location);
 diagnostic_nghost=p.nghost;
 eq(s,rank,q.field_name,k,pi,ri,"kind",(int)k,(int)a.kind);
 eq(s,rank,q.field_name,k,pi,ri,"stage",(int)stage(p.stage),(int)a.stage);
 eq(s,rank,q.field_name,k,pi,ri,"location",(int)loc(p.location),(int)a.location);
 eq(s,rank,q.field_name,k,pi,ri,"ghost_layers",p.nghost,a.ghost_layers);
 eq(s,rank,q.field_name,k,pi,ri,"local.rank",e.this_rank,a.local_block.rank);
 eq(s,rank,q.field_name,k,pi,ri,"local.block",e.this_block,a.local_block.local_block);
 eq(s,rank,q.field_name,k,pi,ri,"peer.rank",e.neighbor_rank,a.peer_block.rank);
 eq(s,rank,q.field_name,k,pi,ri,"peer.block",e.neighbor_block,a.peer_block.local_block);
 eq(s,rank,q.field_name,k,pi,ri,"pattern_order",ri,a.pattern_order);
 eq(s,rank,q.field_name,k,pi,ri,"send_box.has_value",has_send(p.kind),a.send_box.has_value());
 eq(s,rank,q.field_name,k,pi,ri,"recv_box.has_value",has_recv(p.kind),a.recv_box.has_value());
 if(a.send_box)box(s,rank,q.field_name,k,pi,ri,"send_box",e.send_box,*a.send_box);
 if(a.recv_box)box(s,rank,q.field_name,k,pi,ri,"recv_box",e.recv_box,*a.recv_box);
 for(int d=0;d<3;++d) {
  eq(s,rank,q.field_name,k,pi,ri,"transform.permutation["+std::to_string(d)+"]",e.trans.perm[d],a.transform->permutation[d]);
  eq(s,rank,q.field_name,k,pi,ri,"transform.sign["+std::to_string(d)+"]",e.trans.sign[d],a.transform->sign[d]);
 }
 eq(s,rank,q.field_name,k,pi,ri,"transform.offset.i",e.trans.offset.i,a.transform->offset.i);
 eq(s,rank,q.field_name,k,pi,ri,"transform.offset.j",e.trans.offset.j,a.transform->offset.j);
 eq(s,rank,q.field_name,k,pi,ri,"transform.offset.k",e.trans.offset.k,a.transform->offset.k);
 const bool send_tagged=p.kind!=LegacyHaloPatternKind::InnerCopy&&has_send(p.kind);
 const bool recv_tagged=p.kind!=LegacyHaloPatternKind::InnerCopy&&has_recv(p.kind);
 eq(s,rank,q.field_name,k,pi,ri,"send_tag.has_value",send_tagged,a.send_tag.has_value());
 eq(s,rank,q.field_name,k,pi,ri,"recv_tag.has_value",recv_tagged,a.recv_tag.has_value());
 if(a.send_tag)eq(s,rank,q.field_name,k,pi,ri,"send_tag",e.send_flag,*a.send_tag);
 if(a.recv_tag)eq(s,rank,q.field_name,k,pi,ri,"recv_tag",e.recv_flag,*a.recv_tag);
}
template<class F> bool rejected(F f) { try { f(); } catch(const std::invalid_argument&) { return true; } return false; }
}

namespace Z0_TEST {
bool RunExchangePlanTests(const Grid& grid,const TOPO::Topology& topology,Field& field,Halo& halo,int rank) {
 int size=0; PARALLEL::mpi_size(&size);
 const auto partition=compat::LegacyMeshPartitionAdapter::build(grid,topology,rank,size);
 Summary s; long long plans=0,regions=0,unsupported=0;
 for(const auto& q:field.halo_requests()) {
  if(!q.do_halo)continue;
  if(q.orientation_aware) {
   ++unsupported;
   eq(s,rank,q.field_name,ExchangeKind::InnerCopy,0,0,"unsupported.orientation_aware",true,
      rejected([&]{compat::LegacyHaloPlanAdapter::build(partition,halo,q);}));
   continue;
  }
  const auto raw=halo.debug_component_copy_patterns(q.location,q.nghost);
  const auto plan=compat::LegacyHaloPlanAdapter::build(partition,halo,q);
  ++plans;
  eq(s,rank,q.field_name,ExchangeKind::InnerCopy,0,0,"plan.location",(int)loc(q.location),(int)plan.location());
  eq(s,rank,q.field_name,ExchangeKind::InnerCopy,0,0,"plan.ghost_layers",q.nghost,plan.ghost_layers());
  size_t flat=0;
  for(size_t pi=0;pi<raw.size();++pi) {
   eq(s,rank,q.field_name,kind(raw[pi].kind),pi,0,"pattern.location",(int)loc(q.location),(int)loc(raw[pi].location));
   eq(s,rank,q.field_name,kind(raw[pi].kind),pi,0,"pattern.nghost",q.nghost,raw[pi].nghost);
   for(size_t ri=0;ri<raw[pi].regions.size();++ri,++flat) {
    if(flat<plan.regions().size()) region(s,rank,q,pi,ri,raw[pi],raw[pi].regions[ri],plan.regions()[flat]);
    else { ++s.checked; ++s.failed; std::cerr<<"[ExchangePlan][rank "<<rank<<"][request "<<q.field_name<<"] missing flat region "<<flat<<"\n"; }
    ++regions;
   }
  }
  eq(s,rank,q.field_name,ExchangeKind::InnerCopy,0,0,"region_count",flat,plan.regions().size());
 }
 FieldHaloRequest missing; missing.field_name="missing-pattern"; missing.location=StaggerLocation::Cell; missing.nghost=987654;
 eq(s,rank,missing.field_name,ExchangeKind::InnerCopy,0,0,"missing_pattern_rejected",true,
    rejected([&]{compat::LegacyHaloPlanAdapter::build(partition,halo,missing);}));

 parallel::ExchangeRegion valid;
 valid.location=parallel::FieldLocation::Cell; valid.ghost_layers=1;
 valid.stage=parallel::HaloStage::Face; valid.kind=ExchangeKind::InnerCopy;
 valid.local_block={rank,0}; valid.peer_block={rank,0};
 valid.send_box=core::LogicalBox{}; valid.recv_box=core::LogicalBox{};
 valid.transform=mesh::IndexTransform{};
 auto bad_address=valid; bad_address.local_block.local_block=-1;
 eq(s,rank,"contract",ExchangeKind::InnerCopy,0,0,"invalid_block_address_rejected",true,
    rejected([&]{parallel::ExchangePlan(parallel::FieldLocation::Cell,1,{bad_address});}));
 auto missing_transform=valid; missing_transform.transform.reset();
 eq(s,rank,"contract",ExchangeKind::InnerCopy,0,0,"missing_transform_rejected",true,
    rejected([&]{parallel::ExchangePlan(parallel::FieldLocation::Cell,1,{missing_transform});}));
 const auto gc=global_sum(s.checked),gf=global_sum(s.failed),gp=global_sum(plans),gr=global_sum(regions),gu=global_sum(unsupported);
 std::ostringstream d; d<<"plans="<<gp<<" regions="<<gr<<" checks="<<gc<<" failed="<<gf
                         <<" unsupported_orientation_requests="<<gu;
 return print_result("ExchangePlan legacy pattern parity",s.failed==0,d.str());
}
}
