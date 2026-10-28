#include "9_runtime/SimulationContext.h"
#include "2_topology/TopologyBuilder.h"
#include <algorithm>
#include <set>
#include <stdexcept>
#include <utility>
namespace
{
void byte(std::uint64_t &h, unsigned char v) { h ^= v; h *= 1099511628211ULL; }
template <typename T> void value(std::uint64_t &h, T v)
{ auto raw = static_cast<std::uint64_t>(v); for (unsigned i=0;i<sizeof(raw);++i) byte(h,(raw>>(8U*i))&0xffU); }
void text(std::uint64_t &h, const std::string &v) { value(h,v.size()); for(auto c:v) byte(h,c); }
bool includes(HaloLevel a, HaloLevel b) { return static_cast<int>(a)>=static_cast<int>(b); }
}
namespace MPCNS
{
void MetricAccess::recompute()
{ if(!enabled_||!field_||!grid_) throw std::logic_error("metric fields disabled"); METRIC::compute_metric_fields(*field_,*grid_); }
METRIC::MetricDiagnostics MetricAccess::diagnostics()
{ if(!enabled_||!field_) throw std::logic_error("metric fields disabled"); return METRIC::diagnose_metric_fields(*field_); }
void DecAccess::d0(const std::string &n,const DEC::EdgeFormNames &e,int nc,int ec)
{ DEC::d0_node_to_edge(*field_,n,e,nc,ec); }
void DecAccess::d1(const DEC::EdgeFormNames &e,const DEC::FaceFormNames &f,int ec,int fc)
{ DEC::d1_edge_to_face(*field_,e,f,ec,fc); }
void DecAccess::d2(const DEC::FaceFormNames &f,const std::string &c,int fc,int cc)
{ DEC::d2_face_to_cell(*field_,f,c,fc,cc); }
void DecAccess::faraday(const DEC::EdgeFormNames &e,const DEC::FaceFormNames &f,double dt,int ec,int fc)
{ DEC::faraday_update_face_2form(*field_,e,f,dt,ec,fc); }
SimulationContext::SimulationContext(Param &p,int r,SimulationConfig c)
 : param_(p),rank_(r),config_(std::move(c)) {}
SimulationContext::~SimulationContext()=default;
void SimulationContext::require_initialized_(const char *op) const
{ if(!initialized_) throw std::logic_error(std::string("SimulationContext::")+op+" before initialize"); }
void SimulationContext::validate_config_() const
{
 std::set<std::string> names;
 for(const auto &d:config_.fields)
 {
  if(d.name.empty()||d.ncomp<=0||d.nghost<0) throw std::invalid_argument("invalid field descriptor");
  if(!names.insert(d.name).second) throw std::invalid_argument("duplicate field: "+d.name);
  if((d.sync.do_halo||d.sync.do_physical||d.sync.do_coupling||d.sync.owner_sync!=OwnerSyncPolicy::None)&&d.sync.group.empty())
   throw std::invalid_argument("synchronized field has empty group: "+d.name);
 }
 for(const auto &p:config_.coupling_pairs) if(p.first.empty()||p.second.empty()||p.first==p.second)
  throw std::invalid_argument("invalid coupling pair");
 auto check=[&](const auto &v){for(const auto &n:v)if(!names.count(n))throw std::invalid_argument("unknown IO field: "+n);};
 check(config_.io.output_fields); check(config_.io.restart_fields);
}
void SimulationContext::validate_collective_configuration_() const
{
 std::uint64_t h=1469598103934665603ULL; value(h,config_.fields.size());
 value(h,config_.geometry_ghost>=0?config_.geometry_ghost:param_.GetInt("ngg"));
 value(h,config_.dimension>=0?config_.dimension:param_.GetInt("dimension"));
 value(h,(int)config_.physical_corner_mode);value(h,config_.io.enabled);
 value(h,config_.io.residual_components);value(h,config_.io.output_fields.size());
 for(const auto &n:config_.io.output_fields)text(h,n);
 value(h,config_.io.restart_fields.size());for(const auto &n:config_.io.restart_fields)text(h,n);
 for(const auto &d:config_.fields)
 {
  text(h,d.name); value(h,(int)d.location); value(h,d.ncomp); value(h,d.nghost); text(h,d.physics);
  text(h,d.sync.group); value(h,d.sync.do_coupling); value(h,d.sync.do_physical); value(h,d.sync.do_halo);
  value(h,(int)d.sync.halo_level); value(h,(int)d.sync.owner_sync); value(h,d.sync.orientation_aware);
  value(h,(int)d.value_kind); value(h,(int)d.storage_role);
 }
 value(h,config_.coupling_pairs.size()); for(const auto &p:config_.coupling_pairs){text(h,p.first);text(h,p.second);}
 value(h,config_.enable_metric_fields);
 unsigned long long local=h,low=0,high=0;
 MPI_Allreduce(&local,&low,1,MPI_UNSIGNED_LONG_LONG,MPI_MIN,MPI_COMM_WORLD);
 MPI_Allreduce(&local,&high,1,MPI_UNSIGNED_LONG_LONG,MPI_MAX,MPI_COMM_WORLD);
 if(low!=high) throw std::runtime_error("SimulationContext configuration differs across ranks");
}
void SimulationContext::initialize()
{
 if(initialized_) throw std::logic_error("initialize called twice");
 int active=0,finalized=0; MPI_Initialized(&active); MPI_Finalized(&finalized);
 if(!active||finalized) throw std::logic_error("active caller-owned MPI required");
 int actual_rank=-1;MPI_Comm_rank(MPI_COMM_WORLD,&actual_rank);
 if(actual_rank!=rank_)throw std::invalid_argument("SimulationContext rank does not match MPI_COMM_WORLD");
 if(config_.dimension>=0&&config_.dimension!=param_.GetInt("dimension"))
  throw std::invalid_argument("SimulationContext dimension differs from Grid Param");
 validate_config_(); validate_collective_configuration_();
 int ghost=config_.geometry_ghost>=0?config_.geometry_ghost:param_.GetInt("ngg");
 int dim=config_.dimension>=0?config_.dimension:param_.GetInt("dimension");
 grid_.Grid_Preprocess(&param_);
 topology_=std::make_unique<TOPO::Topology>(TOPO::build_topology(grid_,rank_,dim));
 if(config_.configure_topology) config_.configure_topology(grid_,*topology_);
 field_=std::make_unique<Field>(&grid_,&param_,ghost);
 if(config_.enable_metric_fields) METRIC::register_metric_fields(*field_,ghost);
 for(const auto &d:config_.fields) field_->register_field(d);
 if(config_.enable_metric_fields) METRIC::compute_metric_fields(*field_,grid_);
 field_->register_declared_coupling_channels(config_.coupling_pairs);
 field_->build_coupling_buffers(*topology_,dim);
 halo_=std::make_unique<Halo>(field_.get(),topology_.get());
 halo_->set_topology_equiv(topology_.get());
 halo_->register_halo_fields(field_->halo_requests());
 halo_->build_registered_patterns(); prepare_oriented_exchanges_();
 boundary_=std::make_unique<BoundaryCore>();
 boundary_->SetUp(&grid_,field_.get(),topology_.get(),&param_,field_->boundary_field_names());
 if(config_.io.enabled)
 {
  io_=std::make_unique<IOModule>();
  io_->Setup(&param_,&grid_,field_.get(),std::max(0,config_.io.residual_components));
  io_->SetTecplotFields(config_.io.output_fields); io_->SetParaViewFields(config_.io.output_fields);
  io_->SetRestartFields(config_.io.restart_fields);
 }
 initialized_=true; ++communication_generation_;
}
void SimulationContext::prepare_oriented_exchanges_()
{
 bool e2=false,e3=false,f2=false,f3=false;
 for(const auto &d:config_.fields)
 {
  if(!d.sync.do_halo||!d.sync.orientation_aware) continue;
  bool at2=includes(d.sync.halo_level,HaloLevel::Corner2D),at3=includes(d.sync.halo_level,HaloLevel::Corner3D);
  if(d.value_kind==FieldValueKind::EdgeCovariant1Form){e2|=at2;e3|=at3;}
  if(d.value_kind==FieldValueKind::FaceContravariant2Form){f2|=at2;f3|=at3;}
 }
 if(grid_.dimension>=2&&e2&&(!halo_->prepare_edge_1form_corner_metadata()||!halo_->activate_prepared_edge_1form_edge_execution())) throw std::runtime_error("Edge1Form edge prepare failed");
 if(grid_.dimension>=3&&e3&&(!halo_->prepare_vertex_edge_1form_metadata()||!halo_->activate_prepared_edge_1form_vertex_execution())) throw std::runtime_error("Edge1Form vertex prepare failed");
 if(grid_.dimension>=2&&f2&&(!halo_->prepare_edge_face_2form_metadata()||!halo_->activate_prepared_face_2form_edge_execution())) throw std::runtime_error("Face2Form edge prepare failed");
 if(grid_.dimension>=3&&f3&&(!halo_->prepare_vertex_face_2form_metadata()||!halo_->activate_prepared_face_2form_vertex_execution())) throw std::runtime_error("Face2Form vertex prepare failed");
}
void SimulationContext::register_physical_handler(const std::string &f,const std::string &bc,BOUND::PhysicalHandler h)
{require_initialized_("register_physical_handler");boundary_->RegisterPhysical(f,bc,std::move(h));}
void SimulationContext::register_coupling_handler(const std::string &s,const std::string &d,StaggerLocation l,const std::string &t,const std::string &f,BOUND::CouplingHandler h)
{require_initialized_("register_coupling_handler");boundary_->RegisterCoupling(s,d,l,t,f,std::move(h));}
void SimulationContext::validate_physical_handlers() const
{require_initialized_("validate_physical_handlers");boundary_->CheckPhysicalHandlers(field_->boundary_field_names());}
SimulationContext::SyncTarget SimulationContext::resolve_target_(const std::string &target,bool group_only) const
{
 SyncTarget out;out.group=target;
 if(!group_only&&field_->has_field(target))
 {
  const auto &d=field_->descriptor(target);out.by_field=true;out.group=d.sync.group.empty()?target:d.sync.group;
  if(d.sync.orientation_aware&&!d.sync.group.empty()) for(const auto &x:field_->descriptors()) if(x.sync.group==d.sync.group) out.fields.push_back(x.name);
  if(out.fields.empty()) out.fields.push_back(target);
 }
 else for(const auto &d:field_->descriptors()) if(d.sync.group==target) out.fields.push_back(d.name);
 if(out.fields.empty()) throw std::invalid_argument("unknown sync target: "+target);
 return out;
}
HaloLevel SimulationContext::target_level_(const SyncTarget &t,const SyncOptions &o) const
{
 HaloLevel l=HaloLevel::Corner1D;for(const auto &n:t.fields){auto c=field_->descriptor(n).sync.halo_level;if((int)c>(int)l)l=c;}
 if(o.through&&(int)*o.through<(int)l)l=*o.through;return l;
}
std::vector<std::string> SimulationContext::physical_fields_(const SyncTarget &t) const
{std::vector<std::string> out;for(const auto &n:t.fields)if(field_->descriptor(n).sync.do_physical)out.push_back(n);return out;}
std::vector<int32_t> SimulationContext::coupling_ids_(const SyncTarget &t,const Field::PairKey &p) const
{
 std::set<std::string> selected;for(const auto &n:t.fields)if(field_->descriptor(n).sync.do_coupling)selected.insert(n);
 std::vector<int32_t> out;if(selected.empty()||!field_->has_coupling_pair(p.first,p.second))return out;
 const auto &c=field_->coupling_pair(p.first,p.second).channels;
 for(int32_t i=0;i<(int32_t)c.size();++i)if(selected.count(c[(size_t)i].tag))out.push_back(i);return out;
}
void SimulationContext::apply_physical_(const std::vector<std::string> &f,HaloLevel s)
{
 if(f.empty())return;
 if(s==HaloLevel::Corner1D){boundary_->ApplyPhysical(f);if(config_.physical_corner_mode==PhysicalCornerMode::MercuryCompatible)boundary_->ApplyPhysicalCornerDefault(f);}
 else if(s==HaloLevel::Corner2D&&config_.physical_corner_mode!=PhysicalCornerMode::SolverManaged)boundary_->ApplyPhysicalEdgeDefault(f);
 else if(s==HaloLevel::Corner3D&&config_.physical_corner_mode==PhysicalCornerMode::StageLocal)boundary_->ApplyPhysicalVertexDefault(f);
}
void SimulationContext::sync_halo_(const SyncTarget &t,HaloLevel s)
{
 bool any=false;for(const auto &n:t.fields)any|=field_->descriptor(n).sync.do_halo;if(!any)return;
 if(t.by_field&&t.fields.size()==1)halo_->sync_field(t.fields.front(),s);else halo_->sync_group(t.group,s);
}
void SimulationContext::sync_coupling_(const SyncTarget &t,HaloLevel stage)
{
 for(const auto &p:config_.coupling_pairs)
 {
  auto ids=coupling_ids_(t,p);if(ids.empty())continue;std::string src=p.first,dst=p.second;
  if(stage==HaloLevel::Corner1D)halo_->coupling_trans_1DCorner(src,dst,ids);
  else if(stage==HaloLevel::Corner2D)halo_->coupling_trans_2DCorner(src,dst,ids);
  else halo_->coupling_trans_3DCorner(src,dst,ids);
  boundary_->ApplyCouplingPair(src,dst,stage,ids);
 }
}
void SimulationContext::sync_owner_(const SyncTarget &t)
{
 bool any=false;for(const auto &n:t.fields)any|=field_->descriptor(n).sync.owner_sync!=OwnerSyncPolicy::None;
 if(!any)return;if(t.by_field&&t.fields.size()==1)halo_->sync_owner_alias_field(t.fields.front());else halo_->sync_owner_alias_group(t.group);
}
void SimulationContext::sync_resolved_(const SyncTarget &t,const SyncOptions &o)
{
 auto max=target_level_(t,o);auto phy=physical_fields_(t);
 for(auto s:{HaloLevel::Corner1D,HaloLevel::Corner2D,HaloLevel::Corner3D})
 {if(!includes(max,s))continue;if(o.physical)apply_physical_(phy,s);if(o.halo)sync_halo_(t,s);if(o.coupling)sync_coupling_(t,s);}
 if(o.owner_alias)sync_owner_(t);
}
void SimulationContext::sync(const std::string &t){sync(t,SyncOptions{});}
void SimulationContext::sync(const std::string &t,const SyncOptions &o){require_initialized_("sync");sync_resolved_(resolve_target_(t,false),o);}
void SimulationContext::sync_group(const std::string &g){sync_group(g,SyncOptions{});}
void SimulationContext::sync_group(const std::string &g,const SyncOptions &o){require_initialized_("sync_group");sync_resolved_(resolve_target_(g,true),o);}
Grid &SimulationContext::grid(){require_initialized_("grid");return grid_;}
const Grid &SimulationContext::grid()const{require_initialized_("grid");return grid_;}
TOPO::Topology &SimulationContext::topology(){require_initialized_("topology");return *topology_;}
const TOPO::Topology &SimulationContext::topology()const{require_initialized_("topology");return *topology_;}
Field &SimulationContext::field(){require_initialized_("field");return *field_;}
const Field &SimulationContext::field()const{require_initialized_("field");return *field_;}
Halo &SimulationContext::halo(){require_initialized_("halo");return *halo_;}
const Halo &SimulationContext::halo()const{require_initialized_("halo");return *halo_;}
BoundaryCore &SimulationContext::boundary(){require_initialized_("boundary");return *boundary_;}
const BoundaryCore &SimulationContext::boundary()const{require_initialized_("boundary");return *boundary_;}
IOModule &SimulationContext::io(){require_initialized_("io");if(!io_)throw std::logic_error("IO disabled");return *io_;}
const IOModule &SimulationContext::io()const{require_initialized_("io");if(!io_)throw std::logic_error("IO disabled");return *io_;}
RunData &SimulationContext::run_data(){return io().Run();}
const RunData &SimulationContext::run_data()const{return io().Run();}
MetricAccess SimulationContext::metric(){require_initialized_("metric");return MetricAccess(field_.get(),&grid_,config_.enable_metric_fields);}
DecAccess SimulationContext::dec(){require_initialized_("dec");return DecAccess(field_.get());}
} // namespace MPCNS
