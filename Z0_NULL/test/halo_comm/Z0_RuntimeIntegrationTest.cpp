#include "Z0_Tests.h"
#include "Z0_TestCommon.h"
#include "9_runtime/SimulationContext.h"
#include "2_topology/TopologyBuildDetail.h"
#include <cstdlib>
#include <iostream>
#include <set>
namespace
{
const std::string src="RUNTIME_SOURCE",dst="RUNTIME_DESTINATION";
bool enabled(){const char *v=std::getenv("Z0_RUNTIME_TEST");return v&&std::string(v)!="0";}
FieldDescriptor desc(const std::string &n,StaggerLocation l,FieldValueKind k,
 const std::string &g,OwnerSyncPolicy o,bool oriented,bool physical=false,bool coupling=true)
{
 FieldDescriptor d;d.name=n;d.location=l;d.value_kind=k;d.ncomp=1;d.nghost=1;
 d.sync.group=g;d.sync.do_halo=true;d.sync.do_physical=physical;d.sync.do_coupling=coupling;
 d.sync.halo_level=HaloLevel::Corner3D;d.sync.owner_sync=o;d.sync.orientation_aware=oriented;return d;
}
MPCNS::SimulationConfig config()
{
 MPCNS::SimulationConfig c;c.geometry_ghost=1;c.io.residual_components=2;c.enable_metric_fields=true;
 auto cell=desc("RuntimeCell",StaggerLocation::Cell,FieldValueKind::CartesianVector,"RuntimeCell",OwnerSyncPolicy::None,false,true);cell.ncomp=3;c.fields.push_back(cell);
 c.fields.push_back(desc("RuntimeNode",StaggerLocation::Node,FieldValueKind::Scalar,"RuntimeNode",OwnerSyncPolicy::NodeOwner,false,false,false));
 c.fields.push_back(desc("RuntimeE_xi",StaggerLocation::EdgeXi,FieldValueKind::EdgeCovariant1Form,"RuntimeE",OwnerSyncPolicy::EdgeOwner,true));
 c.fields.push_back(desc("RuntimeE_eta",StaggerLocation::EdgeEt,FieldValueKind::EdgeCovariant1Form,"RuntimeE",OwnerSyncPolicy::EdgeOwner,true));
 c.fields.push_back(desc("RuntimeE_zeta",StaggerLocation::EdgeZe,FieldValueKind::EdgeCovariant1Form,"RuntimeE",OwnerSyncPolicy::EdgeOwner,true));
 c.fields.push_back(desc("RuntimeB_xi",StaggerLocation::FaceXi,FieldValueKind::FaceContravariant2Form,"RuntimeB",OwnerSyncPolicy::FaceOwner,true));
 c.fields.push_back(desc("RuntimeB_eta",StaggerLocation::FaceEt,FieldValueKind::FaceContravariant2Form,"RuntimeB",OwnerSyncPolicy::FaceOwner,true));
 c.fields.push_back(desc("RuntimeB_zeta",StaggerLocation::FaceZe,FieldValueKind::FaceContravariant2Form,"RuntimeB",OwnerSyncPolicy::FaceOwner,true));
 c.coupling_pairs={{src,dst},{dst,src}};c.io.output_fields={"RuntimeCell"};c.io.restart_fields={"RuntimeCell"};
 c.configure_topology=[](Grid &grid,TOPO::Topology &topology)
 {
  TOPO::PhysicalPatch wall;wall.this_rank=Z0_TEST::rank();wall.this_block=0;
  wall.this_block_name=grid.grids(0).block_name;wall.bc_id=9001;wall.bc_name="RuntimeWall";
  wall.this_box_node=Box3{{0,0,0},{1,grid.grids(0).my+1,grid.grids(0).mz+1}};wall.direction=-1;
  topology.physical_patches.push_back(wall);
  for(auto &p:topology.inner_patches){p.is_coupling=p.this_block!=p.nb_block&&p.this_rank!=2;if(!p.is_coupling)continue;bool f=p.this_block<p.nb_block;p.this_block_name=f?src:dst;p.nb_block_name=f?dst:src;}
  for(auto &p:topology.parallel_patches){p.is_coupling=(p.this_rank==0&&p.nb_rank==1)||(p.this_rank==1&&p.nb_rank==0);if(!p.is_coupling)continue;bool f=p.this_rank==0;p.this_block_name=f?src:dst;p.nb_block_name=f?dst:src;}
  TOPO::detail::build_edge_patches(grid,topology,grid.dimension);TOPO::detail::build_vertex_patches(grid,topology,grid.dimension);
 };return c;
}
void fill(Field &field,const std::vector<FieldDescriptor> &ds)
{
 for(const auto &d:ds)for(int b=0;b<field.num_blocks();++b){auto &u=field.field(d.name,b);if(!u.is_allocated())continue;auto lo=u.get_lo(),hi=u.get_hi();
 for(int i=lo.i;i<hi.i;++i)for(int j=lo.j;j<hi.j;++j)for(int k=lo.k;k<hi.k;++k)for(int m=0;m<d.ncomp;++m)u(i,j,k,m)=Z0_TEST::unique_code(Z0_TEST::rank(),b,i,j,k,m);}
}
std::vector<double> snap(Field &f,const std::string &n)
{
 std::vector<double> v;auto &d=f.descriptor(n);for(int b=0;b<f.num_blocks();++b){auto &u=f.field(n,b);if(!u.is_allocated())continue;auto lo=u.get_lo(),hi=u.get_hi();
 for(int i=lo.i;i<hi.i;++i)for(int j=lo.j;j<hi.j;++j)for(int k=lo.k;k<hi.k;++k)for(int m=0;m<d.ncomp;++m)v.push_back(u(i,j,k,m));}return v;
}
void restore(Field &f,const std::string &n,const std::vector<double> &v)
{
 size_t q=0;auto &d=f.descriptor(n);for(int b=0;b<f.num_blocks();++b){auto &u=f.field(n,b);if(!u.is_allocated())continue;auto lo=u.get_lo(),hi=u.get_hi();
 for(int i=lo.i;i<hi.i;++i)for(int j=lo.j;j<hi.j;++j)for(int k=lo.k;k<hi.k;++k)for(int m=0;m<d.ncomp;++m)u(i,j,k,m)=v.at(q++);}
}
}
namespace Z0_TEST
{
bool RunRuntimeIntegrationTests(Param &param,int rank)
{
 if(!enabled())return true;if(rank==0)std::cout<<"Z0 Solver Runtime integration tests\n";
 auto c=config();auto declared=c.fields;MPCNS::SimulationContext sim(param,rank,std::move(c));sim.initialize();bool ok=true;
 ok&=print_result("Runtime dependency initialization",sim.initialized()&&sim.communication_generation()==1&&sim.field().has_field("RuntimeCell")&&sim.halo().topology_equiv()==&sim.topology(),"generation=1");
 std::set<std::string> bcs;for(const auto &p:sim.topology().physical_patches)bcs.insert(p.bc_name);
 long long pcalls=0,ccalls=0;
 for(const auto &bc:bcs)sim.register_physical_handler("RuntimeCell",bc,[&](FieldBlock &u,Field *f,const BOUND::PhysicalRegion &r,int n){++pcalls;BoundaryCore::DefaultPhysicalCopy(u,f,r,n);});
 sim.validate_physical_handlers();
 for(const auto &pair:sim.field().coupling_pairs())for(const auto &ch:pair.second.channels)
  sim.register_coupling_handler(pair.first.first,pair.first.second,ch.location,ch.tag,ch.tag,
   [&](FieldBlock &u,Field *f,CouplingBufferBlock &b,const std::string &s,const std::string &d,const std::string &tag){++ccalls;BoundaryCore::DefaultCouplingCopy(u,f,b,s,d,tag);});
 fill(sim.field(),declared);auto initial=snap(sim.field(),"RuntimeCell");
 MPCNS::SyncOptions only;only.physical=false;only.coupling=false;only.owner_alias=false;
 sim.sync("RuntimeCell",only);auto runtime=snap(sim.field(),"RuntimeCell");restore(sim.field(),"RuntimeCell",initial);
 for(auto stage:{HaloLevel::Corner1D,HaloLevel::Corner2D,HaloLevel::Corner3D})sim.halo().sync_field("RuntimeCell",stage);
 auto legacy=snap(sim.field(),"RuntimeCell");
 ok&=print_result("Runtime sync delegates exactly to staged Halo",runtime==legacy,"values="+std::to_string(runtime.size()));
 fill(sim.field(),declared);sim.dec().d0("RuntimeNode",DEC::EdgeFormNames{"RuntimeE_xi","RuntimeE_eta","RuntimeE_zeta"});
 sim.dec().d1(DEC::EdgeFormNames{"RuntimeE_xi","RuntimeE_eta","RuntimeE_zeta"},DEC::FaceFormNames{"RuntimeB_xi","RuntimeB_eta","RuntimeB_zeta"});
 auto diagnostics=sim.metric().diagnostics();(void)diagnostics;sim.run_data().step=7;
 ok&=print_result("Runtime Metric DEC IO access",sim.metric().enabled()&&sim.run_data().step==7&&sim.output_fields()==std::vector<std::string>{"RuntimeCell"}&&sim.restart_fields()==std::vector<std::string>{"RuntimeCell"},"");
 sim.sync("RuntimeCell");sim.sync_group("RuntimeE");sim.sync_group("RuntimeB");sim.sync("RuntimeNode");
 auto stable=sim.halo().debug_prepared_exchange_statistics();auto cbefore=sim.halo().debug_coupling_execution_counters();
 sim.sync("RuntimeCell");sim.sync_group("RuntimeE");sim.sync_group("RuntimeB");sim.sync("RuntimeNode");
 auto repeated=sim.halo().debug_prepared_exchange_statistics();auto cafter=sim.halo().debug_coupling_execution_counters();
 long long gp=global_sum(pcalls),gc=global_sum(ccalls);
 ok&=print_result("Runtime handlers and declarative groups execute",gp>0&&gc>0,"physical="+std::to_string(gp)+" coupling="+std::to_string(gc));
 ok&=print_result("Runtime stable sync performs no metadata rebuild or workspace growth",stable.workspace_growth_events==repeated.workspace_growth_events&&cbefore.initialization_collectives==cafter.initialization_collectives&&sim.communication_generation()==1,"growth="+std::to_string(repeated.workspace_growth_events));
 ok&=print_result("Runtime oriented plans are eagerly active",sim.halo().debug_prepared_edge_1form_edge_execution_active()&&sim.halo().debug_prepared_edge_1form_vertex_execution_active()&&sim.halo().debug_prepared_face_2form_edge_execution_active()&&sim.halo().debug_prepared_face_2form_vertex_execution_active(),"");
 return ok;
}
}
