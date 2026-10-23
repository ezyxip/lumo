#include "lumo_initialguess.h"
#include "ctest.h"

static struct lumo_mesh_operating operating(void){
 struct lumo_mesh_operating op={0};size_t f;
 op.distillate_mol_s=50;op.reflux_ratio=2.5;
 for(f=0;f<LUMO_MESH_FEEDS;++f){op.feeds[f].flow_mol_s=100;op.feeds[f].mole_fractions[0]=0.5;op.feeds[f].mole_fractions[1]=0.5;op.feeds[f].enthalpy_J_mol=2060.414799029395;}
 return op;
}
static int verify(const struct lumo_mesh_operating *op,const double *state){
 double r[LUMO_MESH_SIZE],total=0;size_t i,j,f;
 CTEST_STATUS(lumo_mesh_residual(state,op,r),LUMO_MESH_OK);
 for(f=0;f<LUMO_MESH_FEEDS;++f)total+=op->feeds[f].flow_mol_s;
 for(j=0;j<LUMO_MESH_STAGES;++j){
  size_t b=j*LUMO_MESH_BLOCK;double sx=0,sy=0;
  for(i=0;i<LUMO_MESH_COMPONENTS;++i){
   CTEST_STATUS(state[b+LUMO_MESH_X+i]>=0,1);CTEST_STATUS(state[b+LUMO_MESH_Y+i]>=0,1);
   sx+=state[b+LUMO_MESH_X+i];sy+=state[b+LUMO_MESH_Y+i];
   if(j>0&&j+1<LUMO_MESH_STAGES){
    double t=(double)j/(LUMO_MESH_STAGES-1);
    CTEST_NEAR(state[b+LUMO_MESH_X+i],(1-t)*state[LUMO_MESH_X+i]+t*state[(LUMO_MESH_STAGES-1)*LUMO_MESH_BLOCK+LUMO_MESH_X+i]);
   }
  }
  CTEST_NEAR(sx,1);CTEST_NEAR(sy,1);
  CTEST_STATUS(state[b+LUMO_MESH_T]>=350&&state[b+LUMO_MESH_T]<=390,1);
  if(j>0&&j+1<LUMO_MESH_STAGES)CTEST_NEAR(2*state[b+LUMO_MESH_T],state[b-LUMO_MESH_BLOCK+LUMO_MESH_T]+state[b+LUMO_MESH_BLOCK+LUMO_MESH_T]);
  CTEST_STATUS(state[b+LUMO_MESH_L]>=0,1);
  if(j>0)CTEST_STATUS(state[b+LUMO_MESH_V]>0,1);
  {double sum=0;for(i=0;i<LUMO_MESH_COMPONENTS;++i)sum+=r[b+LUMO_MESH_M+i];CTEST_NEAR(sum/total,0);}
 }
 CTEST_NEAR(state[(LUMO_MESH_STAGES-1)*LUMO_MESH_BLOCK+LUMO_MESH_L],total-op->distillate_mol_s);
 CTEST_NEAR(r[LUMO_MESH_CLOSURE],0);
 CTEST_NEAR(r[LUMO_MESH_H]/1e7,0);CTEST_NEAR(r[(LUMO_MESH_STAGES-1)*LUMO_MESH_BLOCK+LUMO_MESH_H]/1e7,0);
 for(i=0;i<LUMO_MESH_COMPONENTS;++i){
  double feed=0;for(f=0;f<LUMO_MESH_FEEDS;++f)feed+=op->feeds[f].flow_mol_s*op->feeds[f].mole_fractions[i];
  CTEST_NEAR((op->distillate_mol_s*state[LUMO_MESH_X+i]+(total-op->distillate_mol_s)*state[(LUMO_MESH_STAGES-1)*LUMO_MESH_BLOCK+LUMO_MESH_X+i]-feed)/total,0);
 }
 return 0;
}
static int profiles(void){
 struct lumo_mesh_operating op=operating();double state[LUMO_MESH_SIZE],again[LUMO_MESH_SIZE];
 lumo_initialguess_fn initialize=lumo_initialguess_linear_T_x;
 CTEST_STATUS(initialize(&op,state),LUMO_THERMO_OK);CTEST_STATUS(verify(&op,state),0);
 CTEST_STATUS(initialize(&op,again),LUMO_THERMO_OK);CTEST_STATUS(memcmp(state,again,sizeof state),0);
 CTEST_STATUS(state[LUMO_MESH_X]>0.5,1);
 return 0;
}
static int zero_reflux(void){struct lumo_mesh_operating op=operating();double state[LUMO_MESH_SIZE];op.reflux_ratio=0;CTEST_STATUS(lumo_initialguess_linear_T_x(&op,state),0);CTEST_NEAR(state[LUMO_MESH_L],0);return verify(&op,state);}
static int pure(void){struct lumo_mesh_operating op=operating();double state[LUMO_MESH_SIZE];size_t f;for(f=0;f<LUMO_MESH_FEEDS;++f){op.feeds[f].mole_fractions[0]=1;op.feeds[f].mole_fractions[1]=0;}CTEST_STATUS(lumo_initialguess_linear_T_x(&op,state),0);return verify(&op,state);}
static int enthalpy(void){struct lumo_mesh_operating op=operating();double state[LUMO_MESH_SIZE];op.feeds[0].enthalpy_J_mol=40000;CTEST_STATUS(lumo_initialguess_linear_T_x(&op,state),0);return verify(&op,state);}
static int distillate_limits(void){
 struct lumo_mesh_operating op=operating();double state[LUMO_MESH_SIZE];
 const double values[]={1e-5,1,99.99,100*LUMO_MESH_FEEDS-1e-8};size_t i;
 for(i=0;i<sizeof values/sizeof values[0];++i){op.distillate_mol_s=values[i];CTEST_STATUS(lumo_initialguess_linear_T_x(&op,state),0);CTEST_STATUS(verify(&op,state),0);}
 return 0;
}
static int errors(void){
 struct lumo_mesh_operating op;double state[LUMO_MESH_SIZE],before[LUMO_MESH_SIZE];int test;
 memset(state,0x34,sizeof state);memcpy(before,state,sizeof state);
 CTEST_STATUS(lumo_initialguess_linear_T_x(NULL,state),LUMO_THERMO_INVALID_INPUT);
 op=operating();CTEST_STATUS(lumo_initialguess_linear_T_x(&op,NULL),LUMO_THERMO_INVALID_INPUT);
 for(test=0;test<8;++test){op=operating();switch(test){
 case 0:op.distillate_mol_s=100*LUMO_MESH_FEEDS;break;case 1:op.reflux_ratio=-1;break;
 case 2:op.feeds[0].flow_mol_s=0;break;case 3:op.feeds[0].mole_fractions[0]=NAN;break;
 case 4:op.feeds[0].mole_fractions[0]=0.6;break;case 5:op.feeds[0].enthalpy_J_mol=INFINITY;break;
 case 6:op.distillate_mol_s=NAN;break;case 7:op.feeds[0].mole_fractions[0]=-0.5;break;}
 CTEST_STATUS(lumo_initialguess_linear_T_x(&op,state),LUMO_THERMO_INVALID_INPUT);CTEST_STATUS(memcmp(state,before,sizeof state),0);}
 return 0;
}
int main(int argc,char **argv){const struct ctest_case cases[]={{"profiles",profiles},{"zero_reflux",zero_reflux},{"pure",pure},{"errors",errors},{"enthalpy",enthalpy},{"distillate_limits",distillate_limits}};return ctest_run(argc,argv,cases,sizeof cases/sizeof cases[0]);}
