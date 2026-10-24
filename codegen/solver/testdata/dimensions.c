#include "lumo_solver.h"
#include "lumo_prepare.h"
#include "lumo_initialguess.h"
#include "ctest.h"
static double workspace[LUMO_SOLVER_WORKSPACE_DOUBLES];
static int dimensions(void){
 struct lumo_prepare_input input={0};struct lumo_mesh_operating op;struct lumo_solver_report report;
 double state[LUMO_MESH_SIZE],r[LUMO_MESH_SIZE],balance[LUMO_MESH_COMPONENTS]={0};size_t f,i,j;
 for(f=0;f<LUMO_MESH_FEEDS;++f){input.feeds[f].flow=100;input.feeds[f].fractions[0]=0.5;input.feeds[f].fractions[1]=0.5;
#if LUMO_MESH_COMPONENTS == 3
 input.feeds[f].fractions[1]=0.3;input.feeds[f].fractions[2]=0.2;
#endif
 }
 input.distillate_flow=50;input.reflux_ratio=2.5;
 CTEST_STATUS(lumo_prepare_operating(&input,&op),0);CTEST_STATUS(lumo_initialguess_linear_T_x(&op,state),0);
 CTEST_STATUS(lumo_solver_solve(&op,state,NULL,workspace,LUMO_SOLVER_WORKSPACE_DOUBLES,&report),LUMO_SOLVER_CONVERGED);
 CTEST_STATUS(lumo_mesh_residual(state,&op,r),0);
 for(j=0;j<LUMO_MESH_STAGES;++j){double sx=0,sy=0;for(i=0;i<LUMO_MESH_COMPONENTS;++i){sx+=state[j*LUMO_MESH_BLOCK+LUMO_MESH_X+i];sy+=state[j*LUMO_MESH_BLOCK+LUMO_MESH_Y+i];balance[i]+=r[j*LUMO_MESH_BLOCK+LUMO_MESH_M+i];}CTEST_NEAR(sx,1);CTEST_NEAR(sy,1);}
 for(i=0;i<LUMO_MESH_COMPONENTS;++i)CTEST_STATUS(fabs(balance[i]/(100*LUMO_MESH_FEEDS))<1e-9,1);
 return 0;
}
int main(int argc,char **argv){const struct ctest_case cases[]={{"dimensions",dimensions}};return ctest_run(argc,argv,cases,1);}
