#include "lumo_solver.h"
#include "lumo_prepare.h"
#include "lumo_initialguess.h"
#include "ctest.h"
static double workspace[LUMO_SOLVER_WORKSPACE_DOUBLES];
#include "reference.h"
#define N LUMO_MESH_SIZE

static int setup(struct lumo_mesh_operating *op,double *state){
 const struct lumo_prepare_input input={.feeds={{.flow=360,.flow_unit=LUMO_FLOW_KMOL_H,
 .composition_basis=LUMO_COMPOSITION_MASS,.fractions={78.11/(78.11+92.14),92.14/(78.11+92.14)},.state=LUMO_FEED_SATURATED_LIQUID}},
 .distillate_flow=180,.distillate_unit=LUMO_FLOW_KMOL_H,.reflux_ratio=2.5};
 CTEST_STATUS(lumo_prepare_operating(&input,op),LUMO_PREPARE_OK);
 CTEST_STATUS(lumo_initialguess_linear_T_x(op,state),LUMO_THERMO_OK);return 0;
}
static double scale(size_t i){if(i==LUMO_MESH_QD||i==LUMO_MESH_QW)return 1e7;if(i%LUMO_MESH_BLOCK==LUMO_MESH_T)return 350;if(i%LUMO_MESH_BLOCK==LUMO_MESH_L||i%LUMO_MESH_BLOCK==LUMO_MESH_V)return 100;return 1;}
static int compare_reference(const double *state){size_t i;for(i=0;i<N;++i)CTEST_STATUS(fabs((state[i]-reference_state[i])/scale(i))<1e-8,1);return 0;}
static int reference(void){
 struct lumo_mesh_operating op;double state[N];struct lumo_solver_report report;
 CTEST_STATUS(setup(&op,state),0);CTEST_STATUS(lumo_solver_solve(&op,state,NULL,workspace,LUMO_SOLVER_WORKSPACE_DOUBLES,&report),LUMO_SOLVER_CONVERGED);
 CTEST_STATUS(report.status,LUMO_SOLVER_CONVERGED);CTEST_STATUS(report.iterations>0,1);CTEST_STATUS(report.residual_norm<=1e-10,1);return compare_reference(state);
}
static int already_converged(void){
 struct lumo_mesh_operating op;double state[N],before[N];struct lumo_solver_report report;struct lumo_solver_options options;
 CTEST_STATUS(setup(&op,state),0);memcpy(state,reference_state,sizeof state);memcpy(before,state,sizeof state);
 lumo_solver_default_options(&options);options.max_iterations=0;
 CTEST_STATUS(lumo_solver_solve(&op,state,&options,workspace,LUMO_SOLVER_WORKSPACE_DOUBLES,&report),LUMO_SOLVER_CONVERGED);CTEST_STATUS(report.iterations,0);CTEST_STATUS(memcmp(before,state,sizeof state),0);return 0;
}
static int limits(void){
 struct lumo_mesh_operating op;double state[N],before[N];struct lumo_solver_report report;struct lumo_solver_options options;
 CTEST_STATUS(setup(&op,state),0);memcpy(before,state,sizeof state);lumo_solver_default_options(&options);options.max_iterations=0;
 CTEST_STATUS(lumo_solver_solve(&op,state,&options,workspace,LUMO_SOLVER_WORKSPACE_DOUBLES,&report),LUMO_SOLVER_MAX_ITERATIONS);CTEST_STATUS(report.iterations,0);CTEST_STATUS(memcmp(before,state,sizeof state),0);
 options.max_iterations=1;CTEST_STATUS(lumo_solver_solve(&op,state,&options,workspace,LUMO_SOLVER_WORKSPACE_DOUBLES,&report),LUMO_SOLVER_MAX_ITERATIONS);CTEST_STATUS(report.iterations,1);CTEST_STATUS(isfinite(report.residual_norm),1);CTEST_STATUS(memcmp(before,state,sizeof state)!=0,1);return 0;
}
static int zero_jacobian(const double *s,const struct lumo_mesh_operating *op,const void *ctx,double *matrix){(void)s;(void)op;(void)ctx;memset(matrix,0,N*N*sizeof(double));return 0;}
static int failed_jacobian(const double *s,const struct lumo_mesh_operating *op,const void *ctx,double *matrix){(void)s;(void)op;(void)ctx;(void)matrix;return LUMO_MESH_TEMPERATURE_RANGE;}
static int scaled_jacobian(const double *s,const struct lumo_mesh_operating *op,const void *ctx,double *matrix){size_t i;int status=lumo_mesh_jacobian_numeric(s,op,NULL,matrix);if(status!=0)return status;for(i=0;i<N*N;++i)matrix[i]*=*(const double *)ctx;return 0;}
static int nan_jacobian(const double *s,const struct lumo_mesh_operating *op,const void *ctx,double *matrix){int status=zero_jacobian(s,op,ctx,matrix);matrix[0]=NAN;return status;}
static int callbacks(void){
 struct lumo_mesh_operating op;double state[N],before[N],factor;struct lumo_solver_report report;struct lumo_solver_options options;
 CTEST_STATUS(setup(&op,state),0);memcpy(before,state,sizeof state);lumo_solver_default_options(&options);
 options.jacobian=zero_jacobian;CTEST_STATUS(lumo_solver_solve(&op,state,&options,workspace,LUMO_SOLVER_WORKSPACE_DOUBLES,&report),LUMO_SOLVER_SINGULAR_JACOBIAN);CTEST_STATUS(memcmp(before,state,sizeof state),0);
 options.jacobian=failed_jacobian;CTEST_STATUS(lumo_solver_solve(&op,state,&options,workspace,LUMO_SOLVER_WORKSPACE_DOUBLES,&report),LUMO_SOLVER_MODEL_ERROR);CTEST_STATUS(report.model_status,LUMO_MESH_TEMPERATURE_RANGE);CTEST_STATUS(memcmp(before,state,sizeof state),0);
 options.jacobian=nan_jacobian;CTEST_STATUS(lumo_solver_solve(&op,state,&options,workspace,LUMO_SOLVER_WORKSPACE_DOUBLES,&report),LUMO_SOLVER_NUMERIC_ERROR);
 options.jacobian=scaled_jacobian;factor=-1;options.jacobian_context=&factor;
 CTEST_STATUS(lumo_solver_solve(&op,state,&options,workspace,LUMO_SOLVER_WORKSPACE_DOUBLES,&report),LUMO_SOLVER_LINE_SEARCH_FAILED);CTEST_STATUS(memcmp(before,state,sizeof state),0);
 factor=0.1;options.max_backtracks=0;CTEST_STATUS(lumo_solver_solve(&op,state,&options,workspace,LUMO_SOLVER_WORKSPACE_DOUBLES,&report),LUMO_SOLVER_LINE_SEARCH_FAILED);CTEST_STATUS(memcmp(before,state,sizeof state),0);
 options.max_backtracks=30;options.max_iterations=200;
 CTEST_STATUS(lumo_solver_solve(&op,state,&options,workspace,LUMO_SOLVER_WORKSPACE_DOUBLES,&report),LUMO_SOLVER_CONVERGED);CTEST_STATUS(report.backtracks>0,1);return compare_reference(state);
}
static int errors(void){
 struct lumo_mesh_operating op;double state[N],before[N],scales[N];struct lumo_solver_report report;struct lumo_solver_options options;int t;size_t i;
 CTEST_STATUS(setup(&op,state),0);memcpy(before,state,sizeof state);
 CTEST_STATUS(lumo_solver_solve(&op,state,NULL,workspace,LUMO_SOLVER_WORKSPACE_DOUBLES,NULL),LUMO_SOLVER_INVALID_INPUT);
 CTEST_STATUS(lumo_solver_solve(NULL,state,NULL,workspace,LUMO_SOLVER_WORKSPACE_DOUBLES,&report),LUMO_SOLVER_INVALID_INPUT);
 CTEST_STATUS(lumo_solver_solve(&op,NULL,NULL,workspace,LUMO_SOLVER_WORKSPACE_DOUBLES,&report),LUMO_SOLVER_INVALID_INPUT);
 for(i=0;i<N;++i)scales[i]=1;
 for(t=0;t<8;++t){lumo_solver_default_options(&options);switch(t){case 0:options.residual_tolerance=0;break;case 1:options.residual_tolerance=NAN;break;case 2:options.max_iterations=-1;break;case 3:options.max_backtracks=-1;break;case 4:options.max_backtracks=1001;break;case 5:options.jacobian_context=&op;break;case 6:options.variable_scales=scales;scales[0]=0;break;case 7:options.residual_scales=scales;scales[0]=NAN;break;}
 CTEST_STATUS(lumo_solver_solve(&op,state,&options,workspace,LUMO_SOLVER_WORKSPACE_DOUBLES,&report),LUMO_SOLVER_INVALID_INPUT);CTEST_STATUS(memcmp(before,state,sizeof state),0);}
 state[LUMO_MESH_X]=-0.1;memcpy(before,state,sizeof state);CTEST_STATUS(lumo_solver_solve(&op,state,NULL,workspace,LUMO_SOLVER_WORKSPACE_DOUBLES,&report),LUMO_SOLVER_INVALID_INPUT);CTEST_STATUS(memcmp(before,state,sizeof state),0);
 CTEST_STATUS(setup(&op,state),0);
 state[LUMO_MESH_T]=340;memcpy(before,state,sizeof state);CTEST_STATUS(lumo_solver_solve(&op,state,NULL,workspace,LUMO_SOLVER_WORKSPACE_DOUBLES,&report),LUMO_SOLVER_MODEL_ERROR);CTEST_STATUS(report.model_status,LUMO_MESH_TEMPERATURE_RANGE);CTEST_STATUS(memcmp(before,state,sizeof state),0);
 return 0;
}
static int scales(void){
 struct lumo_mesh_operating op;double state[N],xs[N],rs[N];struct lumo_solver_report report;struct lumo_solver_options options;size_t i;
 CTEST_STATUS(setup(&op,state),0);lumo_solver_default_options(&options);
 for(i=0;i<N;++i){xs[i]=scale(i);rs[i]=i==LUMO_MESH_CLOSURE||i%LUMO_MESH_BLOCK<LUMO_MESH_COMPONENTS?100:(i%LUMO_MESH_BLOCK==LUMO_MESH_H?1e7:1);}
 options.variable_scales=xs;options.residual_scales=rs;
 CTEST_STATUS(lumo_solver_solve(&op,state,&options,workspace,LUMO_SOLVER_WORKSPACE_DOUBLES,&report),LUMO_SOLVER_CONVERGED);return compare_reference(state);
}
static int operating_variants(void){
 struct lumo_mesh_operating op;double state[N];struct lumo_solver_report report;int trial;
 for(trial=0;trial<3;++trial){CTEST_STATUS(setup(&op,state),0);op.distillate_mol_s=trial==0?40:(trial==1?60:50);op.reflux_ratio=trial==2?4:2.5;
 CTEST_STATUS(lumo_initialguess_linear_T_x(&op,state),0);CTEST_STATUS(lumo_solver_solve(&op,state,NULL,workspace,LUMO_SOLVER_WORKSPACE_DOUBLES,&report),LUMO_SOLVER_CONVERGED);CTEST_STATUS(report.residual_norm<=1e-10,1);}
 return 0;
}
static int workspace_contract(void){
 struct lumo_mesh_operating op;
 struct lumo_solver_report report;
 double state[N],before[N],saved[N];
 struct {double before;double values[LUMO_SOLVER_WORKSPACE_DOUBLES+3];double after;} guarded;
 double untouched[LUMO_SOLVER_WORKSPACE_DOUBLES+3];
 size_t i;int iterations;
 CTEST_STATUS(lumo_solver_workspace_doubles()==LUMO_SOLVER_WORKSPACE_DOUBLES,1);
 guarded.before=12345;guarded.after=67890;
 for(i=0;i<LUMO_SOLVER_WORKSPACE_DOUBLES+3;++i)guarded.values[i]=NAN;
 memcpy(untouched,guarded.values,sizeof untouched);
 CTEST_STATUS(setup(&op,state),0);memcpy(before,state,sizeof state);
 CTEST_STATUS(lumo_solver_solve(&op,state,NULL,NULL,LUMO_SOLVER_WORKSPACE_DOUBLES,&report),LUMO_SOLVER_INVALID_WORKSPACE);
 CTEST_STATUS(report.status,LUMO_SOLVER_INVALID_WORKSPACE);CTEST_STATUS(report.iterations,0);
 CTEST_STATUS(lumo_solver_solve(&op,state,NULL,guarded.values,LUMO_SOLVER_WORKSPACE_DOUBLES-1,&report),LUMO_SOLVER_INVALID_WORKSPACE);
 CTEST_STATUS(memcmp(state,before,sizeof state),0);CTEST_STATUS(memcmp(guarded.values,untouched,sizeof untouched),0);
 CTEST_STATUS(lumo_solver_solve(&op,state,NULL,guarded.values,LUMO_SOLVER_WORKSPACE_DOUBLES+3,&report),LUMO_SOLVER_CONVERGED);
 CTEST_STATUS(compare_reference(state),0);memcpy(saved,state,sizeof state);iterations=report.iterations;
 CTEST_NEAR(guarded.before,12345);CTEST_NEAR(guarded.after,67890);
 for(i=LUMO_SOLVER_WORKSPACE_DOUBLES;i<LUMO_SOLVER_WORKSPACE_DOUBLES+3;++i)CTEST_STATUS(isnan(guarded.values[i]),1);
 /* Reuse previously overwritten memory without clearing or initializing it. */
 CTEST_STATUS(setup(&op,state),0);
 CTEST_STATUS(lumo_solver_solve(&op,state,NULL,guarded.values,LUMO_SOLVER_WORKSPACE_DOUBLES,&report),LUMO_SOLVER_CONVERGED);
 CTEST_STATUS(report.iterations,iterations);CTEST_STATUS(memcmp(state,saved,sizeof state),0);
 CTEST_NEAR(guarded.before,12345);CTEST_NEAR(guarded.after,67890);
 return 0;
}
int main(int argc,char **argv){const struct ctest_case cases[]={{"workspace_contract",workspace_contract},{"reference",reference},{"already_converged",already_converged},{"limits",limits},{"callbacks",callbacks},{"errors",errors},{"scales",scales},{"operating_variants",operating_variants}};return ctest_run(argc,argv,cases,sizeof cases/sizeof cases[0]);}
