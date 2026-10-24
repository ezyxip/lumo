#include "lumo_solver.h"
#include "lumo_prepare.h"
#include "lumo_initialguess.h"
#include "ctest.h"
#define N LUMO_MESH_SIZE
int off_solve(const struct lumo_mesh_operating *,double *,const struct lumo_solver_options *,double *,size_t,struct lumo_solver_report *);
static void count_event(const struct lumo_solver_iteration *event,void *context){size_t *count=context;(void)event;++*count;}
static int equivalence(void){
    const struct lumo_prepare_input input={.feeds={{.flow=100,.fractions={0.5,0.5}}},.distillate_flow=50,.reflux_ratio=2.5};
    struct lumo_mesh_operating op;struct lumo_solver_options options;struct lumo_solver_report report,baseline;
    double initial[N],state[N],expected[N],workspace[LUMO_SOLVER_WORKSPACE_DOUBLES];
    size_t count=0;int mode;
    CTEST_STATUS(lumo_prepare_operating(&input,&op),0);CTEST_STATUS(lumo_initialguess_linear_T_x(&op,initial),0);
    memcpy(expected,initial,sizeof initial);
    CTEST_STATUS(off_solve(&op,expected,NULL,workspace,LUMO_SOLVER_WORKSPACE_DOUBLES,&baseline),0);
    for(mode=0;mode<3;++mode){
        memcpy(state,initial,sizeof initial);lumo_solver_default_options(&options);count=0;
        if(mode>0){options.observer=count_event;options.observer_context=&count;}
        if(mode==2)CTEST_STATUS(off_solve(&op,state,&options,workspace,LUMO_SOLVER_WORKSPACE_DOUBLES,&report),0);
        else CTEST_STATUS(lumo_solver_solve(&op,state,&options,workspace,LUMO_SOLVER_WORKSPACE_DOUBLES,&report),0);
        CTEST_STATUS(memcmp(expected,state,sizeof state),0);
        CTEST_STATUS(report.status,baseline.status);CTEST_STATUS(report.model_status,baseline.model_status);
        CTEST_STATUS(report.iterations,baseline.iterations);CTEST_STATUS(report.backtracks,baseline.backtracks);
        CTEST_STATUS(report.residual_norm==baseline.residual_norm,1);CTEST_STATUS(report.step_norm==baseline.step_norm,1);
        CTEST_STATUS(count,mode==1?(size_t)report.iterations+1:0);
    }
    return 0;
}
int main(int argc,char **argv){const struct ctest_case cases[]={{"equivalence",equivalence}};return ctest_run(argc,argv,cases,1);}
