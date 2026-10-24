#include "lumo_solver.h"
#include "lumo_prepare.h"
#include "lumo_initialguess.h"
#include "ctest.h"
#define N LUMO_MESH_SIZE
static double workspace[LUMO_SOLVER_WORKSPACE_DOUBLES];
struct capture {
    struct lumo_solver_iteration records[101];
    size_t count;
    int bad;
    const double *state;
    const struct lumo_mesh_operating *operating;
    double xs[N], rs[N], previous[N];
};
static void observe(const struct lumo_solver_iteration *record,void *context){
    struct capture *c=context;
    double r[N],worst_abs=fabs(record->worst_residual_value),step=0;
    size_t i;
    if(c->count>=101){c->bad=1;return;}
    if(record->iteration!=(int)c->count || record->worst_residual_index>=N){c->bad=1;return;}
    if(lumo_mesh_residual(c->state,c->operating,r)!=0){c->bad=1;return;}
    for(i=0;i<N;++i){
        double value=r[i]/c->rs[i];
        if(fabs(value)>worst_abs)c->bad=1;
        if(i<record->worst_residual_index && fabs(value)==worst_abs)c->bad=1;
        if(i==record->worst_residual_index && value!=record->worst_residual_value)c->bad=1;
        step=fmax(step,fabs((c->state[i]-c->previous[i])/c->xs[i]));
    }
    if(worst_abs!=record->residual_norm)c->bad=1;
    if(c->count==0){
        if(record->event!=LUMO_SOLVER_INITIAL || record->step_norm!=0 || record->step_length!=0 || record->rejected_trials!=0 || record->total_rejected_trials!=0)c->bad=1;
    }else{
        const struct lumo_solver_iteration *previous=&c->records[c->count-1];
        if(record->event!=LUMO_SOLVER_ITERATION || record->residual_norm>=previous->residual_norm || record->step_length<=0 || record->step_length>1 || record->step_norm!=step)c->bad=1;
        if(record->total_rejected_trials!=previous->total_rejected_trials+record->rejected_trials)c->bad=1;
    }
    c->records[c->count++]=*record;
    memcpy(c->previous,c->state,sizeof c->previous);
}
static int setup(struct lumo_mesh_operating *op,double *state,struct lumo_solver_options *options,struct capture *c){
    const struct lumo_prepare_input input={.feeds={{.flow=100,.fractions={0.5,0.5}}},.distillate_flow=50,.reflux_ratio=2.5};
    size_t i;
    CTEST_STATUS(lumo_prepare_operating(&input,op),0);
    CTEST_STATUS(lumo_initialguess_linear_T_x(op,state),0);
    memset(c,0,sizeof *c);c->state=state;c->operating=op;
    lumo_solver_default_options(options);options->observer=observe;options->observer_context=c;
    for(i=0;i<N;++i){size_t offset=i%LUMO_MESH_BLOCK;c->xs[i]=1;c->rs[i]=1;
        if(i==LUMO_MESH_QD||i==LUMO_MESH_QW)c->xs[i]=1e7;
        else if(offset==LUMO_MESH_T)c->xs[i]=350;
        else if(offset==LUMO_MESH_L||offset==LUMO_MESH_V)c->xs[i]=100;
        if(i==LUMO_MESH_CLOSURE||offset<LUMO_MESH_COMPONENTS)c->rs[i]=100;
        else if(offset==LUMO_MESH_H)c->rs[i]=1e7;
    }
    options->variable_scales=c->xs;options->residual_scales=c->rs;
    return 0;
}
static int observation(void){
    struct lumo_mesh_operating op;struct lumo_solver_report report;struct lumo_solver_options options;
    struct capture c;double state[N],saved[N];
    CTEST_STATUS(setup(&op,state,&options,&c),0);
    CTEST_STATUS(lumo_solver_solve(&op,state,&options,workspace,LUMO_SOLVER_WORKSPACE_DOUBLES,&report),0);
    CTEST_STATUS(c.bad,0);
#if LUMO_SOLVER_HAS_DIAGNOSTICS
    CTEST_STATUS(c.count,(size_t)report.iterations+1);
    CTEST_NEAR(c.records[c.count-1].residual_norm,report.residual_norm);
    CTEST_NEAR(c.records[c.count-1].step_norm,report.step_norm);
    CTEST_STATUS(c.records[c.count-1].total_rejected_trials,report.backtracks);
#else
    CTEST_STATUS(c.count,0);
#endif
    memcpy(saved,state,sizeof saved);
    /* Already converged still has exactly one initial event. Caller resets storage. */
    c.count=0;c.bad=0;
    CTEST_STATUS(lumo_solver_solve(&op,state,&options,workspace,LUMO_SOLVER_WORKSPACE_DOUBLES,&report),0);
    CTEST_STATUS(report.iterations,0);CTEST_STATUS(c.count,LUMO_SOLVER_HAS_DIAGNOSTICS);
    CTEST_STATUS(memcmp(saved,state,sizeof saved),0);CTEST_STATUS(c.bad,0);
    /* NULL observer ignores even a non-NULL context. */
    options.observer=NULL;options.observer_context=(void *)&report;
    CTEST_STATUS(lumo_solver_solve(&op,state,&options,workspace,LUMO_SOLVER_WORKSPACE_DOUBLES,&report),0);
    CTEST_STATUS(c.count,LUMO_SOLVER_HAS_DIAGNOSTICS);
    return 0;
}
static int zero_jacobian(const double *s,const struct lumo_mesh_operating *op,const void *ctx,double *matrix){(void)s;(void)op;(void)ctx;memset(matrix,0,N*N*sizeof(double));return 0;}
static int enlarged_step(const double *s,const struct lumo_mesh_operating *op,const void *ctx,double *matrix){size_t i;int status;(void)ctx;status=lumo_mesh_jacobian_numeric(s,op,NULL,matrix);if(status)return status;for(i=0;i<N*N;++i)matrix[i]*=0.1;return 0;}
static int observation_failures(void){
    struct lumo_mesh_operating op;struct lumo_solver_report report;struct lumo_solver_options options;
    struct capture c;double state[N];
    CTEST_STATUS(setup(&op,state,&options,&c),0);options.max_iterations=0;
    CTEST_STATUS(lumo_solver_solve(&op,state,&options,workspace,LUMO_SOLVER_WORKSPACE_DOUBLES,&report),LUMO_SOLVER_MAX_ITERATIONS);
    CTEST_STATUS(c.count,LUMO_SOLVER_HAS_DIAGNOSTICS);CTEST_STATUS(c.bad,0);
    CTEST_STATUS(setup(&op,state,&options,&c),0);options.max_iterations=1;
    CTEST_STATUS(lumo_solver_solve(&op,state,&options,workspace,LUMO_SOLVER_WORKSPACE_DOUBLES,&report),LUMO_SOLVER_MAX_ITERATIONS);
    CTEST_STATUS(c.count,2*LUMO_SOLVER_HAS_DIAGNOSTICS);CTEST_STATUS(c.bad,0);
    CTEST_STATUS(setup(&op,state,&options,&c),0);options.jacobian=zero_jacobian;
    CTEST_STATUS(lumo_solver_solve(&op,state,&options,workspace,LUMO_SOLVER_WORKSPACE_DOUBLES,&report),LUMO_SOLVER_SINGULAR_JACOBIAN);
    CTEST_STATUS(c.count,LUMO_SOLVER_HAS_DIAGNOSTICS);CTEST_STATUS(c.bad,0);
    CTEST_STATUS(setup(&op,state,&options,&c),0);state[LUMO_MESH_T]=340;
    CTEST_STATUS(lumo_solver_solve(&op,state,&options,workspace,LUMO_SOLVER_WORKSPACE_DOUBLES,&report),LUMO_SOLVER_MODEL_ERROR);CTEST_STATUS(c.count,0);
    CTEST_STATUS(setup(&op,state,&options,&c),0);options.residual_tolerance=-1;
    CTEST_STATUS(lumo_solver_solve(&op,state,&options,workspace,LUMO_SOLVER_WORKSPACE_DOUBLES,&report),LUMO_SOLVER_INVALID_INPUT);CTEST_STATUS(c.count,0);
    CTEST_STATUS(setup(&op,state,&options,&c),0);
    CTEST_STATUS(lumo_solver_solve(&op,state,&options,NULL,0,&report),LUMO_SOLVER_INVALID_WORKSPACE);CTEST_STATUS(c.count,0);
    return 0;
}
static int backtracking_events(void){
    struct lumo_mesh_operating op;struct lumo_solver_report report;struct lumo_solver_options options;
    struct capture c;double state[N];
    CTEST_STATUS(setup(&op,state,&options,&c),0);options.jacobian=enlarged_step;
    CTEST_STATUS(lumo_solver_solve(&op,state,&options,workspace,LUMO_SOLVER_WORKSPACE_DOUBLES,&report),0);
    CTEST_STATUS(report.backtracks>0,1);CTEST_STATUS(c.bad,0);
    CTEST_STATUS(c.count,LUMO_SOLVER_HAS_DIAGNOSTICS*((size_t)report.iterations+1));
    CTEST_STATUS(setup(&op,state,&options,&c),0);options.jacobian=enlarged_step;options.max_backtracks=0;
    CTEST_STATUS(lumo_solver_solve(&op,state,&options,workspace,LUMO_SOLVER_WORKSPACE_DOUBLES,&report),LUMO_SOLVER_LINE_SEARCH_FAILED);
    CTEST_STATUS(c.count,LUMO_SOLVER_HAS_DIAGNOSTICS);CTEST_STATUS(c.bad,0);
    return 0;
}
int main(int argc,char **argv){const struct ctest_case cases[]={{"observation",observation},{"observation_failures",observation_failures},{"backtracking_events",backtracking_events}};return ctest_run(argc,argv,cases,sizeof cases/sizeof cases[0]);}
