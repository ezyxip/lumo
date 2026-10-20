#include "lumo_prepare.h"
#include "ctest.h"
#include <float.h>

static struct lumo_prepare_input input_default(void) {
    struct lumo_prepare_input input = {0};
    size_t f;
    for(f=0;f<LUMO_MESH_FEEDS;++f){
        input.feeds[f].flow=100;
        input.feeds[f].flow_unit=LUMO_FLOW_MOL_S;
        input.feeds[f].composition_basis=LUMO_COMPOSITION_MOLE;
        input.feeds[f].fractions[0]=0.5;input.feeds[f].fractions[1]=0.5;
        input.feeds[f].state=LUMO_FEED_SATURATED_LIQUID;
    }
    input.distillate_flow=50;input.distillate_unit=LUMO_FLOW_MOL_S;input.reflux_ratio=2.5;
    return input;
}

static int same_operating(const struct lumo_mesh_operating *a,const struct lumo_mesh_operating *b){
    size_t f,i;
    CTEST_NEAR(a->distillate_mol_s,b->distillate_mol_s);CTEST_NEAR(a->reflux_ratio,b->reflux_ratio);
    for(f=0;f<LUMO_MESH_FEEDS;++f){
        CTEST_NEAR(a->feeds[f].flow_mol_s,b->feeds[f].flow_mol_s);
        CTEST_NEAR(a->feeds[f].enthalpy_J_mol,b->feeds[f].enthalpy_J_mol);
        for(i=0;i<LUMO_MESH_COMPONENTS;++i)CTEST_NEAR(a->feeds[f].mole_fractions[i],b->feeds[f].mole_fractions[i]);
    }
    return 0;
}

static int reference(void){
    struct lumo_prepare_input input=input_default();struct lumo_mesh_operating op;
    size_t f;
    CTEST_STATUS(lumo_prepare_operating(&input,&op),LUMO_PREPARE_OK);
    for(f=0;f<LUMO_MESH_FEEDS;++f){
        CTEST_NEAR(op.feeds[f].flow_mol_s,100);
        CTEST_NEAR(op.feeds[f].enthalpy_J_mol,2060.414799029395);
        CTEST_NEAR(op.feeds[f].mole_fractions[0],0.5);
    }
    CTEST_NEAR(op.distillate_mol_s,50);CTEST_NEAR(op.reflux_ratio,2.5);
    return 0;
}

static int conversions(void){
    struct lumo_prepare_input input=input_default();struct lumo_mesh_operating expected,actual;
    size_t f;int basis,unit;
    CTEST_STATUS(lumo_prepare_operating(&input,&expected),LUMO_PREPARE_OK);
    for(basis=0;basis<2;++basis)for(unit=0;unit<3;++unit){
        input=input_default();input.distillate_flow=180;input.distillate_unit=LUMO_FLOW_KMOL_H;
        for(f=0;f<LUMO_MESH_FEEDS;++f){
            struct lumo_prepare_feed *feed=&input.feeds[f];
            feed->flow_unit=(enum lumo_flow_unit)unit;
            feed->flow=unit==0?100:(unit==1?360:100*(78.11+92.14)*0.5*3.6);
            if(basis){
                feed->composition_basis=LUMO_COMPOSITION_MASS;
                feed->fractions[0]=78.11/(78.11+92.14);feed->fractions[1]=92.14/(78.11+92.14);
            }
        }
        CTEST_STATUS(lumo_prepare_operating(&input,&actual),LUMO_PREPARE_OK);
        CTEST_STATUS(same_operating(&actual,&expected),0);
    }
    /* Near-normalized user fractions are normalized consistently. */
    input=input_default();input.feeds[0].fractions[0]+=5e-11;
    CTEST_STATUS(lumo_prepare_operating(&input,&actual),LUMO_PREPARE_OK);
    CTEST_NEAR(actual.feeds[0].mole_fractions[0]+actual.feeds[0].mole_fractions[1],1);
    return 0;
}

static int feed_states(void){
    struct lumo_prepare_input input=input_default();struct lumo_mesh_operating op;
    double betas[]={0,0.4,1};size_t j;
    input.feeds[0].state=LUMO_FEED_TEMPERATURE;input.feeds[0].state_value=360;
    CTEST_STATUS(lumo_prepare_operating(&input,&op),LUMO_PREPARE_OK);
    CTEST_NEAR(op.feeds[0].enthalpy_J_mol,1350);
    input.feeds[0].state_value=390;
    CTEST_STATUS(lumo_prepare_operating(&input,&op),LUMO_PREPARE_OK);
    CTEST_NEAR(op.feeds[0].enthalpy_J_mol,35500);
    for(j=0;j<3;++j){
        struct lumo_flash_result flash;
        const double z[LUMO_MESH_COMPONENTS]={0.5,0.5};
        input.feeds[0].state=LUMO_FEED_VAPOR_FRACTION;input.feeds[0].state_value=betas[j];
        CTEST_STATUS(lumo_prepare_operating(&input,&op),LUMO_PREPARE_OK);
        CTEST_STATUS(lumo_thermo_flash(LUMO_FEED_VAPOR_FRACTION,betas[j],LUMO_MESH_PRESSURE_PA,z,&flash),LUMO_THERMO_OK);
        CTEST_NEAR(op.feeds[0].enthalpy_J_mol,flash.enthalpy_J_mol);
    }
    input=input_default();input.feeds[0].fractions[0]=1;input.feeds[0].fractions[1]=0;
    CTEST_STATUS(lumo_prepare_operating(&input,&op),LUMO_PREPARE_OK);
    CTEST_NEAR(op.feeds[0].mole_fractions[0],1);
    return 0;
}

static int errors(void){
    struct lumo_prepare_input input;struct lumo_mesh_operating op,before;
    int trial;
    char name[64];
    memset(&op,0x34,sizeof op);memcpy(&before,&op,sizeof op);
    CTEST_STATUS(lumo_prepare_operating(NULL,&op),LUMO_PREPARE_INVALID_INPUT);
    input=input_default();CTEST_STATUS(lumo_prepare_operating(&input,NULL),LUMO_PREPARE_INVALID_INPUT);
    for(trial=0;trial<21;++trial){
        int expected=LUMO_PREPARE_INVALID_INPUT;input=input_default();
        (void)snprintf(name,sizeof name,"errors trial=%d",trial);ctest_set_name(name);
        switch(trial){
        case 0:input.distillate_unit=LUMO_FLOW_KG_H;expected=LUMO_PREPARE_UNSUPPORTED_SPECIFICATION;break;
        case 1:input.distillate_unit=(enum lumo_flow_unit)99;break;
        case 2:input.distillate_flow=100*LUMO_MESH_FEEDS;break;
        case 3:input.distillate_flow=0;break;
        case 4:input.reflux_ratio=-1;break;
        case 5:input.reflux_ratio=NAN;break;
        case 6:input.feeds[0].flow=0;break;
        case 7:input.feeds[0].flow=NAN;break;
        case 8:input.feeds[0].flow_unit=(enum lumo_flow_unit)99;break;
        case 9:input.feeds[0].composition_basis=(enum lumo_composition_basis)99;break;
        case 10:input.feeds[0].fractions[0]=-0.5;break;
        case 11:input.feeds[0].fractions[0]=NAN;break;
        case 12:input.feeds[0].fractions[0]=0.6;break;
        case 13:input.feeds[0].state=(enum lumo_feed_state)99;break;
        case 14:input.feeds[0].state=LUMO_FEED_TEMPERATURE;input.feeds[0].state_value=340;expected=LUMO_PREPARE_TEMPERATURE_RANGE;break;
        case 15:input.feeds[0].state=LUMO_FEED_VAPOR_FRACTION;input.feeds[0].state_value=1.2;break;
        case 16:input.distillate_flow=INFINITY;break;
        case 17:input.feeds[0].flow=INFINITY;break;
        case 18:input.feeds[0].state=LUMO_FEED_TEMPERATURE;input.feeds[0].state_value=NAN;break;
        case 19:{size_t f;for(f=0;f<LUMO_MESH_FEEDS;++f){input.feeds[f].flow=DBL_MIN;input.feeds[f].flow_unit=LUMO_FLOW_KG_H;}input.distillate_flow=DBL_MIN;break;}
        case 20:input.feeds[0].state=LUMO_FEED_VAPOR_FRACTION;input.feeds[0].state_value=NAN;break;
        }
        CTEST_STATUS(lumo_prepare_operating(&input,&op),expected);
        CTEST_STATUS(memcmp(&op,&before,sizeof op),0);
    }
    return 0;
}

static int atomic_failure(void){
    struct lumo_prepare_input input=input_default();struct lumo_mesh_operating op,before;
    memset(&op,0x23,sizeof op);memcpy(&before,&op,sizeof op);
    input.feeds[LUMO_MESH_FEEDS-1].state=LUMO_FEED_TEMPERATURE;
    input.feeds[LUMO_MESH_FEEDS-1].state_value=400;
    CTEST_STATUS(lumo_prepare_operating(&input,&op),LUMO_PREPARE_TEMPERATURE_RANGE);
    CTEST_STATUS(memcmp(&op,&before,sizeof op),0);
    return 0;
}

static int mesh(void){
    struct lumo_prepare_input input=input_default();struct lumo_mesh_operating op;
    double state[LUMO_MESH_SIZE]={0},r[LUMO_MESH_SIZE];size_t j;
    CTEST_STATUS(lumo_prepare_operating(&input,&op),LUMO_PREPARE_OK);
    for(j=0;j<LUMO_MESH_STAGES;++j){
        size_t b=j*LUMO_MESH_BLOCK;
        state[b+LUMO_MESH_X]=0.5;state[b+LUMO_MESH_X+1]=0.5;
        state[b+LUMO_MESH_Y]=0.5;state[b+LUMO_MESH_Y+1]=0.5;
        state[b+LUMO_MESH_T]=370;state[b+LUMO_MESH_L]=100;state[b+LUMO_MESH_V]=150;
    }
    CTEST_STATUS(lumo_mesh_residual(state,&op,r),LUMO_MESH_OK);
    return 0;
}

int main(int argc,char **argv){
    const struct ctest_case cases[]={{"reference",reference},{"conversions",conversions},{"feed_states",feed_states},{"errors",errors},{"atomic_failure",atomic_failure},{"mesh",mesh}};
    return ctest_run(argc,argv,cases,sizeof cases/sizeof cases[0]);
}
