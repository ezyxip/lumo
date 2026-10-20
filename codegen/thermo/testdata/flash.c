#include "lumo_thermo.h"
#include "ctest.h"
#define N LUMO_THERMO_COMPONENTS
#define P 101325.0

static void composition(double *z) {
    z[0] = 0.5;
#if LUMO_THERMO_COMPONENTS == 2
    z[1] = 0.5;
#else
    z[1] = 0.3; z[2] = 0.2;
#endif
}

static int reference(void) {
    struct lumo_flash_result r;
    double z[N];
    composition(z);
    CTEST_STATUS(lumo_thermo_flash(LUMO_FEED_TEMPERATURE,350.0,P,z,&r),LUMO_THERMO_OK);
#if LUMO_THERMO_COMPONENTS == 2
    /* Binary Rachford-Rice has the exact solution beta=1/2 at K=(2,1/2). */
    CTEST_NEAR(r.vapor_fraction,0.5);
    CTEST_NEAR(r.x[0],1.0/3.0); CTEST_NEAR(r.x[1],2.0/3.0);
    CTEST_NEAR(r.y[0],2.0/3.0); CTEST_NEAR(r.y[1],1.0/3.0);
    CTEST_NEAR(r.enthalpy_J_mol,10000.0);
#else
    /* Published example in chemicals.rachford_rice, independent of our code.
     * https://chemicals.readthedocs.io/chemicals.rachford_rice.html
     */
    const double x[] = {0.33940869696634357,0.3650560590371706,0.2955352439964858};
    const double y[] = {0.5719036543882889,0.27087159580558057,0.15722474980613044};
    size_t i;
    CTEST_NEAR(r.vapor_fraction,0.6907302627738544);
    for(i=0;i<N;++i){CTEST_NEAR(r.x[i],x[i]);CTEST_NEAR(r.y[i],y[i]);}
    CTEST_NEAR(r.enthalpy_J_mol,20000*0.6907302627738544);
#endif
    return 0;
}

static int state_roundtrip(void) {
    struct lumo_flash_result r, tp, saturated;
    double z[N], e[N], hl, hv;
    const double fractions[] = {0,0.2,0.5,0.9,1};
    size_t i,j;
    composition(z);
    for(j=0;j<sizeof fractions/sizeof fractions[0];++j){
        double sx=0,sy=0;
        CTEST_STATUS(lumo_thermo_flash(LUMO_FEED_VAPOR_FRACTION,fractions[j],P,z,&r),LUMO_THERMO_OK);
        CTEST_NEAR(r.vapor_fraction,fractions[j]);
        CTEST_STATUS(lumo_thermo_eval(r.temperature_K,P,r.x,r.y,e,&hl,&hv),LUMO_THERMO_OK);
        for(i=0;i<N;++i){sx+=r.x[i];sy+=r.y[i];CTEST_NEAR(e[i],0);CTEST_NEAR((1-fractions[j])*r.x[i]+fractions[j]*r.y[i],z[i]);}
        CTEST_NEAR(sx,1);CTEST_NEAR(sy,1);
        CTEST_NEAR(r.enthalpy_J_mol,(1-fractions[j])*hl+fractions[j]*hv);
        CTEST_STATUS(lumo_thermo_flash(LUMO_FEED_TEMPERATURE,r.temperature_K,P,z,&tp),LUMO_THERMO_OK);
        CTEST_STATUS(fabs(tp.vapor_fraction-fractions[j])<1e-10,1);
        CTEST_NEAR(tp.enthalpy_J_mol,r.enthalpy_J_mol);
        if(j==0){
            CTEST_STATUS(lumo_thermo_flash(LUMO_FEED_SATURATED_LIQUID,NAN,P,z,&saturated),LUMO_THERMO_OK);
            CTEST_NEAR(saturated.temperature_K,r.temperature_K);CTEST_NEAR(saturated.enthalpy_J_mol,r.enthalpy_J_mol);
        }
    }
    return 0;
}

static int single_phase(void) {
    struct lumo_flash_result r;
    double z[N];size_t i;
    composition(z);
    CTEST_STATUS(lumo_thermo_flash(LUMO_FEED_TEMPERATURE,330,P,z,&r),LUMO_THERMO_OK);
    CTEST_NEAR(r.vapor_fraction,0);CTEST_NEAR(r.enthalpy_J_mol,-200);
    for(i=0;i<N;++i)CTEST_NEAR(r.x[i],z[i]);
    CTEST_STATUS(lumo_thermo_flash(LUMO_FEED_TEMPERATURE,400,P,z,&r),LUMO_THERMO_OK);
    CTEST_NEAR(r.vapor_fraction,1);CTEST_NEAR(r.enthalpy_J_mol,21000);
    for(i=0;i<N;++i)CTEST_NEAR(r.y[i],z[i]);
    return 0;
}

static int pure_component(void) {
    struct lumo_flash_result r,before;
    double z[N]={0};z[1]=1;
    CTEST_STATUS(lumo_thermo_flash(LUMO_FEED_VAPOR_FRACTION,0.3,P,z,&r),LUMO_THERMO_OK);
    CTEST_NEAR(r.x[1],1);CTEST_NEAR(r.y[1],1);CTEST_NEAR(r.vapor_fraction,0.3);
    before=r;
    CTEST_STATUS(lumo_thermo_flash(LUMO_FEED_TEMPERATURE,r.temperature_K,P,z,&r),LUMO_THERMO_FLASH_AMBIGUOUS);
    CTEST_STATUS(memcmp(&before,&r,sizeof r),0);
    return 0;
}

static int errors(void) {
    struct lumo_flash_result r,before;
    double z[N];size_t trial;
    const double invalid[]={-1,1.1,NAN,INFINITY};
    composition(z);memset(&r,0x45,sizeof r);memcpy(&before,&r,sizeof r);
    CTEST_STATUS(lumo_thermo_flash(LUMO_FEED_TEMPERATURE,350,P,NULL,&r),LUMO_THERMO_INVALID_INPUT);
    CTEST_STATUS(lumo_thermo_flash(LUMO_FEED_TEMPERATURE,350,P,z,NULL),LUMO_THERMO_INVALID_INPUT);
    CTEST_STATUS(lumo_thermo_flash((enum lumo_feed_state)99,0,P,z,&r),LUMO_THERMO_INVALID_INPUT);
    for(trial=0;trial<sizeof invalid/sizeof invalid[0];++trial){
        CTEST_STATUS(lumo_thermo_flash(LUMO_FEED_VAPOR_FRACTION,invalid[trial],P,z,&r),LUMO_THERMO_INVALID_INPUT);
    }
    CTEST_STATUS(lumo_thermo_flash(LUMO_FEED_TEMPERATURE,NAN,P,z,&r),LUMO_THERMO_INVALID_INPUT);
    CTEST_STATUS(lumo_thermo_flash(LUMO_FEED_TEMPERATURE,320,P,z,&r),LUMO_THERMO_TEMPERATURE_RANGE);
    CTEST_STATUS(lumo_thermo_flash(LUMO_FEED_SATURATED_LIQUID,0,1e9,z,&r),LUMO_THERMO_TEMPERATURE_RANGE);
    CTEST_STATUS(lumo_thermo_flash(LUMO_FEED_SATURATED_LIQUID,0,1,z,&r),LUMO_THERMO_TEMPERATURE_RANGE);
    CTEST_STATUS(lumo_thermo_flash(LUMO_FEED_TEMPERATURE,350,0,z,&r),LUMO_THERMO_INVALID_INPUT);
    z[0]=0;
    CTEST_STATUS(lumo_thermo_flash(LUMO_FEED_TEMPERATURE,350,P,z,&r),LUMO_THERMO_INVALID_INPUT);
    z[0]=NAN;
    CTEST_STATUS(lumo_thermo_flash(LUMO_FEED_TEMPERATURE,350,P,z,&r),LUMO_THERMO_INVALID_INPUT);
    CTEST_STATUS(memcmp(&r,&before,sizeof r),0);
    return 0;
}

int main(int argc,char **argv){
    const struct ctest_case cases[]={{"reference",reference},{"state_roundtrip",state_roundtrip},{"single_phase",single_phase},{"errors",errors},{"pure_component",pure_component}};
    return ctest_run(argc,argv,cases,sizeof cases/sizeof cases[0]);
}
