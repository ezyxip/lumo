#include "lumo_thermo.h"
#include "ctest.h"
#include "case.h"

static int numerical_boundary(void) {
    struct lumo_flash_result result, before;
    const double unit[] = {1,1}, zero[] = {0,0};
    double e[2], hl, hv;
    long double k0, k1, a, b, expected;
    int status;
    size_t i;
    CTEST_STATUS(lumo_thermo_eval(350,101325,unit,zero,e,&hl,&hv),LUMO_THERMO_OK);
    k0=-e[0]; k1=-e[1]; a=k0-1; b=k1-1;
    /* Independent closed-form binary solution, using the actual rounded K. */
    expected=-((long double)z[0]*a+(long double)z[1]*b)/(a*b*((long double)z[0]+z[1]));
    memset(&result,0x34,sizeof result);memcpy(&before,&result,sizeof result);
    status=lumo_thermo_flash(LUMO_FEED_TEMPERATURE,350,101325,z,&result);
#if TRACE_LIQUID
    /* A beta too close to 1 may be unrepresentable with accurate phase sums.
     * An explicit failure is preferable to a successful unnormalized result.
     */
    if(status==LUMO_THERMO_FLASH_NO_CONVERGENCE){
        CTEST_STATUS(memcmp(&result,&before,sizeof result),0);
        return 0;
    }
#endif
    CTEST_STATUS(status,LUMO_THERMO_OK);
    CTEST_STATUS(fabs(result.vapor_fraction-(double)expected)<2e-9,1);
    CTEST_STATUS(fabs(result.x[0]+result.x[1]-1)<1e-12,1);
    CTEST_STATUS(fabs(result.y[0]+result.y[1]-1)<1e-12,1);
    for(i=0;i<2;++i){
        CTEST_STATUS(fabs((1-result.vapor_fraction)*result.x[i]+result.vapor_fraction*result.y[i]-z[i])<1e-12,1);
    }
    CTEST_STATUS(lumo_thermo_eval(350,101325,result.x,result.y,e,&hl,&hv),LUMO_THERMO_OK);
    for(i=0;i<2;++i)CTEST_NEAR(e[i],0);
    CTEST_NEAR(result.enthalpy_J_mol,(1-result.vapor_fraction)*hl+result.vapor_fraction*hv);
    return 0;
}
int main(int argc,char **argv){
    const struct ctest_case cases[]={{"numerical_boundary",numerical_boundary}};
    return ctest_run(argc,argv,cases,1);
}
