#include "lumo_thermo.h"
#include <float.h>
#include "ctest.h"

static int water_reference(void) {
    double x[2] = {0.25, 0.75}, y[2] = {0.6, 0.4};
    double r[2], hl, hv;
    const double temperatures[4] = {344.0, 350.0, 360.0, 373.0};
    /* Reference Psat in Pa from the original NIST bar/K fit, not measurements. */
    const double water_psat[4] = {
        32333.58468078457, 41641.6707822834,
        62134.389249761814, 100786.71200447337
    };
    size_t j;
    for (j = 0; j < 4; ++j) {
        char name[80];
        (void)snprintf(name, sizeof name, "water_reference T=%.2f K", temperatures[j]);
        ctest_set_name(name);
        CTEST_STATUS(lumo_thermo_eval(temperatures[j], 100000.0, x, y, r, &hl, &hv), LUMO_THERMO_OK);
        CTEST_NEAR(r[0], y[0] - x[0] * water_psat[j] / 100000.0);
        CTEST_NEAR(r[1], y[1] - x[1] * 133322.4 / 100000.0);
    }
    return 0;
}

static int enthalpy(void) {
    double x[2] = {0.25, 0.75}, y[2] = {0.6, 0.4};
    double r[2], hl, hv;
    CTEST_STATUS(lumo_thermo_eval(360.0, 100000.0, x, y, r, &hl, &hv), LUMO_THERMO_OK);
    CTEST_NEAR(hl, 1750.0);
    CTEST_NEAR(hv, 34620.0);
    return 0;
}

static int pressure(void) {
    double x[2] = {0.25, 0.75}, y[2] = {0.6, 0.4};
    double r[2], high_pressure[2], hl, hv;
    CTEST_STATUS(lumo_thermo_eval(360.0, 100000.0, x, y, r, &hl, &hv), LUMO_THERMO_OK);
    CTEST_STATUS(lumo_thermo_eval(360.0, 200000.0, x, y, high_pressure, &hl, &hv), LUMO_THERMO_OK);
    CTEST_NEAR(y[0] - high_pressure[0], 0.5 * (y[0] - r[0]));
    CTEST_NEAR(y[1] - high_pressure[1], 0.5 * (y[1] - r[1]));
    CTEST_NEAR(hl, 1750.0);
    CTEST_NEAR(hv, 34620.0);
    return 0;
}

static int reference_temperature(void) {
    double x[2] = {0.25, 0.75}, y[2] = {0.6, 0.4};
    double r[2], hl, hv;
    CTEST_STATUS(lumo_thermo_eval(350.0, 100000.0, x, y, r, &hl, &hv), LUMO_THERMO_OK);
    CTEST_NEAR(hl, 0.0);
    CTEST_NEAR(hv, 34000.0);
    return 0;
}

static int newton_iterates(void) {
    double x[2] = {-0.1, 1.2}, y[2] = {0.6, 0.4};
    double r[2], hl, hv;
    CTEST_STATUS(lumo_thermo_eval(360.0, 100000.0, x, y, r, &hl, &hv), LUMO_THERMO_OK);
    CTEST_NEAR(r[0], 0.6 + 0.1 * 62134.389249761814 / 100000.0);
    CTEST_NEAR(hl, 2300.0);
    return 0;
}

#define CHECK_UNCHANGED() do { \
    CTEST_NEAR(r[0], 42.0); CTEST_NEAR(r[1], 42.0); \
    CTEST_NEAR(hl, 42.0); CTEST_NEAR(hv, 42.0); \
} while (0)

static int temperature_range(void) {
    double x[2] = {0.25, 0.75}, y[2] = {0.6, 0.4};
    double r[2] = {42.0, 42.0}, hl = 42.0, hv = 42.0;
    CTEST_STATUS(lumo_thermo_eval(343.0, 100000.0, x, y, r, &hl, &hv), LUMO_THERMO_TEMPERATURE_RANGE);
    CHECK_UNCHANGED();
    CTEST_STATUS(lumo_thermo_eval(374.0, 100000.0, x, y, r, &hl, &hv), LUMO_THERMO_TEMPERATURE_RANGE);
    CHECK_UNCHANGED();
    return 0;
}

static int invalid_input(void) {
    double x[2] = {0.25, 0.75}, y[2] = {0.6, 0.4};
    double r[2] = {42.0, 42.0}, hl = 42.0, hv = 42.0;
    CTEST_STATUS(lumo_thermo_eval(NAN, 100000.0, x, y, r, &hl, &hv), LUMO_THERMO_INVALID_INPUT);
    CHECK_UNCHANGED();
    CTEST_STATUS(lumo_thermo_eval(360.0, 0.0, x, y, r, &hl, &hv), LUMO_THERMO_INVALID_INPUT);
    CHECK_UNCHANGED();
    CTEST_STATUS(lumo_thermo_eval(360.0, INFINITY, x, y, r, &hl, &hv), LUMO_THERMO_INVALID_INPUT);
    CHECK_UNCHANGED();
    CTEST_STATUS(lumo_thermo_eval(360.0, 100000.0, NULL, y, r, &hl, &hv), LUMO_THERMO_INVALID_INPUT);
    CHECK_UNCHANGED();
    return 0;
}

static int nonfinite_composition(void) {
    double x[2] = {0.25, 0.75}, y[2] = {0.6, NAN};
    double r[2] = {42.0, 42.0}, hl = 42.0, hv = 42.0;
    CTEST_STATUS(lumo_thermo_eval(360.0, 100000.0, x, y, r, &hl, &hv), LUMO_THERMO_INVALID_INPUT);
    CHECK_UNCHANGED();
    return 0;
}

static int numeric_overflow(void) {
    double x[2] = {0.25, 0.75}, y[2] = {0.6, DBL_MAX};
    double r[2] = {42.0, 42.0}, hl = 42.0, hv = 42.0;
    CTEST_STATUS(lumo_thermo_eval(360.0, 100000.0, x, y, r, &hl, &hv), LUMO_THERMO_NUMERIC_ERROR);
    CHECK_UNCHANGED();
    y[1] = 0.4;
    CTEST_STATUS(lumo_thermo_eval(360.0, DBL_MIN, x, y, r, &hl, &hv), LUMO_THERMO_NUMERIC_ERROR);
    CHECK_UNCHANGED();
    return 0;
}

int main(int argc, char **argv) {
    const struct ctest_case cases[] = {
        {"water_reference", water_reference},
        {"enthalpy", enthalpy},
        {"pressure", pressure},
        {"reference_temperature", reference_temperature},
        {"newton_iterates", newton_iterates},
        {"temperature_range", temperature_range},
        {"invalid_input", invalid_input},
        {"nonfinite_composition", nonfinite_composition},
        {"numeric_overflow", numeric_overflow}
    };
    return ctest_run(argc, argv, cases, sizeof cases / sizeof cases[0]);
}
