#include "lumo_mesh.h"
#include "ctest.h"
#include <float.h>

static void setup(double *s, struct lumo_mesh_operating *p) {
    size_t j, f;
    memset(s, 0, LUMO_MESH_SIZE * sizeof *s);
    memset(p, 0, sizeof *p);
    p->distillate_mol_s = 4.0;
    p->reflux_ratio = 2.0;
    for (f = 0; f < LUMO_MESH_FEEDS; ++f) {
        p->feeds[f].flow_mol_s = 10.0 / LUMO_MESH_FEEDS;
        p->feeds[f].mole_fractions[0] = 0.4;
        p->feeds[f].mole_fractions[1] = 0.6;
        p->feeds[f].enthalpy_J_mol = 160.0;
    }
    for (j = 0; j < LUMO_MESH_STAGES; ++j) {
        size_t b = j * LUMO_MESH_BLOCK;
        double liquid = 8.0;
        if (j >= LUMO_MESH_FEED_0_STAGE) liquid += p->feeds[0].flow_mol_s;
#if LUMO_MESH_FEEDS > 1
        if (j >= LUMO_MESH_FEED_1_STAGE) liquid += p->feeds[1].flow_mol_s;
#endif
        s[b + LUMO_MESH_X] = s[b + LUMO_MESH_Y] = 0.4;
        s[b + LUMO_MESH_X + 1] = s[b + LUMO_MESH_Y + 1] = 0.6;
        s[b + LUMO_MESH_T] = 350.0;
        s[b + LUMO_MESH_L] = j + 1 == LUMO_MESH_STAGES ? 6.0 : liquid;
        s[b + LUMO_MESH_V] = j == 0 ? 19800.0 : 12.0;
    }
    s[LUMO_MESH_QW] = 19800.0;
}

static int zero_residual(void) {
    double s[LUMO_MESH_SIZE], r[LUMO_MESH_SIZE];
    struct lumo_mesh_operating p;
    size_t i;
    setup(s, &p);
    CTEST_STATUS(lumo_mesh_residual(s, &p, r), LUMO_MESH_OK);
    for (i = 0; i < LUMO_MESH_SIZE; ++i) {
        char name[80];
        (void)snprintf(name, sizeof name, "zero_residual row=%lu", (unsigned long)i);
        ctest_set_name(name);
        CTEST_NEAR(r[i], 0.0);
    }
    return 0;
}

static int hand_calculated(void) {
#if LUMO_MESH_STAGES == 3 && LUMO_MESH_COMPONENTS == 2 && LUMO_MESH_FEEDS == 1
    /* Direct substitution into the stage balances documented in README.
     * All K=1. Pure-component CpL=(2,4), CpV=(3,5), deltaHvap=(1000,2000).
     */
    const double s[LUMO_MESH_SIZE] = {
        0.8, 0.2, 0.7, 0.3, 350.0, 8.0, 19800.0,
        0.5, 0.5, 0.6, 0.4, 360.0, 18.0, 12.0,
        0.2, 0.8, 0.3, 0.7, 370.0, 6.0, 12.0,
        19800.0
    };
    const double expected[LUMO_MESH_SIZE] = {
        -2.4, 2.4, -0.1, 0.1, 0.0, 0.0, -1704.0,
        -2.2, 2.2, 0.1, -0.1, 0.0, 0.0, 4280.0,
        4.2, -4.2, 0.1, -0.1, 0.0, 0.0, -2568.0,
        0.0
    };
    double r[LUMO_MESH_SIZE], unused[LUMO_MESH_SIZE];
    struct lumo_mesh_operating p;
    size_t i;
    setup(unused, &p);
    p.feeds[0].enthalpy_J_mol = 200.0;
    CTEST_STATUS(lumo_mesh_residual(s, &p, r), LUMO_MESH_OK);
    for (i = 0; i < LUMO_MESH_SIZE; ++i) {
        char name[80];
        (void)snprintf(name, sizeof name, "hand_calculated row=%lu", (unsigned long)i);
        ctest_set_name(name);
        CTEST_NEAR(r[i], expected[i]);
    }
#endif
    return 0;
}

static int summations(void) {
    double s[LUMO_MESH_SIZE], r[LUMO_MESH_SIZE];
    struct lumo_mesh_operating p;
    setup(s, &p);
    s[LUMO_MESH_X] -= 0.7;
    s[LUMO_MESH_Y] += 0.2;
    CTEST_STATUS(lumo_mesh_residual(s, &p, r), LUMO_MESH_OK);
    CTEST_NEAR(r[LUMO_MESH_SX], -0.7);
    CTEST_NEAR(r[LUMO_MESH_SY], 0.2);
    CTEST_NEAR(r[LUMO_MESH_E], 0.9);
    return 0;
}

static int heat_signs(void) {
    double s[LUMO_MESH_SIZE], r[LUMO_MESH_SIZE], before[LUMO_MESH_SIZE];
    struct lumo_mesh_operating p;
    size_t i;
    setup(s, &p);
    CTEST_STATUS(lumo_mesh_residual(s, &p, before), LUMO_MESH_OK);
    s[LUMO_MESH_QD] += 7.0;
    s[LUMO_MESH_QW] += 11.0;
    CTEST_STATUS(lumo_mesh_residual(s, &p, r), LUMO_MESH_OK);
    for (i = 0; i < LUMO_MESH_SIZE; ++i) {
        double delta = i == LUMO_MESH_H ? -7.0 :
            (i == (LUMO_MESH_STAGES - 1) * LUMO_MESH_BLOCK + LUMO_MESH_H ? 11.0 : 0.0);
        CTEST_NEAR(r[i] - before[i], delta);
    }
    return 0;
}

static int overall_balances(void) {
    double s[LUMO_MESH_SIZE], r[LUMO_MESH_SIZE];
    double hl0, hv0, hln, hvn, e[LUMO_MESH_COMPONENTS];
    struct lumo_mesh_operating p;
    size_t j, i, f, last = (LUMO_MESH_STAGES - 1) * LUMO_MESH_BLOCK;
    double energy = 0.0, input_energy = 0.0;
    setup(s, &p);
    /* Deliberately not a solution: internal terms must cancel anyway. */
    for (j = 0; j < LUMO_MESH_STAGES; ++j) {
        size_t b = j * LUMO_MESH_BLOCK;
        s[b + LUMO_MESH_X] = 0.2 + 0.05 * j;
        s[b + LUMO_MESH_X + 1] = 1.0 - s[b + LUMO_MESH_X];
        s[b + LUMO_MESH_Y] = 0.4 + 0.03 * j;
        s[b + LUMO_MESH_Y + 1] = 1.0 - s[b + LUMO_MESH_Y];
        s[b + LUMO_MESH_T] = 330.0 + 5.0 * j;
        s[b + LUMO_MESH_L] = 20.0 + j;
        if (j > 0) s[b + LUMO_MESH_V] = 15.0 + j;
    }
    CTEST_STATUS(lumo_mesh_residual(s, &p, r), LUMO_MESH_OK);
    for (i = 0; i < LUMO_MESH_COMPONENTS; ++i) {
        double balance = 0.0, input = 0.0;
        for (j = 0; j < LUMO_MESH_STAGES; ++j) balance += r[j * LUMO_MESH_BLOCK + LUMO_MESH_M + i];
        for (f = 0; f < LUMO_MESH_FEEDS; ++f) input += p.feeds[f].flow_mol_s * p.feeds[f].mole_fractions[i];
        CTEST_NEAR(balance, input - p.distillate_mol_s * s[LUMO_MESH_X + i] - s[last + LUMO_MESH_L] * s[last + LUMO_MESH_X + i]);
    }
    CTEST_STATUS(lumo_thermo_eval(s[LUMO_MESH_T], LUMO_MESH_PRESSURE_PA,
        s + LUMO_MESH_X, s + LUMO_MESH_Y, e, &hl0, &hv0), LUMO_THERMO_OK);
    CTEST_STATUS(lumo_thermo_eval(s[last + LUMO_MESH_T], LUMO_MESH_PRESSURE_PA,
        s + last + LUMO_MESH_X, s + last + LUMO_MESH_Y, e, &hln, &hvn), LUMO_THERMO_OK);
    for (j = 0; j < LUMO_MESH_STAGES; ++j) energy += r[j * LUMO_MESH_BLOCK + LUMO_MESH_H];
    for (f = 0; f < LUMO_MESH_FEEDS; ++f) input_energy += p.feeds[f].flow_mol_s * p.feeds[f].enthalpy_J_mol;
    CTEST_NEAR(energy, input_energy - p.distillate_mol_s * hl0 - s[last + LUMO_MESH_L] * hln - s[LUMO_MESH_QD] + s[LUMO_MESH_QW]);
    return 0;
}

static int closure(void) {
    double s[LUMO_MESH_SIZE], r[LUMO_MESH_SIZE];
    struct lumo_mesh_operating p;
    setup(s, &p);
    p.reflux_ratio = 3.0;
    CTEST_STATUS(lumo_mesh_residual(s, &p, r), LUMO_MESH_OK);
    CTEST_NEAR(r[LUMO_MESH_CLOSURE], -4.0);
    return 0;
}

#define EXPECT_ERROR(expr, status) do { \
    CTEST_STATUS(expr, status); \
    for (i = 0; i < LUMO_MESH_SIZE; ++i) { CTEST_NEAR(r[i], 42.0); } \
} while (0)

static int errors(void) {
    double s[LUMO_MESH_SIZE], r[LUMO_MESH_SIZE];
    struct lumo_mesh_operating p;
    size_t i;
    setup(s, &p);
    for (i = 0; i < LUMO_MESH_SIZE; ++i) r[i] = 42.0;
    EXPECT_ERROR(lumo_mesh_residual(NULL, &p, r), LUMO_MESH_INVALID_INPUT);
    EXPECT_ERROR(lumo_mesh_residual(s, NULL, r), LUMO_MESH_INVALID_INPUT);
    s[LUMO_MESH_BLOCK + LUMO_MESH_T] = 290.0;
    EXPECT_ERROR(lumo_mesh_residual(s, &p, r), LUMO_MESH_TEMPERATURE_RANGE);
    setup(s, &p);
    s[LUMO_MESH_QW] = NAN;
    EXPECT_ERROR(lumo_mesh_residual(s, &p, r), LUMO_MESH_INVALID_INPUT);
    setup(s, &p);
    p.feeds[0].mole_fractions[0] = 0.8;
    EXPECT_ERROR(lumo_mesh_residual(s, &p, r), LUMO_MESH_INVALID_INPUT);
    setup(s, &p);
    p.distillate_mol_s = 10.0;
    EXPECT_ERROR(lumo_mesh_residual(s, &p, r), LUMO_MESH_INVALID_INPUT);
    setup(s, &p);
    p.reflux_ratio = INFINITY;
    EXPECT_ERROR(lumo_mesh_residual(s, &p, r), LUMO_MESH_INVALID_INPUT);
    setup(s, &p);
    p.feeds[0].enthalpy_J_mol = DBL_MAX;
    EXPECT_ERROR(lumo_mesh_residual(s, &p, r), LUMO_MESH_NUMERIC_ERROR);
    return 0;
}

int main(int argc, char **argv) {
    const struct ctest_case cases[] = {
        {"zero_residual", zero_residual}, {"hand_calculated", hand_calculated},
        {"summations", summations}, {"heat_signs", heat_signs},
        {"overall_balances", overall_balances}, {"closure", closure}, {"errors", errors}
    };
    return ctest_run(argc, argv, cases, sizeof cases / sizeof cases[0]);
}
