#include "lumo_jacobian.h"
#include "ctest.h"
#define N LUMO_JACOBIAN_SIZE

/* A synthetic square system with exact known derivatives, no MESH/thermo. */
static int residual(const double *state, const void *context, double *out) {
    double a = context == NULL ? 1.0 : *(const double *)context;
    size_t i;
    if (state[0] < 0.0 || state[0] > 2.0) return LUMO_JACOBIAN_DOMAIN_ERROR;
    for (i = 0; i < N; ++i) out[i] = a * state[i] * state[i] + 2.0 * state[(i+1)%N];
    return LUMO_JACOBIAN_OK;
}

static int polynomial(void) {
    double state[N], original[N], scales[N], matrix[N*N], a = 3.0;
    struct lumo_jacobian_numeric_options options = {1e-4, scales};
    lumo_jacobian_fn compute = lumo_jacobian_numeric;
    size_t i, j, trial;
    for (i = 0; i < N; ++i) { state[i] = 0.5 + 0.1*i; scales[i] = 2.0; }
    memcpy(original, state, sizeof state);
    for (trial = 0; trial < 2; ++trial) {
        CTEST_STATUS(compute(state, residual, trial ? &a : NULL, trial ? &options : NULL, matrix), LUMO_JACOBIAN_OK);
        for (i = 0; i < N; ++i) for (j = 0; j < N; ++j) {
            double expected = (i == j ? 2.0*(trial ? a : 1.0)*state[i] : 0.0) + ((i+1)%N == j ? 2.0 : 0.0);
            CTEST_NEAR(matrix[i*N+j], expected);
        }
    }
    CTEST_STATUS(memcmp(original, state, sizeof state), 0);
    return 0;
}

static int domain(void) {
    double state[N], matrix[N*N];
    struct lumo_jacobian_numeric_options options = {1e-4, NULL};
    size_t i, trial;
    for (i = 0; i < N; ++i) state[i] = 0.5;
    for (trial = 0; trial < 2; ++trial) {
        state[0] = trial ? 2.0 : 0.0;
        CTEST_STATUS(lumo_jacobian_numeric(state, residual, NULL, &options, matrix), LUMO_JACOBIAN_OK);
        for (i = 0; i < N; ++i) {
            double expected = (i == 0 ? 2.0*state[0] : 0.0) + ((i+1)%N == 0 ? 2.0 : 0.0);
            CTEST_STATUS(fabs(matrix[i*N]-expected) < 1e-3, 1);
        }
    }
    options.relative_step = 4.0;
    CTEST_STATUS(lumo_jacobian_numeric(state, residual, NULL, &options, matrix), LUMO_JACOBIAN_DOMAIN_ERROR);
    return 0;
}

static int failed_residual(const double *state, const void *context, double *out) {
    (void)state; (void)context; (void)out;
    return LUMO_JACOBIAN_NUMERIC_ERROR;
}

static int errors(void) {
    double state[N], matrix[N*N];
    size_t i;
    for (i = 0; i < N; ++i) state[i] = 0.5;
    CTEST_STATUS(lumo_jacobian_numeric(state, NULL, NULL, NULL, matrix), LUMO_JACOBIAN_INVALID_INPUT);
    CTEST_STATUS(lumo_jacobian_numeric(state, failed_residual, NULL, NULL, matrix), LUMO_JACOBIAN_NUMERIC_ERROR);
    state[0] = NAN;
    CTEST_STATUS(lumo_jacobian_numeric(state, residual, NULL, NULL, matrix), LUMO_JACOBIAN_INVALID_INPUT);
    return 0;
}

int main(int argc, char **argv) {
    const struct ctest_case cases[] = {{"polynomial", polynomial}, {"domain", domain}, {"errors", errors}};
    return ctest_run(argc, argv, cases, sizeof cases / sizeof cases[0]);
}
