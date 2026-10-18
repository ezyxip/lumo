/* Test-only finite-difference Newton probe, not a production solver. */
#include "lumo_mesh.h"
#include "ctest.h"
#include "reference.h"

#define N LUMO_MESH_SIZE

static double column_scale(size_t i) {
    if (i == LUMO_MESH_QD || i == LUMO_MESH_QW) return 1e7;
    if (i % LUMO_MESH_BLOCK == LUMO_MESH_T) return 350.0;
    if (i % LUMO_MESH_BLOCK == LUMO_MESH_L || i % LUMO_MESH_BLOCK == LUMO_MESH_V) return 100.0;
    return 1.0;
}

static double row_scale(size_t i) {
    if (i == LUMO_MESH_CLOSURE) return 100.0;
    if (i % LUMO_MESH_BLOCK == LUMO_MESH_H) return 1e7;
    if (i % LUMO_MESH_BLOCK < LUMO_MESH_COMPONENTS) return 100.0;
    return 1.0;
}

static int evaluate(const double *u, double *r) {
    double state[N];
    size_t i;
    int status;
    for (i = 0; i < N; ++i) state[i] = u[i] * column_scale(i);
    status = lumo_mesh_residual(state, &reference_operating, r);
    if (status != LUMO_MESH_OK) return status;
    for (i = 0; i < N; ++i) r[i] /= row_scale(i);
    return LUMO_MESH_OK;
}

static int jacobian(const double *u, double matrix[N][N]) {
    size_t i, k;
    for (k = 0; k < N; ++k) {
        double plus[N], minus[N], rp[N], rm[N];
        double step = 1e-6 * fmax(1.0, fabs(u[k]));
        int status;
        memcpy(plus, u, sizeof plus);
        memcpy(minus, u, sizeof minus);
        plus[k] += step;
        minus[k] -= step;
        status = evaluate(plus, rp);
        if (status != LUMO_MESH_OK) return status;
        status = evaluate(minus, rm);
        if (status != LUMO_MESH_OK) return status;
        for (i = 0; i < N; ++i) matrix[i][k] = (rp[i] - rm[i]) / (2.0 * step);
    }
    return LUMO_MESH_OK;
}

/* Complete-pivot elimination of the scaled Jacobian. A pivot is significant
 * above 1e-9 times the largest entry of the original matrix.
 */
static int matrix_rank(double a[N][N]) {
    size_t i, j, k;
    double largest = 0.0;
    for (i = 0; i < N; ++i) for (j = 0; j < N; ++j) largest = fmax(largest, fabs(a[i][j]));
    for (k = 0; k < N; ++k) {
        size_t row = k, col = k;
        double pivot = 0.0;
        for (i = k; i < N; ++i) for (j = k; j < N; ++j) {
            if (fabs(a[i][j]) > pivot) { pivot = fabs(a[i][j]); row = i; col = j; }
        }
        if (pivot <= 1e-9 * largest) return (int)k;
        for (j = 0; j < N; ++j) { double tmp = a[k][j]; a[k][j] = a[row][j]; a[row][j] = tmp; }
        for (i = 0; i < N; ++i) { double tmp = a[i][k]; a[i][k] = a[i][col]; a[i][col] = tmp; }
        for (i = k+1; i < N; ++i) {
            double factor = a[i][k] / a[k][k];
            for (j = k+1; j < N; ++j) a[i][j] -= factor * a[k][j];
            a[i][k] = 0.0;
        }
    }
    return N;
}

static int reference_residual(void) {
    double r[N];
    size_t i;
    CTEST_STATUS(lumo_mesh_residual(reference_state, &reference_operating, r), LUMO_MESH_OK);
    for (i = 0; i < N; ++i) CTEST_NEAR(r[i] / row_scale(i), 0.0);
    return 0;
}

static int temperature_dependence(void) {
    double x[2] = {0.5, 0.5}, y[2] = {0.5, 0.5};
    double cold[2], warm[2], hl, hv;
    CTEST_STATUS(lumo_thermo_eval(360.0, LUMO_MESH_PRESSURE_PA, x, y, cold, &hl, &hv), LUMO_THERMO_OK);
    CTEST_STATUS(lumo_thermo_eval(370.0, LUMO_MESH_PRESSURE_PA, x, y, warm, &hl, &hv), LUMO_THERMO_OK);
    CTEST_STATUS(warm[0] < cold[0], 1);
    CTEST_STATUS(warm[1] < cold[1], 1);
    CTEST_STATUS(fabs(cold[0] - cold[1]) > 0.01, 1);
    return 0;
}

static int jacobian_rank(void) {
    double u[N], a[N][N];
    const double perturbations[] = {0.0, 1e-4, -1e-4};
    size_t i, trial;
    for (trial = 0; trial < 3; ++trial) {
        char name[80];
        (void)snprintf(name, sizeof name, "jacobian_rank perturbation=%g", perturbations[trial]);
        ctest_set_name(name);
        for (i = 0; i < N; ++i) u[i] = reference_state[i] / column_scale(i) + (i % 2 ? 1 : -1) * perturbations[trial];
        CTEST_STATUS(jacobian(u, a), LUMO_MESH_OK);
        CTEST_STATUS(matrix_rank(a), N);
    }
    return 0;
}

static int linear_solve(double a[N][N], double *b, double *x) {
    size_t i, j, k;
    for (k = 0; k < N; ++k) {
        size_t pivot = k;
        for (i = k+1; i < N; ++i) if (fabs(a[i][k]) > fabs(a[pivot][k])) pivot = i;
        if (fabs(a[pivot][k]) < 1e-12) return 1;
        for (j = k; j < N; ++j) { double tmp = a[k][j]; a[k][j] = a[pivot][j]; a[pivot][j] = tmp; }
        { double tmp = b[k]; b[k] = b[pivot]; b[pivot] = tmp; }
        for (i = k+1; i < N; ++i) {
            double factor = a[i][k] / a[k][k];
            for (j = k+1; j < N; ++j) a[i][j] -= factor * a[k][j];
            a[i][k] = 0.0;
            b[i] -= factor * b[k];
        }
    }
    for (i = N; i-- > 0;) {
        double value = b[i];
        for (j = i+1; j < N; ++j) value -= a[i][j] * x[j];
        x[i] = value / a[i][i];
    }
    return 0;
}

static double norm(const double *r) {
    double largest = 0.0;
    size_t i;
    for (i = 0; i < N; ++i) largest = fmax(largest, fabs(r[i]));
    return largest;
}

static int newton_reference(void) {
    double u[N], r[N];
    const double perturbations[] = {0.002, -0.002};
    size_t i, trial;
    char name[80];
    for (trial = 0; trial < 2; ++trial) {
        int iteration;
        (void)snprintf(name, sizeof name, "newton_reference trial=%lu", (unsigned long)trial);
        ctest_set_name(name);
        for (i = 0; i < N; ++i) u[i] = reference_state[i] / column_scale(i) + (i % 2 ? 1 : -1) * perturbations[trial];
        for (iteration = 0; iteration < 20; ++iteration) {
            double a[N][N], rhs[N], delta[N], step;
            int accepted = 0;
            CTEST_STATUS(evaluate(u, r), LUMO_MESH_OK);
            if (norm(r) < 1e-12) break;
            CTEST_STATUS(jacobian(u, a), LUMO_MESH_OK);
            for (i = 0; i < N; ++i) rhs[i] = -r[i];
            CTEST_STATUS(linear_solve(a, rhs, delta), 0);
            for (step = 1.0; step >= 1.0/4096.0; step *= 0.5) {
                double next[N], next_r[N];
                for (i = 0; i < N; ++i) next[i] = u[i] + step * delta[i];
                if (evaluate(next, next_r) == LUMO_MESH_OK && norm(next_r) < norm(r)) {
                    memcpy(u, next, sizeof u);
                    accepted = 1;
                    break;
                }
            }
            CTEST_STATUS(accepted, 1);
        }
        CTEST_STATUS(iteration < 20, 1);
        CTEST_STATUS(evaluate(u, r), LUMO_MESH_OK);
        CTEST_NEAR(norm(r), 0.0);
        for (i = 0; i < N; ++i) {
            (void)snprintf(name, sizeof name, "newton_reference trial=%lu variable=%lu", (unsigned long)trial, (unsigned long)i);
            ctest_set_name(name);
            CTEST_NEAR(u[i], reference_state[i] / column_scale(i));
        }
    }
    return 0;
}

int main(int argc, char **argv) {
    const struct ctest_case cases[] = {
        {"reference_residual", reference_residual}, {"temperature_dependence", temperature_dependence},
        {"jacobian_rank", jacobian_rank}, {"newton_reference", newton_reference}
    };
    return ctest_run(argc, argv, cases, sizeof cases / sizeof cases[0]);
}
