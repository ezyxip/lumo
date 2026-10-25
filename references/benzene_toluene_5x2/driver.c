/* Validation harness only. IO belongs here, never in generated model code. */
#include "lumo_prepare.h"
#include "lumo_initialguess.h"
#include "lumo_solver.h"
#include "input.h"
#include <stdio.h>

static double state[LUMO_MESH_SIZE];
static double workspace[LUMO_SOLVER_WORKSPACE_DOUBLES];
static double residual[LUMO_MESH_SIZE];

static void observe(const struct lumo_solver_iteration *record, void *context) {
    size_t *count = context;
    (void)record;
    ++*count;
}

static void print_array(const double *values, size_t count) {
    size_t i;
    putchar('[');
    for (i = 0; i < count; ++i) printf("%s%.17g", i ? "," : "", values[i]);
    putchar(']');
}

int main(void) {
    struct lumo_mesh_operating operating;
    struct lumo_solver_options options;
    struct lumo_solver_report report;
    size_t events = 0;
    int status = lumo_prepare_operating(&input, &operating);
    if (status != LUMO_PREPARE_OK) {
        fprintf(stderr, "prepare failed: %d\n", status);
        return 1;
    }
    status = lumo_initialguess_linear_T_x(&operating, state);
    if (status != LUMO_MESH_OK) {
        fprintf(stderr, "initial guess failed: %d\n", status);
        return 1;
    }
    lumo_solver_default_options(&options);
    options.observer = observe;
    options.observer_context = &events;
    status = lumo_solver_solve(&operating, state, &options, workspace,
                              LUMO_SOLVER_WORKSPACE_DOUBLES, &report);
    if (status != LUMO_SOLVER_CONVERGED) {
        fprintf(stderr, "solve failed: %d, model status: %d, iterations: %d\n",
                status, report.model_status, report.iterations);
        return 1;
    }
    status = lumo_mesh_residual(state, &operating, residual);
    if (status != LUMO_MESH_OK) return 1;
    printf("{\"status\":%d,\"iterations\":%d,\"backtracks\":%zu,"
           "\"residual_norm\":%.17g,\"events\":%zu,\"diagnostics\":%d,"
           "\"feed_enthalpy_J_mol\":%.17g,\"state\":",
           report.status, report.iterations, report.backtracks,
           report.residual_norm, events, LUMO_SOLVER_HAS_DIAGNOSTICS,
           operating.feeds[0].enthalpy_J_mol);
    print_array(state, LUMO_MESH_SIZE);
    printf(",\"residual\":");
    print_array(residual, LUMO_MESH_SIZE);
    puts("}");
    return 0;
}
