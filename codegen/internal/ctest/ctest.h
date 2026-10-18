/* Shared C99 assertions and named test runner for generated code. */
#ifndef LUMO_CTEST_H
#define LUMO_CTEST_H

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

static const char *ctest_name;

static inline void ctest_set_name(const char *name) {
    ctest_name = name;
}

#define CTEST_STATUS(expr, expected) do { \
    int actual_status = (expr); \
    int expected_status = (expected); \
    if (actual_status != expected_status) { \
        fprintf(stderr, "%s %s:%d: %s actual=%d expected=%d (%s)\n", \
            ctest_name, __FILE__, __LINE__, #expr, actual_status, expected_status, #expected); \
        return 1; \
    } \
} while (0)

#define CTEST_NEAR(expr, expected) do { \
    double actual_value = (expr); \
    double expected_value = (expected); \
    double error = fabs(actual_value - expected_value); \
    double tolerance = 1e-10 * fmax(1.0, fabs(expected_value)); \
    if (!isfinite(actual_value) || !isfinite(expected_value) || error > tolerance) { \
        fprintf(stderr, "%s %s:%d: %s actual=%.17g expected=%.17g error=%.17g tolerance=%.17g\n", \
            ctest_name, __FILE__, __LINE__, #expr, actual_value, expected_value, error, tolerance); \
        return 1; \
    } \
} while (0)

struct ctest_case {
    const char *name;
    int (*run)(void);
};

static inline int ctest_run(int argc, char **argv,
                            const struct ctest_case *cases, size_t count) {
    size_t i;
    if (argc != 2) {
        fprintf(stderr, "usage: %s <test-name>\n", argv[0]);
        return 2;
    }
    for (i = 0; i < count; ++i) {
        if (strcmp(argv[1], cases[i].name) == 0) {
            ctest_set_name(cases[i].name);
            return cases[i].run();
        }
    }
    fprintf(stderr, "unknown C test: %s\n", argv[1]);
    return 2;
}

#endif
