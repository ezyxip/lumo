#include "ctest.h"

static int success(void) {
    CTEST_NEAR(1.0, 1.0);
    CTEST_STATUS(0, 0);
    return 0;
}

static int near_failure(void) {
    CTEST_NEAR(2.0, 1.0);
    return 0;
}

static int status_failure(void) {
    CTEST_STATUS(3, 0);
    return 0;
}

static int nan_failure(void) {
    CTEST_NEAR(NAN, 1.0);
    return 0;
}

int main(int argc, char **argv) {
    const struct ctest_case cases[] = {
        {"success", success},
        {"near_failure", near_failure},
        {"status_failure", status_failure},
        {"nan_failure", nan_failure}
    };
    return ctest_run(argc, argv, cases, sizeof cases / sizeof cases[0]);
}
