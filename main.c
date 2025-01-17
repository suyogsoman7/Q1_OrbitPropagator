#include <stdio.h>
#include <stdlib.h>
#include "routines.h"
#include "system.h"

int main() {
    const char *input_file = "input.txt";
    double t0, tf; // Initialize t0 and tf
    double *x0 = (double *)malloc(SYSTEM_DIMENSION * sizeof(double));
    double *xf = (double *)malloc(SYSTEM_DIMENSION * sizeof(double));

    // Read initial conditions, including t0 and tf, from the input file
    read_initial_conditions(x0, SYSTEM_DIMENSION, input_file, &t0, &tf);

    // Integrate from t0 to tf
    integrate_rkf78(t0, x0, tf, DEFAULT_TOLERANCE, xf, SYSTEM_DIMENSION, rhs);

    // Write the output
    write_output_csv("trajectory.csv", tf, xf, SYSTEM_DIMENSION); // CSV output


    free(x0);
    free(xf);
    return 0;
}