#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "routines.h"
#include "system.h"

// Constants for the integrator
#define ATTEMPTS 12
#define MIN_SCALE_FACTOR 0.125
#define MAX_SCALE_FACTOR 4.0
#define ERR_EXPONENT 1.0/7.0

// RKF78 coefficients
const double c_1_11 = 41.0/840.0;
const double c6 = 34.0/105.0;
const double c_7_8 = 9.0/35.0;
const double c_9_10 = 9.0/280.0;
const double err_factor = -41.0/840.0;

// Helper function to copy arrays
void vector_copy(double *dest, const double *src, int n) {
    for (int i = 0; i < n; i++) {
        dest[i] = src[i];
    }
}

// Helper function to add scaled vectors: result = y + a*x
void vector_add_scaled(double *result, const double *y, const double *x, double a, int n) {
    for (int i = 0; i < n; i++) {
        result[i] = y[i] + a * x[i];
    }
}

// Function to perform one RKF78 step for a vector system
void rkf78_step(double t, const double *x, double h, double *x_next, double *error, int n, RHSFunction rhs) {
    // Allocate memory for temporary arrays
    double *k1 = (double *)malloc(n * sizeof(double));
    double *k2 = (double *)malloc(n * sizeof(double));
    double *k3 = (double *)malloc(n * sizeof(double));
    double *k4 = (double *)malloc(n * sizeof(double));
    double *k5 = (double *)malloc(n * sizeof(double));
    double *k6 = (double *)malloc(n * sizeof(double));
    double *k7 = (double *)malloc(n * sizeof(double));
    double *k8 = (double *)malloc(n * sizeof(double));
    double *k9 = (double *)malloc(n * sizeof(double));
    double *k10 = (double *)malloc(n * sizeof(double));
    double *k11 = (double *)malloc(n * sizeof(double));
    double *k12 = (double *)malloc(n * sizeof(double));
    double *k13 = (double *)malloc(n * sizeof(double));
    double *temp = (double *)malloc(n * sizeof(double));

    // Time coefficients
    const double a2 = 2.0/27.0, a3 = 1.0/9.0, a4 = 1.0/6.0;
    const double a5 = 5.0/12.0, a6 = 1.0/2.0, a7 = 5.0/6.0;
    const double a8 = 1.0/6.0, a9 = 2.0/3.0, a10 = 1.0/3.0;

    // Calculate k1
    rhs(x, t, k1, n);
    for (int i = 0; i < n; i++) k1[i] *= h;

    // Calculate k2
    vector_add_scaled(temp, x, k1, 2.0/27.0, n);
    rhs(temp, t + a2*h, k2, n);
    for (int i = 0; i < n; i++) k2[i] *= h;

    // Calculate k3
    for (int i = 0; i < n; i++) {
        temp[i] = x[i] + (1.0/36.0)*k1[i] + (3.0/36.0)*k2[i];
    }
    rhs(temp, t + a3*h, k3, n);
    for (int i = 0; i < n; i++) k3[i] *= h;

    // Calculate k4
    for (int i = 0; i < n; i++) {
        temp[i] = x[i] + (1.0/24.0)*k1[i] + (3.0/24.0)*k3[i];
    }
    rhs(temp, t + a4*h, k4, n);
    for (int i = 0; i < n; i++) k4[i] *= h;

    // Calculate k5
    for (int i = 0; i < n; i++) {
        temp[i] = x[i] + (20.0/48.0)*k1[i] - (75.0/48.0)*k3[i] + (75.0/48.0)*k4[i];
    }
    rhs(temp, t + a5*h, k5, n);
    for (int i = 0; i < n; i++) k5[i] *= h;

    // Calculate k6
    for (int i = 0; i < n; i++) {
        temp[i] = x[i] + (1.0/20.0)*k1[i] + (5.0/20.0)*k4[i] + (4.0/20.0)*k5[i];
    }
    rhs(temp, t + a6*h, k6, n);
    for (int i = 0; i < n; i++) k6[i] *= h;

    // Calculate k7
    for (int i = 0; i < n; i++) {
        temp[i] = x[i] + (-25.0/108.0)*k1[i] + (125.0/108.0)*k4[i] - 
                  (260.0/108.0)*k5[i] + (250.0/108.0)*k6[i];
    }
    rhs(temp, t + a7*h, k7, n);
    for (int i = 0; i < n; i++) k7[i] *= h;

    // Calculate k8
    for (int i = 0; i < n; i++) {
        temp[i] = x[i] + (31.0/300.0)*k1[i] + (61.0/225.0)*k5[i] - 
                  (2.0/9.0)*k6[i] + (13.0/900.0)*k7[i];
    }
    rhs(temp, t + a8*h, k8, n);
    for (int i = 0; i < n; i++) k8[i] *= h;

    // Calculate k9
    for (int i = 0; i < n; i++) {
        temp[i] = x[i] + 2.0*k1[i] - (53.0/6.0)*k4[i] + (704.0/45.0)*k5[i] - 
                  (107.0/9.0)*k6[i] + (67.0/90.0)*k7[i] + 3.0*k8[i];
    }
    rhs(temp, t + a9*h, k9, n);
    for (int i = 0; i < n; i++) k9[i] *= h;

    // Calculate k10
    for (int i = 0; i < n; i++) {
        temp[i] = x[i] + (-91.0/108.0)*k1[i] + (23.0/108.0)*k4[i] - 
                  (976.0/135.0)*k5[i] + (311.0/54.0)*k6[i] - (19.0/60.0)*k7[i] + 
                  (17.0/6.0)*k8[i] - (1.0/12.0)*k9[i];
    }
    rhs(temp, t + a10*h, k10, n);
    for (int i = 0; i < n; i++) k10[i] *= h;

    // Calculate k11
    for (int i = 0; i < n; i++) {
        temp[i] = x[i] + (2383.0/4100.0)*k1[i] - (341.0/164.0)*k4[i] + 
                  (4496.0/1025.0)*k5[i] - (301.0/82.0)*k6[i] + (2133.0/4100.0)*k7[i] + 
                  (45.0/82.0)*k8[i] + (45.0/164.0)*k9[i] + (18.0/41.0)*k10[i];
    }
    rhs(temp, t + h, k11, n);
    for (int i = 0; i < n; i++) k11[i] *= h;

    // Calculate k12
    for (int i = 0; i < n; i++) {
        temp[i] = x[i] + (3.0/205.0)*k1[i] - (6.0/41.0)*k6[i] - (3.0/205.0)*k7[i] - 
                  (3.0/41.0)*k8[i] + (3.0/41.0)*k9[i] + (6.0/41.0)*k10[i];
    }
    rhs(temp, t, k12, n);
    for (int i = 0; i < n; i++) k12[i] *= h;

    // Calculate k13
    for (int i = 0; i < n; i++) {
        temp[i] = x[i] + (-1777.0/4100.0)*k1[i] - (341.0/164.0)*k4[i] + 
                  (4496.0/1025.0)*k5[i] - (289.0/82.0)*k6[i] + (2193.0/4100.0)*k7[i] + 
                  (51.0/82.0)*k8[i] + (33.0/164.0)*k9[i] + (12.0/41.0)*k10[i] + k12[i];
    }
    rhs(temp, t + h, k13, n);
    for (int i = 0; i < n; i++) k13[i] *= h;

    // Calculate final result
    for (int i = 0; i < n; i++) {
        x_next[i] = x[i] + (c_1_11 * (k1[i] + k11[i]) + c6 * k6[i] + 
                           c_7_8 * (k7[i] + k8[i]) + c_9_10 * (k9[i] + k10[i]));
    }

    // Calculate error
    *error = 0.0;
    for (int i = 0; i < n; i++) {
        //double local_error = fabs(err_factor * (k1[i] + k11[i] - k12[i] - k13[i]));
        double local_error = err_factor * (k1[i] + k11[i] - k12[i] - k13[i]);
        if (local_error > *error) *error = local_error;
    }

    // Free allocated memory
    free(k1); free(k2); free(k3); free(k4); free(k5); free(k6); free(k7);
    free(k8); free(k9); free(k10); free(k11); free(k12); free(k13); free(temp);
}

// Main integration function for vector system
void integrate_rkf78(double t0, const double *x0, double tf, double tolerance,
                     double *xf, int n, RHSFunction rhs) {
    double t = t0;
    double *x = (double *)malloc(n * sizeof(double));
    double *x_next = (double *)malloc(n * sizeof(double));
    vector_copy(x, x0, n);

    double h = (tf - t0) / 100000.0; // Initial step size

    // Open output file for writing
    const char *output_file = "trajectory.csv";
    remove(output_file); // Ensure a fresh file each run

    // Save initial conditions
    write_output_csv(output_file, t, x, n);

    while (t < tf) {
        if (t + h > tf) {
            h = tf - t; // Adjust final step size
        }

        double error;
        int attempts;
        double scale = 1.0;

        // Try to take a step
        for (attempts = 0; attempts < ATTEMPTS; attempts++) {
            rkf78_step(t, x, h, x_next, &error, n, rhs);

            if (error == 0.0) {
                scale = MAX_SCALE_FACTOR;
                break;
            }

            double max_x = tolerance;
            for (int i = 0; i < n; i++) {
                if (fabs(x[i]) > max_x) max_x = fabs(x[i]);
            }

            scale = 0.8 * pow(tolerance * max_x / error, ERR_EXPONENT);
            scale = fmin(fmax(scale, MIN_SCALE_FACTOR), MAX_SCALE_FACTOR);

            if (error < tolerance * max_x) {
                break;
            }

            h *= scale;
        }

        if (attempts < ATTEMPTS) {
            vector_copy(x, x_next, n);
            t += h;
            h *= scale;

            // Write the state to the CSV file
            write_output_csv(output_file, t, x, n);
        } else {
            h *= MIN_SCALE_FACTOR;
        }
    }

    // Copy final result
    vector_copy(xf, x, n);

    free(x);
    free(x_next);
}
