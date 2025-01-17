#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "system.h"

// Read initial conditions from a file 
void read_initial_conditions(double *x0, int n, const char *filename, double *t0, double *tf) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Error opening input file");
        exit(EXIT_FAILURE);
    }

    // Read t0 and tf
    if (fscanf(file, "%lf %lf", t0, tf) != 2) {
        fprintf(stderr, "Error reading t0 and tf from input file\n");
        fclose(file);
        exit(EXIT_FAILURE);
    }
    printf("Read t0 = %f, tf = %f\n", *t0, *tf);

    // Read initial conditions
    for (int i = 0; i < n; i++) {
        if (fscanf(file, "%lf", &x0[i]) != 1) {
            fprintf(stderr, "Error reading x[%d] from input file\n", i);
            fclose(file);
            exit(EXIT_FAILURE);
        }
        printf("x[%d] = %f\n", i, x0[i]);
    }

    fclose(file);
}

// Write states to a file
void write_output_csv(const char *filename, double t, const double *x, int n) {
    static int header_written = 0;  // To ensure the header is written only once
    FILE *file = fopen(filename, "a");
    if (!file) {
        perror("Error opening output file");
        exit(EXIT_FAILURE);
    }

    if (!header_written) {
        fprintf(file, "Time"); // Write header for time
        for (int i = 0; i < n; i++) {
            fprintf(file, ",x[%d]", i); // Header for each x component
        }
        fprintf(file, "\n");
        header_written = 1;
    }

    // Write time and state vector
    fprintf(file, "%.16f", t);
    for (int i = 0; i < n; i++) {
        fprintf(file, ",%.16f", x[i]);
    }
    fprintf(file, "\n");

    fclose(file);
}

// Define the right-hand side of differential equations
// This is where you specify actual system

// Define constants for the Earth and J2 perturbation
#define MU 398600.4418             // Earth's gravitational parameter [km^3/s^2]
#define R_EARTH 6378.137           // Earth's equatorial radius [km]
#define J2 1.08263e-3              // Earth's J2 value

void rhs(const double *x, double t, double *dxdt, int n) {
    // Extract equinoctial elements
    double p = x[0];     
    double f = x[1];     
    double g = x[2];     
    double h = x[3];     
    double k = x[4];     
    double L = x[5];     

    // Calculate intermediate parameters
    double sqrt_p_mu = sqrt(p / MU);
    double s_sq = 1 + h * h + k * k;
    double w = 1 + f * cos(L) + g * sin(L);
    double r = p / w;
    double h_cos_L_plus_k_sin_L = h * cos(L) + k * sin(L);
    double h_sin_L_minus_k_cos_L = h * sin(L) - k * cos(L);

    // J2 Perturbations
    double factor = MU * J2 * R_EARTH * R_EARTH / pow(r, 4);
    double delta_r = -factor * (3.0/2.0) * (1.0 - 12.0 * h_sin_L_minus_k_cos_L * h_sin_L_minus_k_cos_L / (s_sq * s_sq));
    double delta_t = -factor * 12.0 * h_sin_L_minus_k_cos_L * h_cos_L_plus_k_sin_L / (s_sq * s_sq);
    double delta_n = -factor * 6.0 * (1.0 - h * h - k * k) * h_sin_L_minus_k_cos_L / (s_sq * s_sq);

    // Compute derivatives
    dxdt[0] = 2 * sqrt_p_mu * (p / w) * delta_t; // dp/dt
    dxdt[1] = sqrt_p_mu * (delta_r * sin(L) + 
                           ((w + 1) * cos(L) + f) * delta_t / w - 
                           h_sin_L_minus_k_cos_L * g * delta_n / w); // df/dt
    dxdt[2] = sqrt_p_mu * (-delta_r * cos(L) + 
                           ((w + 1) * sin(L) + g) * delta_t / w + 
                           h_sin_L_minus_k_cos_L * g * delta_n / w); // dg/dt
    dxdt[3] = sqrt_p_mu * (s_sq * delta_n * cos(L)) / (2.0 * w); // dh/dt
    dxdt[4] = sqrt_p_mu * (s_sq * delta_n * sin(L)) / (2.0 * w); // dk/dt
    dxdt[5] = sqrt(MU * p) * pow(w / p, 2) + 
              sqrt_p_mu * h_sin_L_minus_k_cos_L * delta_n / w; // dL/dt
}

// Define constant for the Earth for 2BP without perturbations
/*#define MU 398600.4418             // Earth's gravitational parameter [km^3/s^2]

void rhs(const double *x, double t, double *dxdt, int n) {
    // Extract equinoctial elements
    double p = x[0];     
    double f = x[1];     
    double g = x[2];     
    double h = x[3];     
    double k = x[4];     
    double L = x[5];     

    // Calculate intermediate parameters
    double w = 1 + f * cos(L) + g * sin(L);

    // Compute derivatives
    dxdt[0] = 0; // dp/dt
    dxdt[1] = 0; // df/dt
    dxdt[2] = 0; // dg/dt
    dxdt[3] = 0; // dh/dt
    dxdt[4] = 0; // dk/dt
    dxdt[5] = sqrt(MU*p)*(w/p)*(w/p); // dL/dt
}*/
