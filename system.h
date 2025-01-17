#ifndef SYSTEM_H
#define SYSTEM_H

// System parameters and configuration
#define SYSTEM_DIMENSION 6
#define DEFAULT_TOLERANCE 1e-15

// Function declarations
void rhs(const double *x, double t, double *dxdt, int n);
void read_initial_conditions(double *x0, int n, const char *filename, double *t0, double *tf);
void write_output_csv(const char *filename, double t, const double *x, int n);

#endif
