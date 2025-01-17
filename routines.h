// routines.h
#ifndef ROUTINES_H
#define ROUTINES_H

// Function pointer type for RHS
typedef void (*RHSFunction)(const double *x, double t, double *dxdt, int n);

// Function declarations for the integrator
void vector_copy(double *dest, const double *src, int n);
void vector_add_scaled(double *result, const double *y, const double *x, double a, int n);
void rkf78_step(double t, const double *x, double h, double *x_next, double *error, int n, RHSFunction rhs);
void integrate_rkf78(double t0, const double *x0, double tf, double tolerance, double *xf, int n, RHSFunction rhs);

#endif