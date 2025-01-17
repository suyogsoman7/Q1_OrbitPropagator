# Orbit Propagation using RKF78 Integration

This repository contains code for orbit propagation using the Runge-Kutta-Fehlberg 7(8) (RKF78) numerical integration method in C language, and plotting and comparison using MATLAB.

## Table of Contents
- [Introduction](#introduction)
- [Usage](#usage)
- [References](#references)


## Introduction
The RKF78 integration method is an adaptive step-size numerical technique used to solve ordinary differential equations. This code is specifically designed for orbit propagation where the dynamics considered are two-body problem with J2 perturbation. The code is written in a way that it can be easily modified to include other perturbations. The equations are in modified equinoctial orbital elements. The file "system.c" can be modified to include high-fidelity dynamics equations.  

The code written in C language is for propagation of initial equinoctial state as given in "input.txt" file. The output is written in "trajectory.csv" file. This file contains the time steps and corresponding modified equinoctial states. This file is taken as input to MATLAB code for plotting the trajectory and comparing with GMAT results.

## Usage
1. If the initial state is in classical orbital elements, then convert them to modified equinoctial elements using the MATLAB code "classical_to_equinoctial.m". 

Open "input.txt" file and enter the initial state along with the time of propagation. The format is,

```
t0 tf
p f g h k L
```
where, t0 is initial time in seconds, tf is final time in seconds, p (km), f, g, h, k, L (rad.) are the modified equinoctial orbital elements.

2. Compile the code using
```
make
```
3. Run the code using
```
./solver
```

4. The output is written in "trajectory.csv" file. This file contains the time steps and corresponding modified equinoctial states.

5. Run the MATLAB code "meqcplot.m" to get the trajectory in cartesian coordinates, and plot it along with animation.

6. To compare results with GMAT, paste the results from GMAT in cartesian coordinates in variable "gmat_end" in "diff_gmat.m" file. Run this code to get the error in position and velocity.

## References
1. [Modified Equinoctal Elements - JPL](https://spsweb.fltops.jpl.nasa.gov/portaldataops/mpg/MPG_Docs/Source%20Docs/EquinoctalElements-modified.pdf)
2. [RKF78 MATLAB code by Meysam Mahooti](https://in.mathworks.com/matlabcentral/fileexchange/61130-runge-kutta-fehlberg-rkf78)

