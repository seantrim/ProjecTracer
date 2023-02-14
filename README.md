# ProjecTracer

## Contents
1. [Background](#background)
2. [Variable Definitions](#variable-definitions)
3. [File Description](#file-description)
    1. [H_routines](#H_routines)
4. [How to Use](#how-to-use)
5. [Legal](#legal)

## Background

2D mantle convection code that uses Lagrangian tracer particles for the advection of tempertaure and compostion. Tracer information is projected onto a fixed Eulerian grid to calculate quantites such as thermal diffusion. 

## Variable Definitions
* $\lambda$ = aspect ratio
* $x$ = horizontal position ( $x \in [0,\lambda]$ )
* $z$ = vertical position ( $z \in [0,1]$ , increasing upward)
* $t$ = time
* $C$ = composition
* $T$ = temperature
* $H$ = internal heating rate
* $\eta$ = viscosity

## Governing Equations

## File Description

### [H_routines](/H_routines)
* Contains Fortran source files for the calculation of $H$ for certain problems
* See https://github.com/seantrim/exact-thermochem-solution for details

## How to Use

## Legal

This repository is subject to the GPLv3 [license](/LICENSE).

The routines contained within [xelbdj2_all_routines.f90](/xelbdj2_all_routines.f90) and [xgscd_routines.f90](/xgscd_routines.f90) are adapted from routines by Toshio Fukushima available under the CC BY-SA 4.0 license. Original versions of these routines can be found at http://dx.doi.org/10.13140/RG.2.2.27011.66085 and https://www.researchgate.net/publication/233903220_xgscdtxt_Fortran_program_package_to_compute_the_Jacobian_elliptic_functions_snum_cnum_dnum.
