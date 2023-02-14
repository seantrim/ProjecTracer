# ProjecTracer [![DOI](https://zenodo.org/badge/599500216.svg)](https://zenodo.org/badge/latestdoi/599500216)

## Contents
1. [Background](#background)
2. [Variable Definitions](#variable-definitions)
3. [Governing Equations](#governing-equations)
4. [File Description](#file-description)
    1. [Job Files](#job-files) 
    2. [H_routines](#h_routines)
    3. [ProjecTracer_Source](#projectracer_source)
    4. [plot_scripts](#plot_scripts)
5. [How to Use](#how-to-use)
6. [Legal](#legal)

## Background

2D mantle convection code that uses Lagrangian tracer particles for the advection of tempertaure and compostion. Tracer information is projected onto a fixed Eulerian grid to calculate quantites such as thermal diffusion. 

## Variable Definitions
* $\lambda$ = aspect ratio
* $x$ = horizontal position ( $x \in [0,\lambda]$ )
* $z$ = vertical position ( $z \in [0,1]$ , increasing upward)
* $t$ = time
* $\psi$ = stream function
* $\mathbf{v}$ = velocity
* $C$ = composition
* $T$ = temperature
* $H$ = internal heating rate
* $\eta$ = viscosity

## Governing Equations

### Conservation of Mass and Momentum: the Biharmonic Equation
$\left(\frac{\partial^2}{\partial x^2}-\frac{\partial^2}{\partial z^2}\right) \left[ \eta \left(\frac{\partial^2 \psi}{\partial x^2}-\frac{\partial^2 \psi}{\partial z^2}\right) \right] + 4\frac{\partial^2}{\partial x \partial z} \left[ \eta \frac{\partial^2 \psi}{\partial x \partial z} \right] = Ra_T \frac{\partial T}{\partial x} - Ra_C \frac{\partial C}{\partial x}$

### Conservation of Eneregy
$\frac{\partial T}{\partial t}+\mathbf{v} \cdot \nabla T = \nabla^2 T+H$

### Compositional Transport
$\frac{\partial C}{\partial t}+\mathbf{v} \cdot \nabla C = 0$

## File Description

### Job Files
[compile_script.sh](/compile_script.sh)
* Linux script for compiling Fortran source files

[ProjecTracer](/ProjecTracer)
* Exectutable for ProjecTracer code

[run_job.sh](/run_job.sh)
* Linux script for submitting local OpenMP jobs

### [H_routines](/H_routines)
* Contains Fortran source files for the calculation of $H$ for certain problems
* See https://github.com/seantrim/exact-thermochem-solution for details

### [ProjecTracer_Source](/ProjecTracer_Source)
* Fortran source code files for ProjecTracer code

### [plot_scripts](/plot_scripts)
* Gnuplot scripts for analysis and visualization

## How to Use

1. Specify physical parameters, boundary conditions, and code options (e.g., aspect ratio, Rayleigh numbers, Courant number, etc.) in [inputs.nml](/inputs.nml)
    * Entries in the Fortran namelist [inputs.nml](/inputs.nml) are fully commented
2. Specify initial conditions for $T$ and $C$ in [initial_conditions.f90](/ProjecTracer_Source/initial_conditions.f90)
    * Variable $H$ is specified in `compute_heating` subroutine in [energy.f90](/ProjecTracer_Source/energy.f90)
3. Compile source code using [compile_script.sh](/compile_script.sh)
    * Linux: `$ source compile_script`
        * This produces the executable [ProjecTracer](/ProjecTracer)
    * gfortran 11.3.0 or later is recommended
    * Other compilers may be possible but results should be tested
4. Run [ProjecTracer](/ProjecTracer) from the command line
    * Specify number of OpenMP processes in [run_job.sh](/run_job.sh) then run the script
    * Linux: `$ source run_job.sh`
        * Job will run in the background

## Legal

This repository is subject to the GPLv3 [license](/LICENSE).

The routines contained within [xelbdj2_all_routines.f90](/xelbdj2_all_routines.f90) and [xgscd_routines.f90](/xgscd_routines.f90) are adapted from routines by Toshio Fukushima available under the CC BY-SA 4.0 license. Original versions of these routines can be found at http://dx.doi.org/10.13140/RG.2.2.27011.66085 and https://www.researchgate.net/publication/233903220_xgscdtxt_Fortran_program_package_to_compute_the_Jacobian_elliptic_functions_snum_cnum_dnum.
