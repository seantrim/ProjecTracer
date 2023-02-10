#!/bin/bash

####Set up environoment variables for oneMKL libraries if needed:
#source /opt/intel/oneapi/setvars.sh

####Examples using gfortran

##use LAPACK with dynamic linking (no external routines for computing H)
#gfortran -O3 -fopenmp -flto *.f90 -lblas -llapack -o ProjecTracer ##use LAPACK with dynamic linking

##use LAPACK with static linking and include external routines for computing H (in H_routines directory -- must be supplied by the user)
gfortran -O3 -fopenmp -flto *.f90 H_routines/*.f90 -lblas -llapack -o ProjecTracer

##use LAPACK with static linking and include external routines for computing H (in H_routines directory -- must be supplied by the user)
#gfortran -static -O3 -fopenmp -flto *.f90 H_routines/*.f90 -lblas -llapack -o ProjecTracer

##Enable debug options with dynamic linking, LAPACK, and include routines for computing H (in H_routines directory -- must be supplied by the user)
#gfortran -fcheck=bounds -Og -ffpe-trap=invalid,zero,overflow,underflow -ggdb3 -fopenmp *.f90 H_routines/*.f90 -lblas -llapack -o ProjecTracer

##use oneMKL with dynamic linking and include external routines for computing H (in H_routines directory -- must be supplied by the user)
#gfortran -O3 -fopenmp -flto -m64 *.f90 H_routines/*.f90 -I"${MKLROOT}/include" -L${MKLROOT}/lib/intel64 -lmkl_rt -Wl,--no-as-needed -lpthread -lm -ldl -o ProjecTracer
