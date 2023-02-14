#!/bin/bash

#######set up a job on a machine that does not require a queue

####prep for animation files
mkdir -p movie_files
rm -r movie_files
mkdir  movie_files

export OMP_NUM_THREADS=2

nohup ./ProjecTracer > stdout &
