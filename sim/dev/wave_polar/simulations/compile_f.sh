#!/bin/bash

gfortran -O3 -o $1.o $1.f90 -I/usr/include/hdf5/serial -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_fortran -lhdf5
