#!/bin/bash

gfortran -c NRmin_mod.f90
gfortran -c mod_chi2.f90
gfortran -c NRmin.f90
gfortran NRmin.o NRmin_mod.o mod_chi2.o -o NRmin.bin
