#!/bin/bash

gfortran -c NRmin_mod.f90
gfortran -c NRmin.f90
gfortran NRmin.o NRmin_mod.o -o NRmin.bin
