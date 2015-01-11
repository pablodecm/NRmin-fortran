# NRmin-fortran

Newton-Raphson minimization for non-linear fitting.
Also includes bias and error estimation using simple toy MC.
Fortran project from programming course @ University of Cantabria


## Requirements

* gfortran (modify Makefile for other compilers)

## Compilation and execution

Set execution parameters (i.e. data file, initial values, tolerance or
numebr of simualtions) at *NRmod.f90*.

For compiling, a Makefile is provided in the src folder:

    cd src
    make

This will create a executable, that can be executed as:

    ./NRmin.bin

## IPython Notebook and Report

The report of this project is written as a IPython Notebook,
which can be checked online [here](http://nbviewer.ipython.org/github/pablodecm/NRmin-fortran/blob/master/report/NRmin-report.ipynb).

It can also be opened and executed ( if IPython, NumPy
and Matplotlib are available)by:
    cd report
    ipython notebook NRmin-report.ipynb

The report has been exported to Latex and pdf with the
nbconvert utility:
    cd report
    ipython nbconvert
