PROGRAM NRmin

! Fortran program to perform a non-linear fitting by least-square
! minimization using Newton-Raphson method. In addition, an estimation
! of the bias and uncertainty is carried out by toy MC.

! Most of the code has been implemented as subroutines of NRmin_mod
USE NRmin_mod

IMPLICIT NONE

file_name = "exp.dat"
CALL  read_from_file(file_name, x, y, sigma)

tol = 1e-5
par(1,1) = 1. ! initial value of a
par(2,1) = 1. ! initial value of b
CALL  nr_fit(x, y, sigma, par, tol)

WRITE(*,*) "BEST PARAMETERS: ", par

END PROGRAM NRmin
