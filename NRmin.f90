PROGRAM NRmin

! Fortran program to perform a non-linear fitting by least-square
! minimization using Newton-Raphson method. In addition, an estimation
! of the bias and uncertainty is carried out by toy MC.

! Most of the code has been implemented as subroutines of NRmin_mod
USE NRmin_mod

IMPLICIT NONE

file_name = "exp.dat"
CALL  read_from_file(file_name, x, y, y_err)

END PROGRAM NRmin
