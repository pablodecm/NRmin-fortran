PROGRAM NRmin

! Fortran program to perform a non-linear fitting by least-square
! minimization using Newton-Raphson method. In addition, an estimation
! of the bias and uncertainty is carried out by toy MC.

! Most of the code has been implemented as subroutines of NRmin_mod
USE NRmin_mod
USE mod_chi2

IMPLICIT NONE

! read data from file
file_name = "../data/exp.dat"
CALL  read_from_file(file_name, x, y, sigma)

! estimate best parameters
tol = 1e-5
par(1,1) = 1. ! initial value of a
par(2,1) = 1. ! initial value of b
CALL  nr_fit(x, y, sigma, par, tol, .TRUE.)

! estimate covariance matrix
cov = cov_func(x, y, sigma, par)

WRITE(*,*) "--- NR minimization results --- "
WRITE(*,"(A,F7.5,A,F7.5)") "a parameter: ", par(1,1), " \pm ", SQRT(cov(1,1))
WRITE(*,"(A,F7.5,A,F7.5)") "b parameter: ", par(2,1), " \pm ", SQRT(cov(2,2))
WRITE(*,*) "--- goodness of fit ---"
WRITE(*,"(A,F5.2)") "Chi^2 value ", chi2_func(x, y, sigma, par)
WRITE(*,"(A,I3)") "Degrees of freedom: ", SIZE(x) - SIZE(par)
WRITE(*,"(A,F5.3)") "Probability chi^2(real)>chi^2(obtained): " &
    & , gammq(REAL(SIZE(x)-SIZE(par)), chi2_func(x,y,sigma,par))

WRITE(*,*) "--- MC simulation ---"
n_sim = 100000
WRITE(*,*) "# simulations: ", n_sim
CALL simulate(x, y, sigma, par, n_sim, bias, error )
WRITE(*,"(A,F8.5,A,F7.5)") "a parameter -> bias = ", bias(1,1) & 
    & ," error= ", error(1,1)
WRITE(*,"(A,F8.5,A,F7.5)") "b parameter -> bias = ", bias(2,1) & 
    & ," error= ", error(2,1)

WRITE(*,*) "--- final results with simulation --- "
WRITE(*,"(A,F7.5,A,F7.5)") "a parameter: ", par(1,1)-bias(1,1), " \pm ", error(1,1)
WRITE(*,"(A,F7.5,A,F7.5)") "b parameter: ", par(2,1)-bias(2,1), " \pm ", error(2,1)
WRITE(*,*) "--- end of the program ---"

END PROGRAM NRmin
