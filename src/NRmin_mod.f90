MODULE NRmin_mod

! Module used by the NRmin program to perform a non-linear fitting
! by Newton-Raphson method. Also is used to estimate bias and error
! using toy MC.

USE mod_random

IMPLICIT NONE
CHARACTER (LEN=20) :: file_name ! data file
INTEGER, PARAMETER :: in_unit = 15
INTEGER :: err_open, err_read
REAL, ALLOCATABLE, DIMENSION(:) :: x, y, sigma ! data arrays
REAL, DIMENSION(2,1) :: par ! to keep parameters
REAL, DIMENSION(2,2) :: cov ! to save covariance matrix
REAL :: tol
INTEGER :: n_sim
REAL, DIMENSION(2,1):: bias, error

CONTAINS

    ! Subroutine to read data from file
    SUBROUTINE read_from_file(file_name, x, y, sigma)
        IMPLICIT NONE
        CHARACTER(LEN=20), INTENT(IN) :: file_name
        INTEGER :: data_count ! for  checking number of data points
        REAL :: data_value ! dummy float variable
        INTEGER :: i
        REAL, ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: x, y, sigma


        WRITE(*,*) "--- importing data ---"
        ! Open file
        WRITE(*,*) 'Opening file: ', file_name
        OPEN(in_unit, FILE = file_name, STATUS ='OLD', IOSTAT = err_open)
        IF (err_open > 0 ) STOP 'Error while opening file'

        ! Count number of entries in file
        data_count = 0
        DO
            READ(in_unit, *, IOSTAT = err_read) data_value
            IF (err_read > 0 ) THEN
                STOP 'Error while reading file'
            ELSE IF (err_read < 0) THEN
                EXIT
            ELSE
                data_count = data_count+1
            END IF
        END DO
        WRITE(*,*) '# data entries', data_count

        ! Read data from file
        REWIND(in_unit)
        ALLOCATE(x(data_count))
        ALLOCATE(y(data_count))
        ALLOCATE(sigma(data_count))
        DO i=1,data_count
            READ(in_unit, *, IOSTAT = err_read) x(i), y(i), sigma(i)
            IF(err_read > 0) STOP 'Error whilw reading the file'
        END DO

    END SUBROUTINE read_from_file

    ! Function obtain the chi2 value for the least square fitting
    FUNCTION chi2_func(x,y,sigma, par)
        IMPLICIT NONE
        REAL, ALLOCATABLE, DIMENSION(:), INTENT(IN) :: x, y, sigma
        REAL, DIMENSION(2,1), INTENT(IN) :: par
        REAL :: chi2_func

        chi2_func = SUM((-par(2,1)*x**2 - EXP(par(1,1)*x) + y )**2*sigma**(-2))

    END FUNCTION chi2_func

    ! Function to obtain the Hessian of the function Chi^2
    ! respect to parameters a and b, using analytical derivatives
    FUNCTION hessian_func(x, y, sigma, par)
        IMPLICIT NONE
        REAL, ALLOCATABLE, DIMENSION(:), INTENT(IN) :: x, y, sigma
        REAL, DIMENSION(2,1), INTENT(IN) :: par
        REAL, DIMENSION(2,2) :: hessian_func

        hessian_func(1,1) = SUM(2*sigma**(-2)*x**2 &
            & *((par(2,1)*x**2-y)*EXP(par(1,1)*x)+2*EXP(2*par(1,1)*x)))
        hessian_func(1,2) = SUM(2*sigma**(-2)*x**3*EXP(par(1,1)*x))
        hessian_func(2,2) = SUM(2*sigma**(-2)*x**4)
        hessian_func(2,1) = SUM(2*sigma**(-2)*x**3*EXP(par(1,1)*x))

    END FUNCTION hessian_func

    ! Function to obtain the gradient of the function Chi^2
    ! respect to parameters a and b, using analytical derivatives
    FUNCTION grad_func(x, y, sigma, par)
        IMPLICIT NONE
        REAL, ALLOCATABLE, DIMENSION(:), INTENT(IN) :: x, y, sigma
        REAL, DIMENSION(2,1), INTENT(IN) :: par
        REAL, DIMENSION(2,1) :: grad_func

        grad_func(1,1) = SUM(2*sigma**(-2)*x &
            & *((par(2,1)*x**2-y + EXP(par(1,1)*x))*EXP(par(1,1)*x)))
        grad_func(2,1) = SUM(2*sigma**(-2)*x**2*(par(2,1)*x**2+EXP(par(1,1)*x)-y))

    END FUNCTION grad_func

    ! Function to calculate the inverse of a 2x2 matrix
    FUNCTION inverse_func( mat )
        IMPLICIT NONE
        REAL, DIMENSION(2,2), INTENT(IN):: mat
        REAL :: det
        REAL, DIMENSION(2,2) :: inverse_func

        det = mat(1,1)*mat(2,2) - mat(1,2)*mat(2,1)
        inverse_func(1,1) = mat(2,2)/det
        inverse_func(1,2) = - mat(2,1)/det
        inverse_func(2,2) = mat(1,1)/det
        inverse_func(2,1) = - mat(1,2)/det

    END FUNCTION inverse_func

    ! Subroutine to perform a least squares estimation
    ! of the parameters a and b using Newton-Raphson method
    SUBROUTINE nr_fit(x, y, sigma, par, tol , verbose)
        IMPLICIT NONE
        REAL, ALLOCATABLE, DIMENSION(:), INTENT(IN) :: x, y, sigma
        REAL, DIMENSION(2,1), INTENT(INOUT) :: par ! also as output
        REAL, INTENT(IN) :: tol
        LOGICAL, INTENT(IN) :: verbose
        REAL, DIMENSION(2,2) :: H, H_inv
        REAL, DIMENSION(2,1) :: grad
        REAL, DIMENSION(2,1) :: n_par
        INTEGER :: iter
        REAL :: diff

        IF(verbose) WRITE(*,*) "--- minimization of \Chi^2 by Newton-Raphson ---"
        diff = HUGE(diff) ! be sure of pass condition
        iter = 0 ! init iterator

        ! Loop until difference less than tolerance or not converging
        DO WHILE ( diff > tol )
            IF(verbose) WRITE(*,"(A, I4, A,F8.5,A,F8.5,A,F8.3)") "step " &
                & , iter,  ": a = ", par(1,1), "  b = ", par(2,1) & 
                & , " \Chi^2 = ", chi2_func(x, y, sigma, par)
            iter = iter + 1
            H = hessian_func(x,y,sigma, par) ! obtain hessian
            H_inv = inverse_func(H) ! invert hessian
            grad =  grad_func(x, y, sigma, par) ! get gradient
            n_par = par - MATMUL(H_inv, grad) ! estimate new parameters
            diff = MAXVAL(ABS(par-n_par))
            par = n_par ! update
        END DO

        IF(verbose) WRITE(*,"(A, I4, A,F8.5,A,F8.5,A,F8.3)") "step " &
                & , iter,  ": a = ", par(1,1), "  b = ", par(2,1) & 
                & , " \Chi^2 = ", chi2_func(x, y, sigma, par)
        IF(verbose) WRITE(*,*) "convergence reached!"
    END SUBROUTINE nr_fit

    ! Function to obtain covariance matrix from Hessian
    FUNCTION cov_func(x, y, sigma, par)
        IMPLICIT NONE
        REAL, ALLOCATABLE, DIMENSION(:), INTENT(IN) :: x, y, sigma
        REAL, DIMENSION(2,1), INTENT(IN) :: par
        REAL, DIMENSION(2,2) :: hessian
        REAL, DIMENSION(2,2) :: cov_func

        WRITE(*,*) "--- estimating covariance matrix ---"

        hessian =  hessian_func(x, y, sigma, par)
        cov_func = inverse_func(hessian/2)

    END FUNCTION cov_func

    ! Subroutine to estimate bias and error from MC simulations
    SUBROUTINE simulate(x, y, sigma, par, n_sim, bias, error)
        IMPLICIT NONE
        REAL, ALLOCATABLE, DIMENSION(:), INTENT(IN) :: x, y, sigma
        REAL, DIMENSION(2,1), INTENT(IN) :: par
        INTEGER, INTENT(IN) :: n_sim
        REAL, DIMENSION(2,1), INTENT(OUT) :: bias, error
        REAL, ALLOCATABLE, DIMENSION(:) :: y_sim
        REAL, DIMENSION(2,1) :: t_par
        REAL, DIMENSION(2,n_sim) :: par_sim
        INTEGER :: idum = -347191 ! random seed
        INTEGER :: s, i

        ALLOCATE(y_sim(SIZE(y)))
        DO s=1,n_sim
            y_sim = EXP(par(1,1)*x) + par(2,1)*x**2
            DO i=1,SIZE(y)  ! add noise component
                y_sim(i) =  y_sim(i) + sigma(i)*gasdev(idum)
            END DO
            ! estimate parameters again
            t_par(1,1) = 1.
            t_par(2,1) = 1.
            CALL nr_fit(x, y_sim, sigma, t_par, 1e-5, .FALSE.)
            par_sim(1,s) = t_par(1,1)
            par_sim(2,s) = t_par(2,1)
            END DO

        bias(:,1) = SUM(par_sim,2)/n_sim - par(:,1)
        error(:,1) = SQRT(SUM(par_sim**2, 2)/n_sim &
            &- (SUM(par_sim,2)/n_sim)**2)

    END SUBROUTINE

END MODULE NRmin_mod

