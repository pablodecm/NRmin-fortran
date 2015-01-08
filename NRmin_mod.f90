MODULE NRmin_mod

! Module used by the NRmin program to perform a non-linear fitting
! by Newton-Raphson method. Also is used to estimate bias and error
! using toy MC.

IMPLICIT NONE
CHARACTER (LEN=20) :: file_name ! data file
INTEGER, PARAMETER :: in_unit = 15
INTEGER :: err_open, err_read
REAL, ALLOCATABLE, DIMENSION(:) :: x, y, sigma ! data arrays
REAL :: a,b ! fit parameters
REAL :: chi2 ! value to minimize

CONTAINS

    ! Subroutine to read data from file
    SUBROUTINE read_from_file(file_name, x, y, sigma)
        IMPLICIT NONE
        CHARACTER(LEN=20), INTENT(IN) :: file_name
        INTEGER :: data_count ! for  checking number of data points
        REAL :: data_value ! dummy float variable
        INTEGER :: i
        REAL, ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: x, y, sigma


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
                WRITE(*,*) 'Data file scanned'
                EXIT
            ELSE
                data_count = data_count+1
            END IF
        END DO
        WRITE(*,*) 'The number of entries in the file is ', data_count

        ! Read data from file
        REWIND(in_unit)
        ALLOCATE(x(data_count))
        ALLOCATE(y(data_count))
        ALLOCATE(sigma(data_count))
        DO i=1,data_count
            READ(in_unit, *, IOSTAT = err_read) x(i), y(i), sigma(i)
            IF(err_read > 0) STOP 'Error whilw reading the file'
        END DO
        WRITE(*,*) 'File read'

    END SUBROUTINE read_from_file

    ! Function obtain the chi2 value for the least square fitting
    FUNCTION chi2_func(x,y,sigma,a,b)
        IMPLICIT NONE
        REAL, ALLOCATABLE, DIMENSION(:), INTENT(IN) :: x, y, sigma
        REAL, INTENT(IN) :: a, b
        REAL :: chi2_func

        chi2_func = SUM((-b*x**2 - EXP(a*x) + y )**2*sigma**(-2))

    END FUNCTION chi2_func

    ! Function to obtain the Hessian of the function Chi^2
    ! respect to parameters a and b, using analytical derivatives
    FUNCTION hessian_func(x, y, sigma, par)
        IMPLICIT NONE
        REAL, ALLOCATABLE, DIMENSION(:), INTENT(IN) :: x, y, sigma
        REAL, DIMENSION(2), INTENT(IN) :: par
        REAL, DIMENSION(2,2) :: hessian_func

        hessian_func(1,1) = SUM(2*sigma**(-2)*x**2 &
            & *((b*x**2-y)*EXP(par(1)*x)+2*EXP(2*par(1)*x)))
        hessian_func(1,2) = SUM(2*sigma**(-2)*x**3*EXP(par(1)*x))
        hessian_func(2,2) = SUM(2*sigma**(-2)*x**2)
        hessian_func(2,1) = SUM(2*sigma**(-2)*x**3*EXP(par(1)*x))

    END FUNCTION hessian_func

    ! Function to obtain the gradient of the function Chi^2
    ! respect to parameters a and b, using analytical derivatives
    FUNCTION grad_func(x, y, sigma, par)
        IMPLICIT NONE
        REAL, ALLOCATABLE, DIMENSION(:), INTENT(IN) :: x, y, sigma
        REAL, DIMENSION(2), INTENT(IN) :: par
        REAL, DIMENSION(2) :: grad_func

        grad_func(1) = SUM(2*sigma**(-2)*x &
            & *((b*x**2-y)*EXP(par(1)*x)+2*EXP(2*par(1)*x)))
        grad_func(2) = SUM(2*sigma**(-2)*x**2*(par(2)*x**2+EXP(par(1)*x)-y))

    END FUNCTION grad_func

    ! Function to calculate the inverse of a 2x2 matrix
    FUNCTION inverse_func( mat )
        IMPLICIT NONE
        REAL, DIMENSION(2,2), INTENT(IN):: mat
        REAL :: det
        REAL, DIMENSION(2,2) :: inverse_func

        det = mat(1,1)*mat(2,2) - mat(1,2)*mat(2,1)
        inverse_func(1,1) = mat(2,2)
        inverse_func(1,2) = - mat(2,1)
        inverse_func(2,2) = mat(1,2)
        inverse_func(2,1) = mat(1,1)
        inverse_func = inverse_func/det

    END FUNCTION inverse_func

END MODULE NRmin_mod

