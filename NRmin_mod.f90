MODULE NRmin_mod

! Module used by the NRmin program to perform a non-linear fitting
! by Newton-Raphson method. Also is used to estimate bias and error
! using toy MC.

IMPLICIT NONE
CHARACTER (LEN=20) :: file_name
INTEGER, PARAMETER :: in_unit = 15
INTEGER :: err_open, err_read
REAL, ALLOCATABLE, DIMENSION(:) :: x, y, y_err

CONTAINS

    ! Subroutine to read data from file
    SUBROUTINE read_from_file(file_name, x, y, y_err)
        IMPLICIT NONE
        CHARACTER(LEN=20), INTENT(IN) :: file_name
        INTEGER :: data_count ! for  checking number of data points
        REAL :: data_value ! dummy float variable
        INTEGER :: i
        REAL, ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: x, y, y_err


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
        ALLOCATE(y_err(data_count))
        DO i=1,data_count
            READ(in_unit, *, IOSTAT = err_read) x(i), y(i), y_err(i)
            IF(err_read > 0) STOP 'Error whilw reading the file'
        END DO
        WRITE(*,*) 'File read'

    END SUBROUTINE read_from_file

END MODULE NRmin_mod

