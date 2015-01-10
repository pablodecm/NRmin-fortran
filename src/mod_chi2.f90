module mod_chi2
! Funciones del Numerical Recipes

  contains
! Funcion para calcular la probabilidad de tener un valor mayor que
! "chi2" para una distribucion chi2 con "g" grados de libertad
  real function gammq(g, chi2)
    implicit none
    real, intent(in) :: g,chi2
    real :: a,x
    real :: gammcf,gamser, gln
    a=g/2.0
    x=chi2/2.0
    if ((x .lt. 0.) .or. (a .le. 0.)) stop 'bad arguments in gammq'
    if (x .lt. (a + 1.)) then
       call gser(gamser, a, x, gln)
       gammq = 1. - gamser
    else
       call gcf(gammcf, a, x, gln)
       gammq = gammcf
    end if
    return 
  end function gammq
  
  function gammln(xx)
    implicit none
    real :: gammln
    real, intent(in) :: xx
    integer :: j
    integer, parameter :: dp=kind(1d0)
    real(kind=dp), save :: stp = 2.5066282746310005d0 
    real(kind=dp), save, dimension(6) :: cof=(/ 76.18009172947146d0, &
    -86.50532032941677d0, 24.01409824083091d0, -1.231739572450155d0, &
    .1208650973866179d-2, -.5395239384953d-5 /)
    real(kind=dp) :: ser, tmp, x, y
    x = xx
    y = x
    tmp = x + 5.5d0
    tmp = ((x + 0.5d0) * log(tmp)) - tmp
    ser = 1.000000000190015d0
    do j = 1, 6
       y = y + 1.d0
       ser = ser + (cof(j) / y)
    end do
    gammln = tmp + log((stp * ser) / x)
    return 
  end function gammln

  subroutine gcf(gammcf, a, x, gln)
    implicit none
    integer, parameter :: itmax=100
    real, parameter :: eps = 3.e-7, fpmin = 1.e-30
    real :: a, gammcf, gln, x
    integer :: i
    real :: an, b, c, d, del, h
    gln = gammln(a)
    b = (x + 1.) - a
    c = 1. / fpmin
    d = 1. / b
    h = d
    do i = 1, itmax
       an = - (i * (i - a))
       b = b + 2.
       d = (an * d) + b
       if (abs(d) .lt. fpmin) d = fpmin
       c = b + (an / c)
       if (abs(c) .lt. fpmin) c = fpmin
       d = 1. / d
       del = d * c
       h = h * del
       if (abs(del - 1.) .lt. eps) goto 1
    end do
    stop 'a too large, ITMAX too small in gcf'
1   gammcf = exp(((- x) + (a * log(x))) - gln) * h
    return 
  end subroutine gcf
  
  subroutine gser(gamser, a, x, gln)
    implicit none
    integer :: n
    integer, parameter :: itmax=100
    real, parameter :: eps=3e-7
    real :: a, gamser, gln, x
    real :: ap, del, sum
    gln = gammln(a)
    if (x .le. 0.) then
       if (x .lt. 0.) stop 'x < 0 in gser'
       gamser = 0.
       return 
    end if
    ap = a
    sum = 1. / a
    del = sum
    do n = 1, itmax
       ap = ap + 1.
       del = (del * x) / ap
       sum = sum + del
       if (abs(del) .lt. (abs(sum) * eps)) goto 1
    end do
    stop 'a too large, ITMAX too small in gser'
1   gamser = sum * exp(((- x) + (a * log(x))) - gln)
    return 
  end subroutine gser
end module mod_chi2



