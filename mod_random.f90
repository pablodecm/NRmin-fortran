module mod_random

! Funciones del Numerical Recipes

contains


    real function ran1(idum)

      integer :: idum, ia, im, iq, ir, ntab, ndiv
      real :: am, eps, rnmx
      parameter (ia = 16807, im = 2147483647, am = 1. / im, iq = 127773 &
           &, ir = 2836, ntab = 32, ndiv = 1 + ((im - 1) / ntab), eps = 1.2e-7 &
           &, rnmx = 1. - eps)
      integer j, k, iv(ntab), iy
      save iy, iv
      data iv / 32*0 /
      data iy / 0 /
      
      if ((idum .le. 0) .or. (iy .eq. 0)) then
         idum = max(- idum,1)
         do 11 j = ntab + 8, 1, -1
            k = idum / iq
            idum = (ia * (idum - (k * iq))) - (ir * k)
            if (idum .lt. 0) idum = idum + im
            if (j .le. ntab) iv(j) = idum
11       continue
            iy = iv(1)
      end if
      k = idum / iq
      idum = (ia * (idum - (k * iq))) - (ir * k)
      if (idum .lt. 0) idum = idum + im
      j = 1 + (iy / ndiv)
      iy = iv(j)
      iv(j) = idum
      ran1 = min(am * iy,rnmx)
      return 

   end function ran1


  real function gasdev(idum)

    integer idum
!    real, external :: ran1
    integer iset
    real fac, gset, rsq, v1, v2
    save gset, iset
    data iset / 0 /

    if (iset .eq. 0) then
1      v1 = (2. * ran1(idum)) - 1.
       v2 = (2. * ran1(idum)) - 1.
       rsq = (v1 ** 2) + (v2 ** 2)
       if ((rsq .ge. 1.) .or. (rsq .eq. 0.)) goto 1
       fac = sqrt(- ((2. * log(rsq)) / rsq))
       gset = v1 * fac
       gasdev = v2 * fac
       iset = 1
    else
       gasdev = gset
       iset = 0
    end if
    return 

    end function gasdev




end module mod_random

