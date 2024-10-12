! ############################################################
!            Satellite Data Simulation Unit v2
!                      -  pwvsat -
!
! Function to compute saturation WV pressure in [hPa] for 
! a given air temperature in [K]
!
 Real Function pwvsat ( T )
   Implicit NONE
!
   Real :: x,T
!
! Goff Gratch (1946) equation
!
   x = 373.16/T
   pwvsat = -7.90298 * (x - 1.) + 5.02808 * log10(x) -     &
             1.3816e-7 * (10.**(11.344*(1.-1./x)) - 1.) +  &
             8.1328e-3 * (10.**(-3.49149*(x-1.)) - 1.) + log10(1013.246)
   pwvsat = 10.**pwvsat
!
 End Function pwvsat
