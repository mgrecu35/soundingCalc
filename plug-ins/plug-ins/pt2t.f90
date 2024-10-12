! ############################################################
!            Satellite Data Simulation Unit v2
!                      -  pt2t -
!
! This function calculates air temperature from potential 
! temperature and pressure.
!
!                                 Created by H. M. (Aug, 2009)
!
! Input:
!     pottemp   R   potential temperature [K]
!     p         R   pressure [hPa]
!     dflt      R   default value (for which the same value 
!                   will be returned)
! Output:
!     pt2t      R   temperature [K]
!
 Function pt2t ( pottemp, p, dflt )
   Implicit NONE
!
   Real, intent(in) :: pottemp, p, dflt
   Real :: pt2t
!
   pt2t = dflt 
   If ( pottemp == dflt ) return
   If ( p == dflt ) return
!
   pt2t = pottemp * ( p / 1000. )**0.2856
!
 End Function pt2t
