! ############################################################
!            Satellite Data Simulation Unit v2
!                      -  q2Wkg -
!
! This function calculates rho*q in [kg/m3] from mixing ratio,
! temperature, and pressure.
!
!                                 Created by H. M. (Aug, 2009)
!
! Input:
!     q         R   mixing ratio [kg/kg]
!     T         R   temperature [K]
!     p         R   pressure [hPa]
!     dflt      R   default value (for which the same value 
!                   will be returned)
!    
! Output:
!     q2Wkg     R   W or rho*q [kg/m3]
!
 Function q2Wkg ( q,T,p,dflt )
   Implicit NONE
!
   Real, intent(in) :: q,T,p,dflt
   Real :: q2Wkg
   Real :: rhod
   Real, parameter :: Rgas = 287.
!
   q2Wkg = dflt
   If ( q == dflt ) return
   If ( T == dflt ) return
   If ( p == dflt ) return
!
   rhod = p * 100. / ( Rgas * T )
   q2Wkg = rhod * q
!
 End Function q2Wkg
