! ############################################################
!            Satellite Data Simulation Unit v2
!                      -  q2Wg -
!
! This function calculates rho*q in [g/m3] from mixing ratio,
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
!     q2Wg      R   W or rho*q [g/m3]
!
 Function q2Wg ( q,T,p,dflt )
   Implicit NONE
!
   Real, intent(in) :: q,T,p,dflt
   Real :: q2Wg
   Real :: rhod
   Real, parameter :: Rgas = 287.
!
   q2Wg = dflt
   If ( q == dflt ) return
   If ( T == dflt ) return
   If ( p == dflt ) return
!
   rhod = p * 100. / ( Rgas * T )
   q2Wg = rhod * q * 1.e3
!
 End Function q2Wg
