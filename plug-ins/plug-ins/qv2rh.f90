! ############################################################
!            Satellite Data Simulation Unit v2
!                      -  qv2rh -
!
! This function calculates relative humidity from water vapor
! mixing ratio, temperature, and pressure.
!
!                                 Created by H. M. (Aug, 2009)
!
! Input:
!     qv        R   water vapor mixing ratio [kg/kg]
!     T         R   temperature [K]
!     p         R   pressure [hPa]
!     dflt      R   default value (for which the same value 
!                   will be returned)
!    
! Output:
!     qv2rh     R   relative humidity [%]
!
 Function qv2rh ( qv,T,p,dflt )
   Implicit NONE
!
   Real, intent(in) :: qv,T,p,dflt
   Real :: qv2rh,q2Wkg
   Real :: es,ewv,rhov,pwvsat
   Real, parameter :: Rvap = 462.
!
   qv2rh = dflt
   If ( qv == dflt ) return
   If ( T == dflt ) return
   If ( p == dflt ) return
!
! Saturation vapor pressure
!
   es = pwvsat ( T ) * 100.
!
! vapor pressure [Pa]
!
   rhov = q2Wkg ( qv,T,p,dflt )
   ewv = rhov * Rvap * T
!
! relative humidity
!
   qv2rh = min(ewv / es, 1.0) * 100.
!
 End Function qv2rh
