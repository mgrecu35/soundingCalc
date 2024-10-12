! ############################################################
!            Satellite Data Simulation Unit v2
!                    -  df_mie_melt -
!
! This subroutine synthesizes the radiative properties within
! the melting layer. 
!
!                 Updated from df_mie-melt.f by H. M. (May, 2009)
!
 Subroutine df_mie_melt(freq,lwc,iwc,Nc_r,Nc_s,kext,asym,salb,pbck)
   Implicit NONE
!
! The number fraction of snow flakes to the rain and snow mixture 
! (=f_snow) is assumed to be linearly decreased downward from the 
! 0 C height. The number of sublayers is given by "nsublayer".
!
! I/O parameters
!
! freq : frequency in [GHz]
! lwc  : rain water content in [g/m3]
! iwc  : snow water content in [g/m3]
! Nc_r:  rain number concentration [1/m**3]
! Nc_s:  snow number concentration [1/m**3]
! kext : extinction coefficent in [1/km]
! salb : single scattering albedo
! asym : asymmetry parameter
! pbck : back-scattering phase function/(4*pi)
!
   real, intent(in)  :: freq,lwc,iwc,Nc_s,Nc_r
   real, intent(out) :: kext,asym,salb,pbck
!
! Local variables
!
   real :: kext_s,asym_s,salb_s,pbck_s,f_snow
   integer :: i
   integer, parameter :: nsublyr = 10
!
   kext = 0. ; asym = 0. ; salb = 0. ;  pbck = 0.
   if (lwc < 0.001 .and. iwc < 0.001) return
!
   Do i = 1 , nsublyr
      f_snow = float(i) / nsublyr - .5 / nsublyr
      Call mie_melt(freq, iwc, lwc, Nc_s, Nc_r, f_snow, &
                    kext_s, salb_s, asym_s, pbck_s)
      kext = kext + kext_s
      asym = asym + asym_s
      salb = salb + salb_s
      pbck = pbck + pbck_s
  End do
  kext = kext / nsublyr
  asym = asym / nsublyr
  salb = salb / nsublyr
  pbck = pbck / nsublyr
!
  Return
End Subroutine df_mie_melt
