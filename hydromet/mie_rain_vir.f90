! ############################################################
!            Satellite Data Simulation Unit v2
!                    - mie_rain_vir -
!
! This subroutine computes the extinction and absorption
! coefficients and the phase function for a given water content 
! of rain in [g/m^3]. Microphysical properties are 
! specified in the data structure 'rain'. For use by the 
! visible/IR simulator.
!
!                                 created by H. M. (May, 2009)
!
 Subroutine mie_rain_vir(wl,temp,lwc,Nc,cext,csca,dscs)
   Use VIRparameters, only: knang,nang
   Use MCRPparameters
   Implicit none
!
! I/O parameters
!
! wl: wavelentgh [micron]
! temp: temperature of particles [K]
! lwc:  rain water content [g/m**3]
! Nc:   rain water number concentration [1/m**3]
! cext: extinction cross section [1/km]
! csca: scattering cross section [1/km]
! dscs: differential scattering cross section [1/km/str]
!       (= Csca * phase function)
!
   real(8), intent(in)  :: wl,temp,lwc,Nc
   real(8), intent(out) :: cext,csca
   real(8), dimension(knang), intent(out) :: dscs
!
! Local variables
!
   integer :: i
   real :: lwcs,temps,Ncs,Ds ! single-precision parameters to 
                             ! interface with psd_func.
   real(8) :: wave
   real(8) :: rad,D,num
   real(8) :: qext,qsca
   real(8) :: bext,bsca
   real(8) :: densp,norm
   real(8), dimension(knang) :: pfnc,bdcs
!
   real(8), parameter :: dr = 2.d-2 ! Radius increment [mm] for PSD integration
   real(8), parameter :: rmax = 1.d+1 ! Maximum radius [mm] for PSD integration
   real(8), parameter :: dD = dr * 2
   real(8), parameter :: pid = 3.141592653589793d0
   integer :: imax
!
! Mie parameters
!
   real(8) :: xd
   real(8) :: ewatred,ewatimd
   real(8) :: qscad,qextd
   real(8), dimension(knang) :: phsfd
!
   complex(kind(0d0)) :: ewatd,crefd
!
   imax = nint(rmax/dr)
   wave = wl * 1.d-3
   if (lwc < 1.d-3) then
      cext=0.d0 ;  csca=0.d0 ; dscs=0.d0
      return
   endif
!
   bext=0.d0 ; bsca=0.d0 ; bdcs=0.d0
!
! Loop over particle sizes:
!
! normalize PSD if rain%normaliz = .true.
!
   If ( rain%normaliz ) then
      norm=0.d0
      do i=0,imax
         rad=dr*(.5d0 + 1.d0*i)
         D = 2*rad
         lwcs = lwc ; temps = temp ; Ncs = Nc ; Ds = D
         num = psd_func ( rain,lwcs,temps,Ncs,Ds )
         densp = rain%paramset(1)
         norm = norm + num*densp* &                 ! [1/m4]*[kg/m3]*
                       4*pid*((rad*1.d-3)**3)/3*dD  ! [m3]*[mm]=[g/m3]
      end do
      norm = norm / lwc
   Else
      norm = 1.0d0
   Endif
!
   do i=0,imax
      rad=dr*(.5d0 + 1.d0*i)
      xd=dble(2.*pid*rad/wave)
      D = 2*rad
      lwcs = lwc ; temps = temp ; Ncs = Nc ; Ds = D
      num = psd_func ( rain,lwcs,temps,Ncs,Ds ) / norm
!
!     complex refractive index of liquid water
!
      call watrefind_vir ( wl,crefd )
!
!     call Mie program
! 
      call mie_sphere_vir( xd,crefd,qscad,qextd,phsfd )
!
      qext=qextd
      qsca=qscad
      pfnc(1:nang)=phsfd(1:nang)
!
!       integrate over particle size distribution
!
      bext        =bext        +num*qext        *pid*rad*rad*dD
      bsca        =bsca        +num*qsca        *pid*rad*rad*dD
      bdcs(1:nang)=bdcs(1:nang)+num*pfnc(1:nang)*pid*rad*rad*dD*qsca
   end do
!
!     check for distribution with very small extinction;
!     set parameters to zero to avoid numerical problems
!
   if ( bext > 1.e-6 ) then
!
      cext=bext
      csca=bsca
      dscs(1:nang)=bdcs(1:nang)
!
   else
!
      cext=0.d0
      csca=0.d0
      dscs=0.d0
!
   end if
!
! Unit conversion [mm3/m-4] -> [km-1]
! 
   cext = cext * 1.d-6
   csca = csca * 1.d-6
   dscs = dscs * 1.d-6
! 
   return
 End Subroutine mie_rain_vir
