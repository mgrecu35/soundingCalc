! ############################################################
!            Satellite Data Simulation Unit v2
!                    - mie_hail_vir -
!
! This subroutine computes the extinction and absorption
! coefficients and the phase function for a given water content 
! of hail in [g/m^3]. Microphysical properties are 
! specified in the data structure 'hail'. For use by the 
! visible/IR simulator.
!
!                                 created by H. M. (May, 2009)
!
 Subroutine mie_hail_vir(wl,temp,iwc,Nc,cext,csca,dscs)
   Use VIRparameters, only: knang,nang
   Use MCRPparameters
   Implicit none
!
! I/O parameters
!
! wl: wavelentgh [micron]
! temp: temperature of particles [K]
! iwc:  hail water content [g/m**3]
! Nc:   hail water number concentration [1/m**3]
! cext: extinction cross section [1/km]
! csca: scattering cross section [1/km]
! dscs: differential scattering cross section [1/km/str]
!       (= Csca * phase function)
!
   real(8), intent(in)  :: wl,temp,iwc,Nc
   real(8), intent(out) :: cext,csca
   real(8), dimension(knang), intent(out) :: dscs
!
! Local variables
!
   integer :: i
   real :: iwcs,temps,Ncs,Ds ! single-precision parameters to 
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
   real :: faa
   real(8) :: xd
   real(8) :: qscad,qextd
   real(8), dimension(knang) :: phsfd
!
   complex :: eice,ei
   complex(kind(0d0)) :: crefd,eid
!
   imax = nint(rmax/dr)
   wave = wl * 1.d-3
   if (iwc < 1.d-3) then
      cext=0.d0 ;  csca=0.d0 ; dscs=0.d0
      return
   endif
!
   bext=0.d0 ; bsca=0.d0 ; bdcs=0.d0
!
! Loop over particle sizes:
!
! normalize PSD if hail%normaliz = .true.
!
   If ( hail%normaliz ) then
      norm=0.d0
      do i=0,imax
         rad=dr*(.5d0 + 1.d0*i)
         D = 2*rad
         iwcs = iwc ; temps = temp ; Ncs = Nc ; Ds = D
         num = psd_func ( hail,iwcs,temps,Ncs,Ds )
         densp = hail%paramset(1)
         norm = norm + num*densp* &                 ! [1/m4]*[kg/m3]*
                       4*pid*((rad*1.d-3)**3)/3*dD  ! [m3]*[mm]=[g/m3]
      end do
      norm = norm / iwc
   Else
      norm = 1.0d0
   Endif
!
   do i=0,imax
      rad=dr*(.5d0 + 1.d0*i)
      xd=dble(2.*pid*rad/wave)
      D = 2*rad
      iwcs = iwc ; temps = temp ; Ncs = Nc ; Ds = D
      num = psd_func ( hail,iwcs,temps,Ncs,Ds ) / norm
!
!     complex refractive index of liquid water
!
      call icerefind_vir ( wl,crefd )
!
!     Maxwell-Garnett mixing
!
      eice=crefd*crefd
      densp = hail%paramset(1)
      faa = 1.-(densp/densice)
      call mg_ellips(faa,eice,eair,ei)
      eid=dcmplx(ei)
      crefd=cdsqrt(eid)
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
 End Subroutine mie_hail_vir
