! ############################################################
!            Satellite Data Simulation Unit v2
!                    - mie_hail -
!
! This subroutine computes the extinction, absorption, asymmetry 
! parameter and backscatter for a given water content of hail
! in [g/m^3]. Microphysical properties are specified in the
! data structure 'hail'.
!
!                 Updated from mie-hail_bs.f by H. M. (May, 2009)
!
 Subroutine mie_hail(freq,temp,iwc,Nc,ksca,asca,gsca,pbck)
   Use MCRPparameters
   implicit none
!
! I/O parameters
!
! freq: frequency of radiation [GHz]
! temp: temperature of particles [K]
! iwc:  hail ice water content [g/m**3]
! Nc:   hail number concentration [1/m**3]
! ksca: extinction coefficient [1/km]
! asca: single-scatter albedo
! gsca: asymmetry facto
! pbck: backscatter phase function/(4*pi)
!
   real, intent(in)  :: freq,temp,iwc,Nc
   real, intent(out) :: ksca,asca,gsca,pbck
!
! Local variables
!
   integer :: i
   real :: wave
   real :: densp,rad,D,num
   real :: qext,qsca,asym,qbsca
   real :: bext,bsca,bsym,bq11
   real :: faa,eicere,eiceim
   real :: norm
!
   real, parameter :: dr = 5.e-2 ! Radius increment [mm] for PSD integration
   real, parameter :: rmax = 10. ! Maximum radius [mm] for PSD integration
   real, parameter :: dD = dr * 2
   integer :: imax
!
   real(8) :: xd,freqd,tempd,sald
   real(8) :: ewatred,ewatimd,eicered,eiceimd
   real(8) :: qscad,qextd,asymd,qbscad
!
   complex :: eice,ei
   complex(kind(0d0)) :: ewatd,eid,crefd
!
   imax = nint(rmax/dr)
   wave=300./freq
   if (iwc < 0.001) then
      ksca=0. ; asca=0. ; gsca=0. ; pbck=0.
      return
   endif
!
   bext=0. ; bsca=0. ; bsym=0. ; bq11=0.
!
! Loop over particle sizes:
!
! normalize PSD if hail%normaliz = .true.
!
   If ( hail%normaliz ) then
      norm=0.
      do i=0,imax
         rad=dr*(.5 + float(i))
         D = 2*rad
         num = psd_func ( hail,iwc,temp,Nc,D )
         densp = hail%paramset(1)
         norm = norm + num*densp* &                ! [1/m4]*[kg/m3]*
                       4*pi*((rad*1.e-3)**3)/3*dD  ! [m3]*[mm]=[g/m3]
      end do
      norm = norm / iwc
   Else
      norm = 1.0
   Endif
!
   do i=0,imax
      rad=dr*(.5 + float(i))
      xd=dble(2.*pi*rad/wave)
      D = 2*rad
      num = psd_func ( hail,iwc,temp,Nc,D ) / norm
!
!       complex refractive index of hail
!
      freqd=dble(freq)
      tempd=dble(temp)
      call iceoptic(freqd,tempd,eicered,eiceimd)
      eicere=real(eicered)
      eiceim=real(eiceimd)
      eice=cmplx(eicere,eiceim)
!
!       calculate dielectric constant of hail as
!       ice matrix with air inclusions, using
!       Maxwell-Garnett mixing
!
      densp = hail%paramset(1)
      faa = 1.-(densp/densice)
      call mg_ellips(faa,eice,eair,ei)
      eid=dcmplx(ei)
      crefd=cdsqrt(eid)
!
!       call Mie program
! 
      call mie_sphere(xd,crefd,qscad,qextd,asymd,qbscad)
!
      qext=real(qextd)
      qsca=real(qscad)
      asym=real(asymd)
      qbsca=real(qbscad)
!     
!       integrate over particle size distribution
!
      bext=bext+num*qext     *pi*rad*rad*dD
      bsca=bsca+num*qsca     *pi*rad*rad*dD
      bsym=bsym+num*qsca*asym*pi*rad*rad*dD
      bq11=bq11+num*qbsca    *pi*rad*rad*dD
!
   end do
!
!     check for distribution with very small extinction;
!     set parameters to zero to avoid numerical problems
!
    if ( bext > 1.e-6 ) then
!
      ksca=bext
      asca=bsca/bext
      gsca=bsym/bsca
      pbck=bq11/bsca
!
   else
!
      ksca=0.
      asca=0.
      gsca=0.
      pbck=0.
!
   end if
!
! Unit conversion [mm3/m-4] -> [km-1]
! 
   ksca = ksca  * 1.e-6
!
   return
 end Subroutine mie_hail
