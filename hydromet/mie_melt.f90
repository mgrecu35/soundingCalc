! ############################################################
!            Satellite Data Simulation Unit v2
!                    - mie_melt -
!
! This subroutine computes the extinction, absorption, asymmetry 
! parameter and backscatter for a given water content of melting 
! particles in [g/m^3] synthesized from rain and snow.
!
!                 Updated from mie-melt.f by H. M. (May, 2009)
!
 Subroutine mie_melt(freq,iwc,lwc,Nc_s,Nc_r,f_snow,ksca,asca,gsca,pbck)
   Use MCRPparameters
   implicit none
!
! I/O parameters
!
! freq: frequency of radiation [GHz]
! iwc:  snow ice water content [g/m**3]
! lwc:  rain water content [g/m**3]
! Nc_s: snow number concentration [1/m**3]
! Nc_r: rain number concentration [1/m**3]
! f_snow: volume fraction of snow in a melting particle to synthesize
! ksca: extinction coefficient [1/km]
! asca: single-scatter albedo
! gsca: asymmetry facto
! pbck: backscatter phase function/(4*pi)
!
   real, intent(in)  :: freq,iwc,lwc,f_snow,Nc_s,Nc_r
   real, intent(out) :: ksca,asca,gsca,pbck
!
! Local variables
!
   integer :: i
   real :: wave,temp
   real :: densp,rad,D,num,num_s,num_r
   real :: qext,qsca,asym,qbsca
   real :: bext,bsca,bsym,bq11
   real :: faa,eicere,eiceim
   real :: norm_s,norm_r
!
   real, parameter :: dr = 5.e-2 ! Radius increment [mm] for PSD integration
   real, parameter :: rmax = 10. ! Maximum radius [mm] for PSD integration
   real, parameter :: dD = dr * 2
   integer :: imax
!
   real(8) :: xd,freqd,freqhzd,tempd,sald
   real(8) :: ewatred,ewatimd,eicered,eiceimd
   real(8) :: qscad,qextd,asymd,qbscad
!
   complex :: eice,ei,ewat
   complex(kind(0d0)) :: ewatd,eid,crefd,crefd_s,crefd_r
!
   imax = nint(rmax/dr)
   temp=273.16
   wave=300./freq
   if ( iwc < 0.001 .and. lwc < 0.001 ) then
      ksca=0. ; asca=0. ; gsca=0. ; pbck=0.
      return
   endif
!
   bext=0. ; bsca=0. ; bsym=0. ; bq11=0.
!
! Loop over particle sizes:
!
! normalize PSD if snow/rain%normaliz = .true.
!
   If ( snow%normaliz ) then
      norm_s=0.
      do i=0,imax
         rad=dr*(.5 + float(i))
         D = 2*rad
         num = psd_func ( snow,iwc+lwc,temp,Nc_s,D )
         densp = snow%paramset(1)
         norm_s = norm_s + num*densp* &                ! [1/m4]*[kg/m3]*
                           4*pi*((rad*1.e-3)**3)/3*dD  ! [m3]*[mm]=[g/m3]
      end do
      norm_s = norm_s / (iwc+lwc)
   Else
      norm_s = 1.0
   Endif
!
   If ( rain%normaliz ) then
      norm_r=0.
      do i=0,imax
         rad=dr*(.5 + float(i))
         D = 2*rad
         num = psd_func ( rain,iwc+lwc,temp,Nc_r,D )
         densp = rain%paramset(1)
         norm_r = norm_r + num*densp* &                ! [1/m4]*[kg/m3]*
                           4*pi*((rad*1.e-3)**3)/3*dD  ! [m3]*[mm]=[g/m3]
      end do
      norm_r = norm_r / (iwc+lwc)
   Else
      norm_r = 1.0
   Endif
!
   do i=0,imax
      rad=dr*(.5 + float(i))
      xd=dble(2.*pi*rad/wave)
      D = 2*rad
      num_s = psd_func ( snow,iwc+lwc,temp,Nc_s,D ) / norm_s
      num_r = psd_func ( rain,iwc+lwc,temp,Nc_r,D ) / norm_r
!
!       complex refractive index of snow
!
      freqd=dble(freq)
      tempd=dble(temp)
      call iceoptic(freqd,tempd,eicered,eiceimd)
      eicere=real(eicered)
      eiceim=real(eiceimd)
      eice=cmplx(eicere,eiceim)
!
!       calculate dielectric constant of snow as
!       ice matrix with air inclusions, using
!       Maxwell-Garnett mixing
!
      densp = snow%paramset(1)
      faa = 1.-(densp/densice)
      call mg_ellips(faa,eice,eair,ei)
      eid=dcmplx(ei)
      crefd_s=cdsqrt(eid)
!
!       complex refractive index of liquid water
!
      freqhzd=dble(freq*1.e+9)
      sald=0.d0
      call watoptic(freqhzd,tempd,sald,ewatred,ewatimd)
      ewatd=dcmplx(ewatred,ewatimd)
      crefd_r=cdsqrt(ewatd)
!
!       define melting particle properties
!
      num   = f_snow * num_s + (1. - f_snow) * num_r
!
! The old formula
!      crefd = .5d0 * (crefd_s + crefd_r)
! has been replacd with a Mawell-Garnett model
! with a water matrix and ice inclusions.
!
      eice=crefd_s*crefd_s
      ewat=cmplx(ewatred,ewatimd)
      faa = f_snow
      call mg_ellips(faa,ewat,eice,ei)
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
 end Subroutine mie_melt
