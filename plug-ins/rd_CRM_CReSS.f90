! ############################################################
!            Satellite Data Simulation Unit v2
!                      -  rd_CRM_CReSS -
!
! This subroutine reads cloud-resolving model profiles for the
! Cloud Resolving Storm Simulator (CReSS).
!
!                  ////// INSTRUCTIONS //////
!
! To read in CReSS simulations, you need to
!  1. Make sure you have your CReSS simulated data accessible
!     in "SDSU-v2/CRM/".
!  2. Edit this file so it properly interfaces your CReSS data.
!     The CReSS output format changes from a simulation to another, 
!     so this subroutine only offers a template and may not be 
!     applicable to all CReSS outputs as it is.
!  3. Replace "rd_CRM" with "rd_CRM_CReSS" in "SDSU_V2.f90". 
!  4. Edit "def_flie_list_CReSS.f90" so the input file list
!     matches your CReSS data files.
!  5. Replace "def_file_list" with "def_file_list_CReSS" in "SDSU_V2.f90". 
!  6. The SDSU output is named automatically from the input file
!     list. If your input file list includes subdirectory paths, 
!     the same subdirectories must exist in the output directory 
!     "SDSU-v2/outputs/".
!
!                   Updated from rd_CRM.f by H. M. (Aug, 2009)
!
 Subroutine rd_CRM_CReSS ( ifile )
   Use MAINparameters
   Implicit NONE
!
! I/O parameters
!
   Integer, intent(in) :: ifile
!
! Local varaibles
!
   Integer :: i,j,k,l,irec
   Real    :: wind_lev
   character :: input_file*90
   character :: pathandfile*90
!
! CReSS outputs
!
   Real, dimension(:,:), allocatable :: us    ! u wind at 10m 
   Real, dimension(:,:), allocatable :: vs    ! v wind at 10m
   Real, dimension(:,:), allocatable :: land, xlon, ylat
   Real, dimension(:,:,:), allocatable :: mean
   Real, dimension(:,:,:), allocatable :: pt ! base state potential temperature [K]
   Real, dimension(:,:,:), allocatable :: qv    ! water vapor mixing ratio [kg/kg]
   Real, dimension(:,:,:), allocatable :: qc    ! cloud water mixing ratio [kg/kg]
   Real, dimension(:,:,:), allocatable :: qr    ! rain water mixing ratio [kg/kg]
   Real, dimension(:,:,:), allocatable :: qi    ! cloud ice mixing ratio [kg/kg]
   Real, dimension(:,:,:), allocatable :: qs    ! snow mixing ratio [kg/kg]
   Real, dimension(:,:,:), allocatable :: qg    ! graupel mixing ratio [kg/kg]
   Real, dimension(:,:,:), allocatable :: nci   ! concentrations of cloud ice [1/kg]
   Real, dimension(:,:,:), allocatable :: ncs   ! concentrations of snow [1/kg]
   Real, dimension(:,:,:), allocatable :: ncg   ! concentrations of graupel [1/kg]
!
   Real :: pt2t, q2Wg, q2Wkg, qv2rh
   Real, parameter :: dflt = -1.e35
!
! CRM grid size
!  ngridx and ngridy: #s of grid points along x and y directions.
!  nlyr: # of vertical layers.
!  deltax and deltay: grid intervals along x and y directions.
!   (deltax and deltay will be used only for beam convolution)
!
! These values must be defined consistently with the CRM
! output dimensions. Check with the CRM experimental setups.
!
   ngridx = 160 ; ngridy = 160 ; nlyr = 45
   deltax = 4.0 ; deltay = 4.0
!      
! Define the input file name
!
   input_file = CRM_file_list(ifile)
   pathandfile = CRMdir//trim(input_file)
   open(unit = 14, file = pathandfile, access = 'DIRECT' , &
        recl = ngridx*ngridy*4, convert='BIG_ENDIAN' )
!
   print *, 'Input: '//trim(input_file)
!
! Allocate 2D variables
!
   allocate ( ics(ngridx,ngridy), tskin(ngridx,ngridy) )
   allocate ( rainrate_sfc(ngridx,ngridy), sfc_type(ngridx,ngridy) )
   allocate ( wind_sfc(ngridx,ngridy), soilH2O(ngridx,ngridy) )
   allocate ( albd_visir(ngridx,ngridy,mxwavel) )
!
! Allocate 3D variables defined at layer interfaces
!
   allocate (  hght_lev(ngridx,ngridy,0:nlyr) ) 
   allocate (  pres_lev(ngridx,ngridy,0:nlyr) ) 
   allocate (  temp_lev(ngridx,ngridy,0:nlyr) ) 
!
! Allocate 3D variables defined as layer averages
!
   allocate (  rlhm_avg(ngridx,ngridy,nlyr)   ) 
   allocate ( Wcldw_avg(ngridx,ngridy,nlyr)   ) 
   allocate ( Wrain_avg(ngridx,ngridy,nlyr)   ) 
   allocate ( Wcldi_avg(ngridx,ngridy,nlyr)   ) 
   allocate ( Wsnow_avg(ngridx,ngridy,nlyr)   ) 
   allocate ( Wgrpl_avg(ngridx,ngridy,nlyr)   ) 
   allocate ( Whail_avg(ngridx,ngridy,nlyr)   ) 
   allocate ( Ncldw_avg(ngridx,ngridy,nlyr)   ) 
   allocate ( Nrain_avg(ngridx,ngridy,nlyr)   ) 
   allocate ( Ncldi_avg(ngridx,ngridy,nlyr)   ) 
   allocate ( Nsnow_avg(ngridx,ngridy,nlyr)   ) 
   allocate ( Ngrpl_avg(ngridx,ngridy,nlyr)   ) 
   allocate ( Nhail_avg(ngridx,ngridy,nlyr)   ) 
!
! Allocate CReSS parameters
!
   allocate ( us   (ngridx,ngridy) )
   allocate ( vs   (ngridx,ngridy) )
   allocate ( xlon (ngridx,ngridy) )
   allocate ( ylat (ngridx,ngridy) )
   allocate ( land (ngridx,ngridy) )
   allocate ( mean (ngridx,ngridy,nlyr) )
   allocate ( qv   (ngridx,ngridy,nlyr) )
   allocate ( qc   (ngridx,ngridy,nlyr) )
   allocate ( qr   (ngridx,ngridy,nlyr) )
   allocate ( qi   (ngridx,ngridy,nlyr) )
   allocate ( qs   (ngridx,ngridy,nlyr) )
   allocate ( qg   (ngridx,ngridy,nlyr) )
   allocate ( nci  (ngridx,ngridy,nlyr) )
   allocate ( ncs  (ngridx,ngridy,nlyr) )
   allocate ( ncg  (ngridx,ngridy,nlyr) )
   allocate ( pt   (ngridx,ngridy,0:nlyr) )
!
! --------------- Read CRM profiles -----------------
!
! Vertical levels in [km]:
! These numbers must be defined consistently with the CRM
! output dimensions. Check with the CRM experimental setups.
!
   do i = 1, ngridx
      do j = 1, ngridy
         hght_lev(i,j,0:nlyr) = (/   0.000000E+00, &
         0.400000E-01, 0.141397E+00, 0.283215E+00, 0.460842E+00, 0.669948E+00, &
         0.906481E+00, 0.116667E+01, 0.144702E+01, 0.174433E+01, 0.205565E+01, &
         0.237833E+01, 0.271000E+01, 0.304856E+01, 0.339220E+01, 0.373937E+01, &
         0.408883E+01, 0.443959E+01, 0.479096E+01, 0.514251E+01, 0.549411E+01, &
         0.584590E+01, 0.619829E+01, 0.655199E+01, 0.690796E+01, 0.726748E+01, &
         0.763207E+01, 0.800356E+01, 0.838403E+01, 0.877586E+01, 0.918171E+01, &
         0.960451E+01, 0.100475E+02, 0.105141E+02, 0.110081E+02, 0.115337E+02, &
         0.120950E+02, 0.126968E+02, 0.133203E+02, 0.139443E+02, 0.145682E+02, &
         0.151922E+02, 0.158161E+02, 0.164401E+02, 0.170641E+02, 0.176880E+02 /)
      end do
   end do
!
! Surface parameters
!
   Read(14,rec=1)  ((us   (i,j),i=1,ngridx),j=1,ngridy) ! us
   Read(14,rec=2)  ((vs   (i,j),i=1,ngridx),j=1,ngridy) ! vs
   Read(14,rec=3)  ((pres_lev(i,j,0),i=1,ngridx),j=1,ngridy) ! ps
   Read(14,rec=4)  ((pt   (i,j,0),i=1,ngridx),j=1,ngridy) ! pts
   Read(14,rec=6)  ((tskin(i,j),i=1,ngridx),j=1,ngridy) ! tgs
!
! 3D parameters
!
   irec = 16
   Do l = 1, 3
      Do k = 1, nlyr
         irec = irec + 1
      End do
   End do
   Do k = 1, nlyr
      irec = irec + 1
      Read(14,rec=irec) ((mean    (i,j,k),i=1,ngridx),j=1,ngridy) ! pbar
   End do
   Do k = 1, nlyr
      irec = irec + 1
      Read(14,rec=irec) ((pres_lev(i,j,k),i=1,ngridx),j=1,ngridy) ! pp
      where ( mean(:,:,k) /= dflt .and. pres_lev(:,:,k) /= dflt )
         pres_lev(:,:,k) = mean(:,:,k) + pres_lev(:,:,k)
      elsewhere
         pres_lev(:,:,k) = dflt
      end where
   End do
   Do k = 1, nlyr
      irec = irec + 1
      Read(14,rec=irec) ((mean(i,j,k),i=1,ngridx),j=1,ngridy) ! ptbar
   End do
   Do k = 1, nlyr
      irec = irec + 1
      Read(14,rec=irec) ((pt  (i,j,k),i=1,ngridx),j=1,ngridy) ! ptp
      where ( mean(:,:,k) /= dflt .and. pt(:,:,k) /= dflt )
         pt(:,:,k) = mean(:,:,k) + pt(:,:,k)
      elsewhere
         pt(:,:,k) = dflt
      end where
   End do
   Do k = 1, nlyr
      irec = irec + 1
      Read(14,rec=irec) ((qv   (i,j,k),i=1,ngridx),j=1,ngridy) ! qv
   End do
   Do k = 1, nlyr
      irec = irec + 1
      Read(14,rec=irec) ((qc   (i,j,k),i=1,ngridx),j=1,ngridy) ! qc
   End do
   Do k = 1, nlyr
      irec = irec + 1
      Read(14,rec=irec) ((qr   (i,j,k),i=1,ngridx),j=1,ngridy) ! qr
   End do
   Do k = 1, nlyr
      irec = irec + 1
      Read(14,rec=irec) ((qi   (i,j,k),i=1,ngridx),j=1,ngridy) ! qi
   End do
   Do k = 1, nlyr
      irec = irec + 1
      Read(14,rec=irec) ((qs   (i,j,k),i=1,ngridx),j=1,ngridy) ! qs
   End do
   Do k = 1, nlyr
      irec = irec + 1
      Read(14,rec=irec) ((qg   (i,j,k),i=1,ngridx),j=1,ngridy) ! qg
   End do
   Do k = 1, nlyr
      irec = irec + 1
      Read(14,rec=irec) ((nci  (i,j,k),i=1,ngridx),j=1,ngridy) ! nci
   End do
   Do k = 1, nlyr
      irec = irec + 1
      Read(14,rec=irec) ((ncs  (i,j,k),i=1,ngridx),j=1,ngridy) ! ncs
   End do
   Do k = 1, nlyr
      irec = irec + 1
      Read(14,rec=irec) ((ncg  (i,j,k),i=1,ngridx),j=1,ngridy) ! ncg
   End do
!
! pressure in [hPa]
!
   where ( pres_lev /= dflt ) &
   pres_lev = pres_lev * 1.e-2 ! [Pa]->[hPa]
!
! Convective/stratiform separation (temporarily set to be all stratiform)
!
   ics = 2
!
   Do i = 1, ngridx
      Do j = 1, ngridy
!
         Do k = 0, nlyr
!
! temperature [K]
!
            temp_lev(i,j,k) = pt2t( pt(i,j,k),pres_lev(i,j,k),dflt )
!
         End do
!
         Do k = 1, nlyr
!
! relatvie humidity [%]
!
            rlhm_avg(i,j,k) = qv2rh( qv(i,j,k),temp_lev(i,j,k), &
                                     pres_lev(i,j,k),dflt )
!
! hydrometeor water contents in [g/m3]
! 
!
            Wcldw_avg(i,j,k) = q2Wg( qc(i,j,k),temp_lev(i,j,k), &
                                     pres_lev(i,j,k),dflt )
            Wrain_avg(i,j,k) = q2Wg( qr(i,j,k),temp_lev(i,j,k), &
                                     pres_lev(i,j,k),dflt )
            Wcldi_avg(i,j,k) = q2Wg( qi(i,j,k),temp_lev(i,j,k), &
                                     pres_lev(i,j,k),dflt )
            Wgrpl_avg(i,j,k) = q2Wg( qg(i,j,k),temp_lev(i,j,k), &
                                     pres_lev(i,j,k),dflt )
            Wsnow_avg(i,j,k) = q2Wg( qs(i,j,k),temp_lev(i,j,k), &
                                     pres_lev(i,j,k),dflt )
            Ncldi_avg(i,j,k) = q2Wkg( nci(i,j,k),temp_lev(i,j,k), &
                                      pres_lev(i,j,k),dflt )
            Ngrpl_avg(i,j,k) = q2Wkg( ncg(i,j,k),temp_lev(i,j,k), &
                                      pres_lev(i,j,k),dflt )
            Nsnow_avg(i,j,k) = q2Wkg( ncs(i,j,k),temp_lev(i,j,k), &
                                      pres_lev(i,j,k),dflt )
            Whail_avg(i,j,k) = 0.
            Ncldw_avg(i,j,k) = 1.e+8
            Nrain_avg(i,j,k) = 0.
            Nhail_avg(i,j,k) = 0.
!
         end do
!
       end do
    end do
!
! Replace default values with zeros for hydrometeor parameters
!
    Where ( Wcldw_avg == dflt ) Wcldw_avg = 0.
    Where ( Wrain_avg == dflt ) Wrain_avg = 0.
    Where ( Wcldi_avg == dflt ) Wcldi_avg = 0.
    Where ( Wsnow_avg == dflt ) Wsnow_avg = 0.
    Where ( Wgrpl_avg == dflt ) Wgrpl_avg = 0.
    Where ( Ncldi_avg == dflt ) Ncldi_avg = 0.
    Where ( Nsnow_avg == dflt ) Nsnow_avg = 0.
    Where ( Ngrpl_avg == dflt ) Ngrpl_avg = 0.
! 
! remove default values
!
    pres_lev(:,:,nlyr) = pres_lev(:,:,nlyr-1)
    temp_lev(:,:,nlyr) = temp_lev(:,:,nlyr-1)
    rlhm_avg(:,:,nlyr) = rlhm_avg(:,:,nlyr-1)
!
    Wcldw_avg(:,:,nlyr) = 0.
    Wrain_avg(:,:,nlyr) = 0. 
    Wcldi_avg(:,:,nlyr) = 0. ; Ncldi_avg(:,:,nlyr) = 0.
    Wgrpl_avg(:,:,nlyr) = 0. ; Ngrpl_avg(:,:,nlyr) = 0. 
    Wsnow_avg(:,:,nlyr) = 0. ; Nsnow_avg(:,:,nlyr) = 0.
!
    Do k = nlyr, 2, -1
!
! Fill in undefined values in some lower-most layers with the numbers
! from the layer directly above, except for relative humidity for which a
! very low dummy number is inserted. This is only a tentative remedy to
! cope with complicated topography.
!
       where ( pres_lev(:,:,k-1) == dflt ) pres_lev(:,:,k-1) = pres_lev(:,:,k)
       where ( temp_lev(:,:,k-1) == dflt ) temp_lev(:,:,k-1) = temp_lev(:,:,k)
       where ( rlhm_avg(:,:,k-1) == dflt ) rlhm_avg(:,:,k-1) = 1.e-3
!
! Apply interpolation (layer interfaces -> layer averages)
!
       where ( rlhm_avg(:,:,k-1) /= 1.e-3 ) &
       rlhm_avg (:,:,k) = 0.5 * (rlhm_avg (:,:,k) + rlhm_avg (:,:,k-1))
!
       Wcldw_avg(:,:,k) = 0.5 * (Wcldw_avg(:,:,k) + Wcldw_avg(:,:,k-1))
       Wrain_avg(:,:,k) = 0.5 * (Wrain_avg(:,:,k) + Wrain_avg(:,:,k-1))
       Wcldi_avg(:,:,k) = 0.5 * (Wcldi_avg(:,:,k) + Wcldi_avg(:,:,k-1))
       Wgrpl_avg(:,:,k) = 0.5 * (Wgrpl_avg(:,:,k) + Wgrpl_avg(:,:,k-1))
       Wsnow_avg(:,:,k) = 0.5 * (Wsnow_avg(:,:,k) + Wsnow_avg(:,:,k-1))
       Ncldi_avg(:,:,k) = 0.5 * (Ncldi_avg(:,:,k) + Ncldi_avg(:,:,k-1))
       Ngrpl_avg(:,:,k) = 0.5 * (Ngrpl_avg(:,:,k) + Ngrpl_avg(:,:,k-1))
       Nsnow_avg(:,:,k) = 0.5 * (Nsnow_avg(:,:,k) + Nsnow_avg(:,:,k-1))
    End do
!
! Define near-surface wind (applicable when background="water") and 
! soil moisture (when background="land").
!
    wind_sfc(:,:) = sqrt(us(:,:)**2 + vs(:,:)**2)
    soilH2O(:,:) = 0.
!
! Define visible/IR ground albedo (applicable when background="land")
!
    albd_visir = 0.2
!
! / modified for SDSUv2.1.0 /
!
! Surface type: defines local surface type when surface type varies spatially.
!               Activated only when background="mixed" (otherwise ignored).
!               Visit 'param_set.f90' to check the setting of 'background'.
!               (By default background='water'.)
!    sfc_type = 0 for open ocean
!             = 1 for land
!
    close ( 14 )
    do l = 1, 100
       if ( input_file(l:l) == '.' ) exit
    end do
    pathandfile = CRMdir//input_file(1:l)//'geography.united.bin'
    open(unit = 14, file = pathandfile, access = 'DIRECT' , &
         recl = ngridx*ngridy*4, convert='BIG_ENDIAN' )
!
    read(14,rec=2) ((ylat(i,j),i=1,ngridx),j=1,ngridy)
    read(14,rec=3) ((xlon(i,j),i=1,ngridx),j=1,ngridy)
    read(14,rec=6) ((land(i,j),i=1,ngridx),j=1,ngridy)
!
    where ( land < 0. )
       sfc_type = 0
    elsewhere
       sfc_type = 1
    endwhere
!
    close ( 14 )
   deallocate ( mean,us,vs,pt,land,xlon,ylat )
   deallocate ( qv,qc,qr,qi,qs,qg,nci,ncs,ncg )
    
    Return
  End Subroutine rd_CRM_CReSS
