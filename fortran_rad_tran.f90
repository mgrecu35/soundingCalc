subroutine calculate_tb_prof(iprof, incAngle, freqs, temperature_2m, u_10m, v_10m, qv, geopotential, &
      temperature, press, tb_all_freq) 
      ! Inputs
      integer, intent(in) :: iprof
      real, intent(in) :: incAngle
      real, intent(in), dimension(:) :: freqs
      real, intent(in), dimension(:) :: temperature_2m, u_10m, v_10m, press
      real, intent(in), dimension(:,:) :: qv, geopotential, temperature
  
      ! Outputs
      real, intent(out), dimension(size(freqs)) :: tb_all_freq
  
      ! Local variables
      integer :: nlyr, ifreq, k, nk, ireturn
      real :: fisot, umu, emis, ebar, w, btemp
      real, dimension(:), allocatable :: hL, lyrtemp, lyrhgt
      real, dimension(29) :: kext1D, asym1D, salb1D
      real :: rho, absair, abswv, tavg, ts
      integer :: nz(2)
      ! Constants
      fisot = 2.7
      umu = cos(incAngle / 180.0 * acos(-1.0))
      nlyr = 29
  
      ! Initialize output
      tb_all_freq = 0.0
  
      ! Loop over frequencies
      allocate(hL(0:nlyr))
      allocate(lyrtemp(0:nlyr))
      allocate(lyrhgt(0:nlyr))
      nz=shape(temperature)
      do ifreq = 1, size(freqs)
        ! Initialize profile-specific variables
        
  
        hL(0) = 0.0
        lyrtemp(0) = temperature_2m(iprof)
  
        do k = 1, nlyr
          nk = nz(2)+1-k
          ireturn = 0
          rho = press(nk) * 100.0 / (287.05 * temperature(iprof,nk))
          call gasabsr98(freqs(ifreq), temperature(iprof, nk), qv(iprof, nk) * rho, press(nk)*100, absair, abswv, ireturn)
          hL(k) = (geopotential(iprof, nk) + geopotential(iprof, nk-1)) / 2.0 / 9.8
          tavg = (temperature(iprof, nk) + temperature(iprof, nk-1)) / 2.0
          lyrtemp(k) = tavg
  
          kext1D(k) = absair + abswv
          asym1D(k) = 0.0
          salb1D(k) = 0.0
          !print*, kext1D(k), asym1D(k), salb1D(k), hL(k), lyrtemp(k), freqs(ifreq)
        end do
  
        ts = temperature_2m(iprof)
        w = sqrt(u_10m(iprof)**2 + v_10m(iprof)**2)
        call emit(freqs(ifreq), 1, ts, w, umu, emis, ebar)
  
        ! Convert heights to km
        
        lyrhgt = hL / 1.0e3
  
        btemp = temperature_2m(iprof)
        call radtran(umu, nlyr, tb_all_freq(ifreq), btemp, lyrtemp, lyrhgt, kext1D, salb1D, asym1D, fisot, emis, ebar, 0 )
  
        ! Deallocate local arrays
        
      end do
      deallocate(hL, lyrtemp, lyrhgt)
end subroutine calculate_tb_prof
  
subroutine calculate_tb_prof_emiss(iprof, incAngle, freqs, temperature_2m, u_10m, v_10m, &
    qv, geopotential, &
    temperature, press, emiss, press_at_surface, geopotential_at_surface,tb_all_freq) 
    ! Inputs
    integer, intent(in) :: iprof
    real, intent(in) :: incAngle
    real, intent(in), dimension(:) :: freqs
    real, intent(in), dimension(:) :: temperature_2m, u_10m, v_10m, press
    real, intent(in), dimension(:,:) :: qv, geopotential, temperature
    real, intent(in), dimension(:,:) :: emiss
    real, intent(in), dimension(:) :: press_at_surface, geopotential_at_surface

    ! Outputs
    real, intent(out), dimension(size(freqs)) :: tb_all_freq

    ! Local variables
    integer :: nlyr, ifreq, k, nk, ireturn
    real :: fisot, umu, emis, ebar, w, btemp
    real, dimension(:), allocatable :: hL, lyrtemp, lyrhgt
    real, dimension(29) :: kext1D, asym1D, salb1D
    real :: rho, absair, abswv, tavg, ts
    integer :: nz(2)
    integer :: isurface
    ! Constants
    fisot = 2.7
    umu = cos(incAngle / 180.0 * acos(-1.0))
    nlyr = 29

    ! Initialize output
    tb_all_freq = 0.0

    ! Loop over frequencies
    allocate(hL(0:nlyr))
    allocate(lyrtemp(0:nlyr))
    allocate(lyrhgt(0:nlyr))
    nz=shape(temperature)
    isurface = 1
    do while (press(nlyr+1-isurface) > press_at_surface(iprof)/1e2)
        isurface = isurface + 1
    end do
    do ifreq = 1, size(freqs)
      ! Initialize profile-specific variables
      

      hL(0) = 0.0
      lyrtemp(0) = temperature_2m(iprof)

      do k = 1, nlyr
        nk = nz(2)+1-k
        ireturn = 0
        rho = press(nk) * 100.0 / (287.05 * temperature(iprof,nk))
        call gasabsr98(freqs(ifreq), temperature(iprof, nk), qv(iprof, nk) * rho, press(nk)*100, absair, abswv, ireturn)
        hL(k) = (geopotential(iprof, nk) + geopotential(iprof, nk-1)) / 2.0 / 9.8
        tavg = (temperature(iprof, nk) + temperature(iprof, nk-1)) / 2.0
        lyrtemp(k) = tavg

        kext1D(k) = absair + abswv
        asym1D(k) = 0.0
        salb1D(k) = 0.0
        !print*, kext1D(k), asym1D(k), salb1D(k), hL(k), lyrtemp(k), freqs(ifreq)
      end do

      ts = temperature_2m(iprof)
      w = sqrt(u_10m(iprof)**2 + v_10m(iprof)**2)
      !call emit(freqs(ifreq), 1, ts, w, umu, emis, ebar)
      emis = emiss(iprof, ifreq)
      ebar = emis
      ! Convert heights to km
      
      lyrhgt = hL / 1.0e3

      btemp = temperature_2m(iprof)
      call radtran(umu, nlyr+1-isurface, tb_all_freq(ifreq), btemp, lyrtemp(isurface-1), lyrhgt(isurface-1), kext1D(isurface), salb1D(isurface), asym1D(isurface), fisot, emis, ebar, 0 )

      ! Deallocate local arrays
      
    end do
    deallocate(hL, lyrtemp, lyrhgt)
end subroutine calculate_tb_prof_emiss

subroutine calculate_tb_hyd(incAngle, freqs, qv, z, &
  temperature, press, rain_kext, rain_ksca, rain_asym, snow_kext, snow_ksca, snow_asym, &
  tb_all_freq,  nxny, nz, nfreqs) 
  ! Inputs
  integer, intent(in) :: nxny, nz, nfreqs
  real, intent(in) :: incAngle
  real, intent(in), dimension(nfreqs) :: freqs
  real, intent(in), dimension(nz,nxny) :: press
  real, intent(in), dimension(nz,nxny) :: qv, temperature
  real, intent(in), dimension(nz+1,nxny) :: z
  real, intent(in), dimension(nz,nxny,nfreqs) :: rain_kext, rain_ksca, rain_asym, snow_kext, snow_ksca, snow_asym

  ! Outputs
  real, intent(out), dimension(nxny,nfreqs) :: tb_all_freq

  ! Local variables
  integer :: nlyr, ifreq1, k, nk, ireturn
  real :: fisot, umu, emis, ebar, w, btemp
  real, dimension(:), allocatable :: hL, lyrtemp, lyrhgt
  real, dimension(nz) :: kext1D, asym1D, salb1D
  real :: rho, absair, abswv, tavg, ts, tb1
 
  ! Constants
  fisot = 2.7
  umu = cos(incAngle / 180.0 * acos(-1.0))
  nlyr = nz

  ! Initialize output
  tb_all_freq = 0.0

  ! Loop over frequencies
  allocate(hL(0:nlyr))
  allocate(lyrtemp(0:nlyr))
  allocate(lyrhgt(0:nlyr))
 
  do iprof=1,nxny
  
  do ifreq1 = 1, nfreqs
    ! Initialize profile-specific variables
    

    hL(0) = z(1,iprof)
    lyrtemp(0) = temperature(1,iprof)

    do k = 1, nlyr
      ireturn = 0
      rho = press(k,iprof) / (287.05 * temperature(k,iprof))
      call gasabsr98(freqs(ifreq1), temperature(k,iprof), qv(k,iprof) * rho, press(k,iprof), absair, abswv, ireturn)
      !print*, absair, abswv, k, temperature(k,iprof), qv(k,iprof), press(k,iprof), freqs(ifreq1)
      hL(k) = z(k+1,iprof)
      if(k.lt.nz) then
        tavg = (temperature(k,iprof) + temperature(k+1,iprof)) / 2.0
      else
        tavg = temperature(k,iprof)
      end if
      lyrtemp(k) = tavg

      kext1D(k) = absair + abswv
      asym1D(k) = 0.0
      salb1D(k) = 0.0
      
    end do

    ts = temperature(1,iprof)
    

    ! Convert heights to km
    
    lyrhgt = hL 
    emis=0.9
    ebar=0.9
    btemp = ts
    !print*, umu, nlyr, btemp
    !print*, lyrhgt
    !print*, lyrtemp
    !print*, kext1D
    !print*, salb1D
    !print*, asym1D
    !print*, fisot
    !print*, emis, ebar
    call radtran(umu, nlyr, tb1, btemp, lyrtemp, lyrhgt, kext1D, salb1D, asym1D, fisot, emis, ebar, 0 )
    tb_all_freq(iprof,ifreq1) = tb1
    end do
  end do
  deallocate(hL, lyrtemp, lyrhgt)
end subroutine calculate_tb_hyd