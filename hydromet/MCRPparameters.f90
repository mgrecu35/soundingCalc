! ###########################################################
!            Satellite Data Simulation Unit v2
!                   - MCRPparameters -
!
! This module defines microphysical variables. 
!
!                                           H. M. (May, 2009)
!
 Module MCRPparameters
   Implicit NONE
!
! Particle size distributions
!
!  npmax: max number of PSD parameters
!  name: PSD name used for labeling the LUT files
!  functype: PSD function type (see psd_func below)
!  nparams: number of PSD parameters. must be <= npmax
!  paramset: PSD parameter set (see psd_func below)
!  normaliz: =.true. if external normalization is required.
!           (i.e., if m(D)psd(D) cannot be analytically integrated.)
!  indpndtT: =.true. if PSD parameters are independent of temperature.
!
   integer, parameter :: npmax = 10
   type psd_params
      character*8 :: name
      character*6 :: functype
      integer :: nparams
      real, dimension(npmax) :: paramset
      logical :: normaliz, indpndtT, twomomnt
   end type psd_params
   type(psd_params) :: rain,snow,cldw,cldi,grpl,hail
!
! mathematical and physical constants
!
   real, parameter :: pi = 3.14159265
   real, parameter :: densliq = 1.0e+3    ! liquid water density [kg/m3]
   real, parameter :: densice = 0.917e+3  ! frozen water density [kg/m3]
   complex :: eair  ! air refractive index, defined in microp_set.f90
   real :: normLN ! a constant defined in microp_set.f90
!
   Contains
!
! Particle size distributions
! - Rain DSD
!
     Real function psd_func ( hydr,lwc,T,Nc,D )
       Implicit NONE
!
! I/O parameters
!
!  hydr: microphysical parameter set
!  lwc: liquid/ice water content in [g/m^3]
!  T:   temperature in [K] (used for 'HyPl84')
!  Nc:  particle number concentration in [1/m^3] (used for 2-moment options)
!  D:   diameter in [mm]
!
       real, intent(in) :: lwc,T,Nc,D
       type(psd_params), intent(inout) :: hydr
!
! Local variables
!
       real :: density,num
       real :: n0,lam
       real :: Dmode,sigma
       real :: nw,mu,D0
       real :: a,b
       real, external :: gammaf
       integer :: nopt
!
       If ( hydr%nparams > npmax ) then
          print *, 'ERROR (MCRPparameters): npmax must be >= nparamx.'
          stop
       Endif
!      
       Select case ( hydr%functype )
!
! ** Exponential distribution with specified particle density
! ** and intercept parameter. 
!
          case( 'expntl' ) ;&
          hydr%normaliz = .false. ;&
          hydr%indpndtT = .true. ;&
          hydr%twomomnt = .false. ;&
          density  = hydr%paramset(1) ;&
          n0       = hydr%paramset(2) ;&
          lam      = (n0*pi*density/(lwc*1.e-3))**0.25 * 1.e-3 ;& ! [1/mm]
          psd_func = n0 * exp(-lam*D)
!
! ** Gamma disribution with specified particle density, mu,
! ** and normalized intercept parameter. 
!
          case( 'gammwN' ) ;&
          hydr%normaliz = .false. ;&
          hydr%indpndtT = .true. ;&
          hydr%twomomnt = .false. ;&
          density  = hydr%paramset(1) ;&
          nw       = hydr%paramset(2) ;&
          mu       = hydr%paramset(3) ;&
          D0       = 3.67*(lwc*1.e-3/(pi*density*nw))**0.25 * 1.e3  ;& ! [mm]
          psd_func = nw * 6*(3.67+mu)**(mu+4.)/(3.67**4*gammaf(mu+4.)) * &
                     (D/D0)**mu * exp(-(3.67+mu)*D/D0)
!
! ** Gamma disribution with specified particle density, mu,
! ** and median volume diameter.
!
          case( 'gammwD' ) ;&
          hydr%normaliz = .false. ;&
          hydr%indpndtT = .true. ;&
          hydr%twomomnt = .false. ;&
          density  = hydr%paramset(1) ;&
          D0       = hydr%paramset(2) ;&
          mu       = hydr%paramset(3) ;&
          nw       = lwc*1.e-3/(pi*density)*(3.67/(D0*1.e-3))**4  ;&
          psd_func = nw * 6*(3.67+mu)**(mu+4.)/(3.67**4*gammaf(mu+4.)) * &
                     (D/D0)**mu * exp(-(3.67+mu)*D/D0)
!
! ** Log-normal disbribution with specified particle density,
! ** mode radius, and dispersion. 
!
          case( 'lgnorm' ) ;&
          hydr%normaliz = .false. ;&
          hydr%indpndtT = .true. ;&
          hydr%twomomnt = .false. ;&
          density = hydr%paramset(1) ;& 
          Dmode   = hydr%paramset(2) ;&
          sigma   = hydr%paramset(3) ;&
          psd_func = normLN * lwc*1.e-3 / (density*sigma*(Dmode*1.e-3)**3) * &
                     exp( -4.5*sigma*sigma ) * &
                     exp( -0.5*((log(D)-log(Dmode))/sigma)**2 ) / &
                     (D*1.e-3)
!
! ** Exponential PSD with a power-law mass spectrum m(D) with specified
! ** intercept parameter, 'a', and 'b', where m(D) = a*D^b
!     lam=[a*N0*Gamma(b+1)/lwc]^{1/(b+1)}
!
          case( 'expaDb' ) ;&
          hydr%normaliz = .false. ;&
          hydr%indpndtT = .true. ;&
          hydr%twomomnt = .false. ;&
          n0       = hydr%paramset(2) ;&
          a        = hydr%paramset(3) ;&
          b        = hydr%paramset(4) ;&
          density  = 6./pi * a * (D*1.e-3)**(b-3.) ;&
          hydr%paramset(1) = density ;&
          lam      = (a*n0*gammaf(b+1.)/(lwc*1.e-3))**(1./(b+1.)) &
                      * 1.e-3 ;& ! [1/mm]
          psd_func = n0 * exp(-lam*D)
!
! ** Modified gamma distribution, suitable for a liquid cloud DSD, 
! ** with specified particle density and mode radius.
!
          case( 'modGam' ) ;&
          hydr%normaliz = .false. ;&
          hydr%indpndtT = .true. ;&
          hydr%twomomnt = .false. ;&
          density = hydr%paramset(1) ;&
          Dmode   = hydr%paramset(2) ;&
          mu      = hydr%paramset(3) ;&
          psd_func = 6 * lwc*1.e-3 * (mu/Dmode)**(mu+4.) * D**mu * 1.e+12 &
                       * exp(-mu*D/Dmode) / (pi*density*gammaf(mu+4.))
!
! **  The Heymsfield and Platt (1984) distribution, given by the fit to 
! **  observed cloud ice distributions.
! **    Two options are available:
! **    nopt=0:   ice particle mass distributed in spherical
! **              volume with diameter equal to maximum crystal
! **              dimension (l).
! **    nopt=1:   ice particle described as pure ice sphere
! **              with same mass as elongated crystal.
!
          case ( 'HyPl84' ) ;&
          hydr%normaliz = .true. ;&
          hydr%indpndtT = .false. ;&
          hydr%twomomnt = .false. ;&
          nopt = nint(hydr%paramset(2)) ;&
          call heymplatt(nopt,T,lwc,D*0.5,density,num,1.) ;&
          psd_func = num ;&
          hydr%paramset(1) = density * 1.e+3 ! [g/cm3] -> [kg/m3]
!
! ** Exponential distribution for a 2-moment bulk microphysical scheme
! ** with specified particle density.
!
          case( 'exp2mo' ) ;&
          hydr%normaliz = .false. ;&
          hydr%indpndtT = .true. ;&
          hydr%twomomnt = .true. ;&
          density  = hydr%paramset(1) ;&
          lam      = (Nc*pi*density/(lwc*1.e-3))**0.333333 * 1.e-3 ;& ! [1/mm]
          n0       = Nc * (lam*1.e3) ;&
          psd_func = n0 * exp(-lam*D)
!
! ** Log-normal disbribution for a 2-moment bulk microphysical scheme
! ** with specified particle density and dispersion. 
!
          case( 'lgn2mo' ) ;&
          hydr%normaliz = .false. ;&
          hydr%indpndtT = .true. ;&
          hydr%twomomnt = .true. ;&
          density = hydr%paramset(1) ;& 
          sigma   = hydr%paramset(2) ;&
          Dmode   = ( 6*lwc*1.e-3 / (pi*density*Nc) )**0.333333 * &
                     exp( -1.5*sigma*sigma ) * 1.e+3 ;&   ! [mm]
          psd_func = normLN * lwc*1.e-3 / (density*sigma*(Dmode*1.e-3)**3) * &
                     exp( -4.5*sigma*sigma ) * &
                     exp( -0.5*((log(D)-log(Dmode))/sigma)**2 ) / &
                     (D*1.e-3)
!
! ** Error (undefined PSD)
!
          case default ;&
          print *, 'ERROR (MCRPparameters): '// &
                   'Function type undefined for '//hydr%name ;&
          stop
!
       End Select
!
     End function psd_func
!
 End Module MCRPparameters
