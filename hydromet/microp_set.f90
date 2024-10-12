! ############################################################
!            Satellite Data Simulation Unit v2
!                    -  microp_set -
!
! This subroutine defines microphysical setups.
!
!                                            H. M. (May, 2009)
!
 Subroutine microp_set
   Use MCRPparameters
   Implicit NONE
! 
   real :: init
!
! ** Exponential distribution with specified particle density and
! ** intercept parameter.
! **               N(D)dD = N0 * exp(-lambda*D)dD,
! ** where lambda is determined from the given water content.
!
!   functype = 'expntl', nparams = 2,
!   paramset(1) = particle density in [kg/m^3]
!   paramset(2) = Intercept parameter (N0) in [1/m^4]
!
! ** Gamma distribution with specified particle density, mu, and
! ** normalized intercept parameter, Nw.
! **    N(D)dD = Nw * f(mu) * (D/D0)^mu * exp(-(3.67+mu)*D/D0)dD
! ** where D0 is determined from the given water content.
! ** The "normalized" intercept parameter, Nw, equals N0 of the
! ** exponential DSD that has the same water content and D0. 
! ** See Eq. (7.62) of Bringi and Chandrasekar (2001).
!
!   functype = 'gammwN', nparams = 3,
!   paramset(1) = particle density in [kg/m^3]
!   paramset(2) = normalized intercept parameter (Nw) in [1/m^4]
!   paramset(3) = mu
!
! ** Gamma distribution with specified particle density, mu, and
! ** median volume diameter, D0.
! **    N(D)dD = Nw * f(mu) * (D/D0)^mu * exp(-(3.67+mu)*D/D0)dD
! ** where Nw is determined from the given water content.
!
!   functype = 'gammwD', nparams = 3,
!   paramset(1) = particle density in [kg/m^3]
!   paramset(2) = median volume diameter in [mm]
!   paramset(3) = mu
!
! ** Log-normal distribution with specified particle density,
! ** mode diameter(Dm), and dispersion (sigma). 
! **  N(D)dD = Nc/sqrt(2*pi)/sigma * exp[-((ln(D)-ln(Dm))/sigma)**2/2]dD/D
! ** where Nc is determined numerically from the given water content.
!
!   functype = 'lgnorm', nparams = 3,
!   paramset(1) = particle density in [kg/m^3]
!   paramset(2) = mode diameter in [mm]
!   paramset(2) = dispersion in terms of log(D)
!
! ** Exponential PSD with a power-law mass spectrum m(D) with specified
! ** intercept parameter, 'a', and 'b', where m(D) = a*D^b
! **            N(D)dD = N0 * exp(-lambda*D)dD
! ** Particle density is 6/pi * a * D^(b-3) and lambda is determined
! ** from the given water content
!
!   functype = 'expaDb', nparams = 4,
!   paramset(1) = not used (overwritten by calculated particle density)
!   paramset(2) = Intercept parameter (N0) in [1/m^4]
!   paramset(3) = 'a' in [kg/m^b]
!   paramset(4) = 'b'
!
! ** Modified gamma distribution, suitable for a liquid cloud DSD,
! ** with specified particle density and mode radius (rm).
! **        N(D)dD = N*6^6/(5!*Dm)*(D/Dm)^6*exp(-6*D/Dm)dD
! ** where N is determined from the given water content.
! ** See Eq. (4.3.3) of Liou (1992).
!
!   functype = 'modGam', nparams = 2,
!   paramset(1) = particle density in [kg/m^3]
!   paramset(2) = mode diameter in [mm]
!
! **  The Heymsfield and Platt (1984) distribution, given by the fit to 
! **  observed cloud ice distributions. The scaling factor is
! **  determined numerically from the given water content.
! **    Two options are available:
! **    nopt=0:   ice particle mass distributed in spherical
! **              volume with diameter equal to maximum crystal
! **              dimension (l).
! **    nopt=1:   ice particle described as pure ice sphere
! **              with same mass as elongated crystal.
!
!   functype = 'HyPl84', nparams = 2,
!   paramset(1) = not used (overwritten by calculated particle density)
!   paramset(2) = nopt (0. or 1.)
!
! ** 2-moment exponential distribution with a specified particle density
! **               N(D)dD = N0 * exp(-lambda*D)dD,
! **               N0 = Nc / lambda
! ** where lambda and N0 are determined from the given water content
! ** and hydrometeor number concentration, Nc.
!
!   functype = 'exp2mo', nparams = 1,
!   paramset(1) = particle density in [kg/m^3]
!
!        //////////    Default Setup    ///////////
!
! - rain DSD
!
   rain%name = 'rain'
   rain%functype = 'expntl' ; rain%nparams = 2
   rain%paramset(1:rain%nparams) = (/ densliq, 2.2e+7 /)
!
! - snow PSD
!
   snow%name = 'snow'
   snow%functype = 'expntl' ; snow%nparams = 2
   snow%paramset(1:snow%nparams) = (/  0.1e+3, 1.0e+8 /)
!
! - graupel PSD
!
   grpl%name = 'grpl'
   grpl%functype = 'expntl' ; grpl%nparams = 2
   grpl%paramset(1:grpl%nparams) = (/  0.4e+3, 4.0e+6 /)
!
! - hail PSD
!
   hail%name = 'hail'
   hail%functype = 'expntl' ; hail%nparams = 2
   hail%paramset(1:hail%nparams) = (/ densice, 4.0e+6 /)
!
! - Cloud water DSD
!
   cldw%name = 'cldw'
   cldw%functype = 'lgnorm' ; cldw%nparams = 3
   cldw%paramset(1:cldw%nparams) = (/ densliq, 2.e-2, 0.35 /)
!
! - Cloud ice PSD
!
   cldi%name = 'cldi'
   cldi%functype = 'HyPl84' ; cldi%nparams = 2
   cldi%paramset(1:cldi%nparams) = (/ densice, 1. /)
!
!       //////////    CReSS microphysics    ///////////
!
! - rain DSD
!
!   rain%name = 'rainCReS'
!   rain%functype = 'expntl' ; rain%nparams = 2
!   rain%paramset(1:rain%nparams) = (/ densliq, 8.0e+6 /)
!
! - snow PSD
!
!   snow%name = 'snowCReS'
!   snow%functype = 'exp2mo' ; snow%nparams = 1
!   snow%paramset(1:snow%nparams) = (/  8.4e+1 /)
!
! - graupel PSD
!
!   grpl%name = 'grplCReS'
!   grpl%functype = 'exp2mo' ; grpl%nparams = 1
!   grpl%paramset(1:grpl%nparams) = (/  3.0e+2 /)
!
! - hail PSD (not used for CReSS)
!
!   hail%name = 'hail'
!   hail%functype = 'expntl' ; hail%nparams = 2
!   hail%paramset(1:hail%nparams) = (/ densice, 4.0e+6 /)
!
! - Cloud water DSD
!
!   cldw%name = 'cldwCReS'
!   cldw%functype = 'lgn2mo' ; cldw%nparams = 2
!   cldw%paramset(1:cldw%nparams) = (/ densliq, 0.35 /)
!
! - Cloud ice PSD
!
!   cldi%name = 'cldiCReS'
!   cldi%functype = 'lgn2mo' ; cldi%nparams = 2
!   cldi%paramset(1:cldi%nparams) = (/ densice, 0.35 /)
!
!
! - /////////  Backup /////////
!
!  / alternative rain DSD with a gamma distribution /
!   rain%name = 'rainBC01'
!   rain%functype = 'gammwN' ; rain%nparams = 3
!   rain%paramset(1:rain%nparams) = (/ densliq, 1.5e+7, 3. /)
!
!  / alternative snow PSD with the Grabowski (1998) scheme /
!   snow%name = 'snowG98'
!   snow%functype = 'expaDb' ; snow%nparams = 4
!   snow%paramset(1:snow%nparams) = (/ densice, 1.e+7, 2.5e-2, 2. /)
!
!  / alternative cloud water DSD with a modified gamma distribution /
!   cldw%name = 'cldwL92'
!   cldw%functype = 'modGam' ; cldw%nparams = 3
!   cldw%paramset(1:cldw%nparams) = (/ densliq, 8.e-3, 6. /)
!
!  / alternative cloud ice PSD with a gamma distribution /
!   cldi%name = 'cldiES95'
!   cldi%functype = 'gammwD' ; cldi%nparams = 3
!   cldi%paramset(1:cldi%nparams) = (/ densice, 1.e-1, 1. /)
!
!
! ####### initialization #######
!   ///// constants /////
!
   eair = cmplx(1.0006,0.0)     ! air refractive index
   normLN = 6./(pi*sqrt(2*pi))  ! a PSD normalization factor
!
!   ///// activate hydr%normaliz /////
!
   init = psd_func ( cldw,1.,300.,0.,1. )
   init = psd_func ( cldi,1.,300.,0.,1. )
   init = psd_func ( rain,1.,300.,0.,1. )
   init = psd_func ( snow,1.,300.,0.,1. )
   init = psd_func ( grpl,1.,300.,0.,1. )
   init = psd_func ( hail,1.,300.,0.,1. )
!
 End Subroutine microp_set
