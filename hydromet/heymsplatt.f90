 Subroutine heymplatt(nopt,t,iwc,r,den,numden,norm)
   Implicit NONE
!
!  Return particle density and number density
!  of ice crystals, given temperature, ice water content of
!  ice crystal distribution, and maximum particle 
!  dimension.  Follows the empirical relation of 
!  Heymsfield and Platt (1984).

!  Input
!  Two options are available:
!  nopt=0:   ice particle mass distributed in spherical
! 	 	volume with diameter equal to maximum crystal
! 		dimension (l).
!  nopt=1:	ice particle described as pure ice sphere
! 		with same mass as elongated crystal.	  
!  t		temperature [K]
!  iwc	equivalent water content of cloud ice [g/m**3]
!  r		particle radius [mm]
!
!  Output
!  den	particle density [g/cm**3]
!  numden	particle number density [1/m**4]
!  norm	normalization factor required to contrain
! 		integrated distribution to equal the equivalent
! 		ice water content.  heymplatt must be called
! 		first with norm set equal to 1 to give
!            proper normalization
!
!  Bill Olson  April, 1998; modified by H. Masunaga, May, 2009.
!
   integer, intent(in) :: nopt
   real, intent(in)  :: t,iwc,r
   real, intent(out) :: den,numden
   real, intent(inout) :: norm
!
! Local variables
!
   integer :: n
   real, parameter :: lmin = 20.
   real, parameter :: lambda = 5.2e-2
   real :: l_lmin
   real, dimension(0:4) :: a = &
         (/-1.1430e+1,-7.3892e-1,-1.8647e-2,-1.4045e-4,0.0e+0/)
   real, dimension(0:4) :: b = &
         (/1.8940e+1,2.7658e0,1.2833e-1,2.7750e-3,2.2994e-5/)
   real, dimension(0:4) :: c = &
         (/-1.0159e+1,-1.4538e0,-1.3511e-2,1.1318e-3,2.2360e-5/)
   real, dimension(0:4) :: d = &
         (/1.6764e+1,-1.5072e-1,-1.9713e-2,-3.5051e-4,-1.6727e-6/)
   real, dimension(0:4) :: e = &
         (/1.5508e+2,1.8377e+1,8.5312e-1,1.6879e-2,1.1873e-4/)
   real, parameter :: pi = 3.14159265
   real :: suma,sumb,sumc,sumd,sume
   real :: tc,mass,l,l0,b1,b2,a1,a2
   real :: n100diwc,n1000diwc,lc,dc,dldr
   real, parameter :: densice = 0.917
   real :: dfac
!
   dfac = densice*3.*sqrt(3.)/2.
!
!  Temperature [C]
!
   tc=t-273.16
!
!  Note: empirical formulae only apply to range
!  -60 < tc < -20
!
   tc=max(min(-20.,tc),-60.)
!
!  Calculate maximum particle dimension [microns]
!  ( Rewritten by H. M. )
!
   mass=densice*4.*pi*((r*1.e-1)**3)/3.
! - for npot=1
   if(2.*r <= 0.3) then
      l = (mass/(dfac*((0.5/2.)**2)*1.e-3))**.33333
      dldr = (4./3.)  *(1.e-3)*densice*pi*r*r * l / mass
      l = l * 1.e+3
   else
      l = (mass/(dfac*((0.2/2.)**2)*1.e-3))**.54945
      dldr = (4./1.82)*(1.e-3)*densice*pi*r*r * l / mass
      l = l * 1.e+3
   end if
! - overridden when npot=0
   if ( nopt == 0 ) then
      l = 2.*r*1.e+3
      dldr = 2.0
   endif
!
! Extrapolation to small particles. See below. (H. M.)
!
   l_lmin = l - lmin
   l = max(l, lmin)
!
   suma=0. ; sumb=0. ; sumc=0. ; sumd=0. ; sume=0.
   do n=0,4
      suma=suma+a(n)*(tc**n)
      sumb=sumb+b(n)*(tc**n)
      sumc=sumc+c(n)*(tc**n)
      sumd=sumd+d(n)*(tc**n)
      sume=sume+e(n)*(tc**n)
   end do
   b1=suma ; b2=sumb
!
!   Liou's fit of b2 fails at low temperature; 
!   since Heymsfield and Platt (1984) data indicate
!   a nearly constant value of -4., we use it here.
!
   b2=-4.
!
   if(tc >= -37.5) then
      n100dIWC=exp(sumc)
   else
      n100dIWC=exp(sumd)
   end if
   n1000dIWC=sume
!
   a1=n100diwc/(100.**b1)
   a2=n1000diwc/(1000.**b2)
   l0=(a2/a1)**(1./(b1-b2))
!
! Rewritten by H. M. (dldr=2 for npot=0 as defined above)
!
   if(l <= l0) then
      numden=(1.e+6)*dldr*a1*(l**b1)*iwc/norm
   else
      numden=(1.e+6)*dldr*a2*(l**b2)*iwc/norm
   end if
!
   numden = numden/2 !! n(r)dr -> n(D)dD (added by H. M.)
!
! N is extrapolated by N0*exp(-lambda*D) to l<lmin=20 micron
! which is outside the range studied by HP84.
! lambda is here assumed to be 520 /cm = 5.2e-2 /micron.
! (Mitchell et al, JAS, 1996). Added by H. Masunaga.
!
   If ( l_lmin < 0. ) numden = numden * exp(-lambda*(l_lmin))
!
!  Compute particle density
!  Assume randomly oriented hexagonal column with
!  width/length relationship determined
!  from Heymsfield's empirical relation
!  ( rewritten by H. M. )
!
! for npot=1
   if ( nopt == 1 ) then
      den=0.917
      return
   end if
! for nopt=0
   if(l*1.e-3 <= 0.3) then
      dc=0.5*(l*1.e-3)*1.e-1
      lc=l*1.e-4
   else
      dc=0.2*((l*1.e-3)**0.41)*1.e-1
      lc=l*1.e-4
   end if
   mass=dfac*((dc/2.)**2.)*lc
   den=mass/(4.*pi*((r*1.e-1)**3.)/3.)
!
   return
 end Subroutine heymplatt

