!
! / Gamma function /
!
 Real function gammaf(x)
   real, intent(in) :: x
   real(8) :: gammln,xx
!
   xx = dble(x)
   gammaf = exp(gammln(xx))
!
 End function gammaf
!
 Function gammln(xx)
   real(8) ::  gammln,xx
   integer :: j
   real(8) :: ser,tmp,x,y
   real(8), dimension(6), save :: cof = (/ &
          76.18009172947146d0,-86.50532032941677d0, &
          24.01409824083091d0,-1.231739572450155d0,&
         .1208650973866179d-2,-.5395239384953d-5 /)
   real(8), save :: stp = 2.5066282746310005d0
!
   x=xx
   y=x
   tmp=x+5.5d0
   tmp=(x+0.5d0)*log(tmp)-tmp
   ser=1.000000000190015d0
   do j=1,6
      y=y+1.d0
      ser=ser+cof(j)/y
   end do
   gammln=tmp+Dlog(stp*ser/x)
!
   return
 End Function gammln
