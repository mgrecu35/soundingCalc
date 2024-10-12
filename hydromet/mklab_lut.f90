! ############################################################
!            Satellite Data Simulation Unit v2
!                    - mklab_lut -
!
! This function creates a label attached to the LUT file name
! based on the channel frequencies for use by the microwave 
! radiometer and radar simulators. 
!
!                                 created by H. M. (May, 2009)
!
 Character*4 function mklab_lut ( freq )
   Implicit NONE
!
! freq: channel frequency in [GHz]
!
   real, intent(in) :: freq
   character*3 :: i2a
!
   mklab_lut = i2a( nint(freq),3 )//'G'
!
 End function mklab_lut
!
! ############################################################
!            Satellite Data Simulation Unit v2
!                   - mklab_lut_vir -
!
! This function creates a label attached to the LUT file name
! based on the channel wavelenghts for use by the visible/IR
! simulator. 
!
!                                 created by H. M. (May, 2009)
!
 Character*4 function mklab_lut_vir ( wvln )
   Implicit NONE
!
! freq: channel frequency in [GHz]
!
   real, intent(in) :: wvln
   character*3 :: i2a
!
   If ( wvln < 1.0 ) then
      mklab_lut_vir = i2a( nint(wvln*1.e3),3 )//'N'
   Else
      mklab_lut_vir = i2a( nint(wvln*1.e1),3 )//'M'
   End If
!
 End function mklab_lut_vir
!
! #######################################################
!               !!!  Character Function i2a  !!!
! Convert an integer variable "number" to the equivalent
! character variable. If the total digit is specified by
! "idigit", excessive digits will be filled with 0s
! (e.g., 3 will be 003 for idigit=3). Substitute zero
! for "idigit" to let the routine find the total digit.
!                                  H. Masunaga (02/03/05)
!       modified to f90 fashion by H. Masunaga (06/27/07)
!
 Function i2a ( number, idigit )
   Implicit NONE
!
! * I/O variables
!
   Integer     :: number
   Character*3 :: i2a
!
! * Internal variables
!
   Integer      :: i,ibuf,idigit,totdig
   Character*11 :: cbuf
!
! Find the total digit of a given integer
!
   Do totdig = 1, 11
      If ( number/10**totdig == 0 ) exit
   End do
   If ( idigit == 0 ) then
      idigit = totdig
   Else if ( idigit < totdig ) then
      print *, 'Warning (i2a): idigit too small'
      idigit = totdig
   Endif
!
   cbuf = char(mod(number/10**(idigit-1),10) + ichar('0'))
   Do i = 2, idigit
      ibuf = mod(number/10**(idigit-i),10) + ichar('0')
      cbuf = cbuf(1:i-1)//char(ibuf)
   End do
   i2a = cbuf(1:idigit)
!
   Return
 End Function i2a
