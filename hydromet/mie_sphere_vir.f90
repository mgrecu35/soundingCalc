 Subroutine mie_sphere_vir (X, Mcmp, Qscat, Qexti, phsfn)
!
!    Mie Routine P. Bauer 
!    modified by H. Masunaga so the phase function is computed.
!
   Use VIRparameters, only: knang,nang,ang
   Implicit NONE
!
   integer, parameter ::  LIMITX = 1500
!
! I/O parameters
! 
   real(8), intent(in)  :: X
   real(8), intent(out) :: Qscat, Qexti
   real(8), dimension(knang), intent(out) :: phsfn
   complex(kind(0d0)), intent(in) :: Mcmp
!
! Local variables
!
   real(8) :: MR, MI, N1, N2
   real(8) :: Rfac1, Rfac2
   real(8) :: Rhelp1 (2), Rhelp2 (2)
   real(8), dimension(knang) :: amu, tau
   real(8), dimension(knang,0:LIMITX) :: PIn
!
   complex(kind(0d0)) :: M, MX
   complex(kind(0d0)) :: Chelp1, Chelp2, CFAC1, CFAC2
!
   complex(kind(0d0)), dimension(0:LIMITX)  :: Dn
   complex(kind(0d0)), dimension(-1:LIMITX) :: Wn
   complex(kind(0d0)), dimension(LIMITX)    :: An, Bn
   complex(kind(0d0)), dimension(knang)     :: S1, S2
!
   integer :: NEnd
   integer :: I100, i
   real(8), parameter :: pid = 3.141592653589793d0
!
   equivalence (Chelp1, Rhelp1 (1))
   equivalence (Chelp2, Rhelp2 (1))
!
   M      = conjg (Mcmp)
   Chelp1 = M
   MR     =          Rhelp1 (1)
   MI     = -1.0d0 * Rhelp1 (2)
!      
   MX   = M  * X
   N1   = MR * X
   N2   = MI * X
!
   If (X <= 20000.0)  NEnd = X + 4.00 * X ** (1.0 / 3.0) + 2.0
   If (X <=  4200.0)  NEnd = X + 4.05 * X ** (1.0 / 3.0) + 2.0
   If (X <=     8.0)  NEnd = X + 4.00 * X ** (1.0 / 3.0) + 1.0
   If (NEnd <=    5)  NEnd = 5
   If (NEnd > LIMITX) NEnd = LIMITX
!
! Modified by H. M. 
!
   If ( N2 < 30.d0 ) then
!
      Rfac1      = sin  (N1) * sin  (N1) + sinh (N2) * sinh (N2)
      Rhelp1 (1) = sin  (N1) * cos  (N1) / Rfac1
      Rhelp1 (2) = sinh (N2) * cosh (N2) / Rfac1
!
   Else
!
      Rhelp1 (1) = 0.d0
      Rhelp1 (2) = 1.d0
!
   Endif
!
   Dn (0) = Chelp1
!
   Rhelp1 (1) =            cos (X)
   Rhelp1 (2) = -1.0d+00 * sin (X)
   Rhelp2 (1) =            sin (X)
   Rhelp2 (2) =            cos (X)
!
   Wn (-1) = Chelp1
   Wn ( 0) = Chelp2
!
   Qexti  = 0.d0 ; Qscat  = 0.d0
   amu(1:nang) = cos(ang(1:nang)*pid/180.d0)
   S1 = dcmplx (0.d0,0.d0) ; S2 = dcmplx (0.d0,0.d0) ; PIn = 0.d0
!
   do I100 = 1, NEnd
      Dn (I100) = -1.0d0 * I100 / MX &
                +  1.0d0 / (I100 / MX - Dn (I100 - 1))
      Wn (I100) = Wn (I100 - 1) * (2.0d0 * I100 - 1.0d0) / X &
                - Wn (I100 - 2)
!
      CFAC1 = Dn (I100) / M + I100 / X
      CFAC2 = M * Dn (I100) + I100 / X
!
      Chelp1 = Wn (I100)
      Chelp2 = Wn (I100 - 1)
!
      An (I100) = (CFAC1 * Rhelp1 (1) - Rhelp2 (1)) &
                / (CFAC1 * Chelp1     - Chelp2    )
      Bn (I100) = (CFAC2 * Rhelp1 (1) - Rhelp2 (1)) &
                / (CFAC2 * Chelp1     - Chelp2    )
!
      Chelp1 = An (I100)
      Chelp2 = Bn (I100)
!
      Rfac1 = Rhelp1 (1) + Rhelp2 (1)
      Rfac2 = cdabs (An (I100)) * cdabs (An (I100)) &
            + cdabs (Bn (I100)) * cdabs (Bn (I100))
!
      Qexti  = Qexti  + (2.0d0 * I100 + 1.0d0) * Rfac1
      Qscat  = Qscat  + (2.0d0 * I100 + 1.0d0) * Rfac2
!
! / Phase function (added by H. M.) /
! based on Eqs. (4.47), (4.74), and (4.77) from Bohren and Huffman (1983)
! phsf = S_11 divided by normalization.
!
      Rfac1 = (2.d0 * I100 + 1.d0) / (I100 * (1.d0 * I100 + 1.d0))
!
      select case ( I100 )
         case (1)  ; PIn(1:nang,I100) = 1.d0
         case (2:) ; PIn(1:nang,I100) = &
                          (2.d0 * I100 - 1.d0) / (1.d0 * I100 - 1.d0) * &
                            PIn(1:nang,I100-1) * amu(1:nang) - I100 * &
                            PIn(1:nang,I100-2) / (1.d0 * I100 - 1.d0)
      end select
      tau(1:nang) = I100 * amu(1:nang) * PIn(1:nang,I100) - & 
                            (I100 + 1) * PIn(1:nang,I100-1)
      S1 (1:nang) = S1(1:nang) + Rfac1 * (An(I100) * PIn(1:nang,I100) + &
                                          Bn(I100) * tau(1:nang))
      S2 (1:nang) = S2(1:nang) + Rfac1 * (An(I100) * tau(1:nang) + &
                                          Bn(I100) * PIn(1:nang,I100))
!
   end do
!
   phsfn(1:nang) = 0.5d0 * (cdabs(S1(1:nang))**2 + cdabs(S2(1:nang))**2)
!
   Qscat = min(Qscat,Qexti)
   phsfn = phsfn / (2 * pid * Qscat)
!
   Qexti = Qexti * 2.0d0 / (X * X)
   Qscat = Qscat * 2.0d0 / (X * X)
!
    return
 end Subroutine mie_sphere_vir
