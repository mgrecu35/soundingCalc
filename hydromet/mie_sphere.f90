 Subroutine mie_sphere (X, Mcmp, Qscat, Qexti, asym, Qbscat)
!
!    Mie Routine P. Bauer 
!
   Implicit NONE
!
   integer, parameter ::  LIMITX = 1500
!
! I/O parameters
! 
   real(8), intent(in)  :: X
   real(8), intent(out) :: Qscat, Qexti, asym, Qbscat
   complex(kind(0d0)), intent(in) :: Mcmp
!
! Local variables
!
   real(8) :: MR, MI, N1, N2
   real(8) :: Rfac1, Rfac2
   real(8) :: Rhelp1 (2), Rhelp2 (2)
!
   complex(kind(0d0)) :: M, MX
   complex(kind(0d0)) :: Chelp1, Chelp2, CFAC1, CFAC2, Cbscat
!
   complex(kind(0d0)), dimension(0:LIMITX)  :: Dn
   complex(kind(0d0)), dimension(-1:LIMITX) :: Wn
   complex(kind(0d0)), dimension(LIMITX)    :: An, Bn
!
   integer :: NEnd
   integer :: I100, I101
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
!   Rfac1      = sin  (N1) * sin  (N1) + sinh (N2) * sinh (N2)
!   Rhelp1 (1) = sin  (N1) * cos  (N1) / Rfac1
!   Rhelp1 (2) = sinh (N2) * cosh (N2) / Rfac1
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
   Qexti  = 0.d0 ; Qscat  = 0.d0 ; Qbscat = 0.d0 ; asym   = 0.d0 
   Cbscat = dcmplx (0.0,0.0)
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
      Cbscat = Cbscat + (2.0d0 * I100 + 1.0d0) * (-1.0d0) ** I100 &
             * (An (I100) - Bn (I100))
!
      If (I100 == 1) cycle
!
      Chelp1 = An (I100 - 1) * conjg (An (I100)) &
             + Bn (I100 - 1) * conjg (Bn (I100))
      Chelp2 = An (I100 - 1) * conjg (Bn (I100 - 1))
!
      I101 = I100 - 1
      Rfac1  = I101 * (I101 + 2) / (I101 + 1.0d0)
      Rfac2  = (2.0d0 * I101 + 1.0d0) / (I101 * (I101 + 1.0d0))
!
      asym = asym + Rfac1 * Rhelp1 (1) + Rfac2 * Rhelp2 (1)
   end do
!
   Qexti  = Qexti * 2.0d0 / (X * X)
   Qscat  = Qscat * 2.0d0 / (X * X)
   asym   = asym  * 4.0d0 / (X * X * Qscat)
   Qbscat = cdabs (Cbscat) * cdabs (Cbscat) / (X * X)
   If (Qscat > Qexti) Qscat = Qexti
!
   return
 end Subroutine mie_sphere
