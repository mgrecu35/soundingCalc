 SUBROUTINE ICEOPTIC (FGHZ, T, EPSRI, EPSII)
!
!  Hufford (1991), see Brussard and Watson (1995), p.297
!
   IMPLICIT NONE
!
!  SAVE
   DOUBLE PRECISION :: FGHZ, EPSRI, EPSII, A, B, THETA, T, TK
!
   EPSRI = 3.15D+00
!
   TK = T
   IF (TK > 273.16) TK = 273.16
!
   THETA = 300.0 / TK
   A     = 1.0D-04 * (50.4 + 62.0 * (THETA - 1.0)) &
                   * EXP (-22.1 * (THETA - 1.0))
   B     = 1.0D-04 * (0.633 / THETA - 0.131) &
       + (7.36D-04 * THETA / (THETA - 0.9927)) &
       * (7.36D-04 * THETA / (THETA - 0.9927))
   EPSII = A / FGHZ + B * FGHZ
!
   RETURN
 END SUBROUTINE ICEOPTIC
