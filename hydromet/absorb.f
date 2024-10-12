      subroutine absorb(freq,temp,pres,relhum,water,kabs)
 
      real  kabs, kabs_H2O, kabs_O2, kabs_N2, kabs_clw
      real  pwvsat
      
      ! freq   [GHz]
      ! pres   [mb]
      ! RH     [%]
      ! temp   [K]
      ! water  [g/m^3]
!
! / update ES for SDSU v2.1.0 /
!
      ES = pwvsat ( temp )
      vapor_pressure =  relhum*ES/100.      
      rho = vapor_pressure*100*18/(8.314*temp)

      call abs_H2O(temp, pres, rho, freq, kabs_H2O)      
 
      call abs_O2(temp, pres, rho, freq, kabs_O2)
 
      call abs_N2(temp, pres, freq, kabs_N2)

      call absorb_clw(freq,temp,water,kabs_clw) 

      kabs = kabs_H2O + kabs_O2 + kabs_N2 + kabs_clw

      return
      end
      
      Subroutine abs_H2O(T,P,RHO,F,ABH2O)
C
C     C. Kummerow, 8/2003.  Changed function to subroutine     
C     Copyright (c) 2002 Massachusetts Institute of Technology
C
C  NAME- ABH2O    LANGUAGE- FORTRAN 77
C
C PURPOSE- COMPUTE ABSORPTION COEF IN ATMOSPHERE DUE TO WATER VAPOR
C 
      IMPLICIT NONE
C  CALLING SEQUENCE PARAMETERS-
C    SPECIFICATIONS
      REAL T,P,RHO,F,ABH2O
C      NAME    UNITS    I/O  DESCRIPTON            VALID RANGE
C      T       KELVIN    I   TEMPERATURE
C      P       MILLIBAR  I   PRESSURE              .1 TO 1000
C      RHO     G/M**3    I   WATER VAPOR DENSITY
C      F       GHZ       I   FREQUENCY             0 TO 800
C      ABH2O   NEPERS/KM O   ABSORPTION COEFFICIENT
C
C   REFERENCES-
C   P.W. ROSENKRANZ, RADIO SCIENCE V.33, PP.919-928 (1998); V.34, P.1025 (1999).
C
C   LINE INTENSITIES SELECTION THRESHOLD=
C     HALF OF CONTINUUM ABSORPTION AT 1000 MB.
C   WIDTHS MEASURED AT 22,183,380 GHZ, OTHERS CALCULATED.
C     A.BAUER ET AL.ASA WORKSHOP (SEPT. 1989) (380GHz).
c     M. TRETYAKOV et al., J. MOLEC. SPEC. (2003)
C
C   REVISION HISTORY-
C    DATE- OCT.6, 1988  P.W.ROSENKRANZ - EQS AS PUBL. IN 1993.
C          OCT.4, 1995  PWR- USE CLOUGH'S DEFINITION OF LOCAL LINE
C                   CONTRIBUTION,  HITRAN INTENSITIES, ADD 7 LINES.
C          OCT. 24, 95  PWR -ADD 1 LINE.
C          JULY 7, 97   PWR -SEPARATE COEFF. FOR SELF-BROADENING, 
C                       REVISED CONTINUUM.
C        Aug. 28, 2002  PWR - CORRECTED LINE INTENSITIES
C        Mar. 2, 2003   PWR - LINE SHIFT
C
C   LOCAL VARIABLES:
      INTEGER NLINES,I,J
      PARAMETER (NLINES=15)
      REAL DF(2),S1(NLINES),B2(NLINES),W3(NLINES),FL(NLINES),X(NLINES),
     & WS(NLINES),XS(NLINES),SR(NLINES)
      REAL PVAP,PDA,DEN,TI,TI2,SUM,WIDTH,WSQ,S,BASE,RES,CON,SHIFT
C     LINE FREQUENCIES:
      DATA FL/22.2351, 183.3101, 321.2256, 325.1529, 380.1974, 439.1508,
     & 443.0183, 448.0011, 470.8890, 474.6891, 488.4911, 556.9360,
     & 620.7008, 752.0332, 916.1712/
C     LINE INTENSITIES AT 300K:
      DATA S1/ .1314E-13, .2279E-11, .8058E-13, .2701E-11, .2444E-10,
     & .2185E-11, .4637E-12, .2568E-10, .8392E-12, .3272E-11, .6676E-12,
     & .1535E-08, .1711E-10, .1014E-08, .4238E-10/
C     T COEFF. OF INTENSITIES:
      DATA B2/ 2.144, .668, 6.179, 1.541, 1.048, 3.595, 5.048, 1.405,
     & 3.597, 2.379, 2.852, .159, 2.391, .396, 1.441/
C     AIR-BROADENED WIDTH PARAMETERS AT 300K:
      DATA W3/.00281, .00287, .0023, .00278, .00287, .0021, .00186,
     & .00263, .00215, .00236, .0026, .00321, .00244, .00306, .00267/
C     T-EXPONENT OF AIR-BROADENING:
      DATA X/.69, .64, .67, .68, .54, .63, .60, .66, .66, .65, .69, .69,
     & .71, .68, .70/
C     SELF-BROADENED WIDTH PARAMETERS AT 300K:
      DATA WS/.01349, .01491, .0108, .0135, .01541, .0090, .00788,
     & .01275, .00983, .01095, .01313, .01320, .01140, .01253, .01275/
C     T-EXPONENT OF SELF-BROADENING:
      DATA XS/ .61, .85, .54, .74, .89, .52, .50, .67, .65, .64, .72,
     & 1.0, .68, .84, .78/
C     RATIO OF SHIFT TO WIDTH
      DATA SR/ 0., -.017, 13*0./
C
      IF(RHO.LE.0.) THEN
        ABH2O = 0.
        RETURN
      ENDIF
      PVAP = RHO*T/217.
      PDA = P -PVAP
      DEN = 3.335E16*RHO ! const includes isotopic abundance
      TI = 300./T
      TI2 = TI**2.5
C
C      CONTINUUM TERMS
      CON = (5.43E-10*PDA*TI**3 + 1.8E-8*PVAP*TI**7.5)*PVAP*F*F 
C
C      ADD RESONANCES
      SUM = 0.
      DO 30 I=1,NLINES
      WIDTH = W3(I)*PDA*TI**X(I) + WS(I)*PVAP*TI**XS(I)
      SHIFT = SR(I)*WIDTH  ! unknown temperature dependence
      WSQ = WIDTH*WIDTH
      S = S1(I)*TI2*EXP(B2(I)*(1.-TI))
      DF(1) = F - FL(I) - SHIFT
      DF(2) = F + FL(I) + SHIFT
C  USE CLOUGH'S DEFINITION OF LOCAL LINE CONTRIBUTION
      BASE = WIDTH/(562500. + WSQ)
C  DO FOR POSITIVE AND NEGATIVE RESONANCES
      RES = 0.
      DO 20 J=1,2
      IF(ABS(DF(J)).LT.750.) RES = RES + WIDTH/(DF(J)**2+WSQ) - BASE
20    CONTINUE
      SUM = SUM + S*RES*(F/FL(I))**2
30    CONTINUE
      ABH2O = .3183E-4*DEN*SUM + CON
      RETURN
      END
      
      subroutine abs_O2(TEMP,PRES,VAPDEN,FREQ,O2ABS)
C
C     C. Kummerow, 8/2003.  Changed function to subroutine      
C  Copyright (c) 2003 Massachusetts Institute of Technology
C
C     PURPOSE: RETURNS ABSORPTION COEFFICIENT DUE TO OXYGEN IN AIR,
C              IN NEPERS/KM
C
C      5/1/95  P. Rosenkranz 
C      11/5/97  P. Rosenkranz - 1- line modification.
c      12/16/98 pwr - updated submm freq's and intensities from HITRAN96
c      8/21/02  pwr - revised width at 425
c      3/20/03  pwr - 1- line mixing and width revised
C
c     IMPLICIT NONE
C
C     ARGUMENTS:
      REAL TEMP,PRES,VAPDEN,FREQ
C
C     NAME    UNITS    DESCRIPTION        VALID RANGE
C
C     TEMP    KELVIN   TEMPERATURE        UNCERTAIN, but believed to be
c                                          valid for atmosphere
C     PRES   MILLIBARS PRESSURE           3 TO 1000
C     VAPDEN  G/M**3   WATER VAPOR DENSITY  (ENTERS LINEWIDTH CALCULATION
C                      DUE TO GREATER BROADENING EFFICIENCY OF H2O)
C     FREQ    GHZ      FREQUENCY          0 TO 900
C
C     REFERENCES FOR EQUATIONS AND COEFFICIENTS:
C     P.W. Rosenkranz, CHAP. 2 and appendix, in ATMOSPHERIC REMOTE SENSING
C      BY MICROWAVE RADIOMETRY (M.A. Janssen, ed., 1993).
C     H.J. Liebe et al, JQSRT V.48, pp.629-643 (1992).
c     M.J. Schwartz, Ph.D. thesis, M.I.T. (1998).
c     A.F. Krupnov et al, J. Mol. Spect. v.215, pp.309-311 (2002).
C     M.Yu. Tretyakov et al, J. Mol. Spect. (2003 preprint).
C     SUBMILLIMETER LINE INTENSITIES FROM HITRAN96.
c
c     This version differs from Liebe's MPM92 in these significant respects:
c     1. The 1- line has the width and mixing coefficient measured by 
c      Tretyakov et al. 
c     2. It modifies the 1- line width temperature dependence to (1/T)**0.9
c     3. It uses the same temperature dependence (X) for submillimeter 
c      line widths as in the 60 GHz band: (1/T)**0.8 
c     4. The 425 GHz line width is from Krupnov et al.
C
c     Local variables:
      REAL TH,TH1,B,PRESWV,PRESDA,DEN,DENS,DFNR,SUM,STR,Y,SF1,SF2,FCEN
      INTEGER K
      REAL X,WB300,W300(40),F(40),Y300(40),S300(40),V(40),BE(40)
      COMMON /O2COM/ X,WB300,W300,F,Y300,S300,V,BE
C      LINES ARE ARRANGED 1-,1+,3-,3+,ETC. IN SPIN-ROTATION SPECTRUM
      DATA F/118.7503, 56.2648, 62.4863, 58.4466, 60.3061, 59.5910,
     2  59.1642, 60.4348, 58.3239, 61.1506, 57.6125, 61.8002,
     3  56.9682, 62.4112, 56.3634, 62.9980, 55.7838, 63.5685,
     4  55.2214, 64.1278, 54.6712, 64.6789, 54.1300, 65.2241,
     5  53.5957, 65.7648, 53.0669, 66.3021, 52.5424, 66.8368,
     6  52.0214, 67.3696, 51.5034, 67.9009, 368.4984, 424.7632,
     7  487.2494, 715.3931, 773.8397, 834.1458/
        DATA S300/.2936E-14,.8079E-15, .2480E-14,.2228E-14,
     &  .3351E-14,.3292E-14, .3721E-14,.3891E-14,
     &  .3640E-14,.4005E-14, .3227E-14,.3715E-14,
     &  .2627E-14,.3156E-14, .1982E-14,.2477E-14,
     &  .1391E-14,.1808E-14, .9124E-15,.1230E-14,
     &  .5603E-15,.7842E-15, .3228E-15,.4689E-15,
     &  .1748E-15,.2632E-15, .8898E-16,.1389E-15,
     &  .4264E-16,.6899E-16, .1924E-16,.3229E-16,
     &  .8191E-17,.1423E-16, .6494E-15, .7083E-14, .3025E-14,
     &  .1835E-14, .1158E-13, .3993E-14/
      DATA BE/.009,.015, .083,.084, 2*.212, 2*.391, 2*.626,
     & 2*.915, 2*1.260, 1.660,1.665, 2.119,2.115, 2.624,2.625,
     & 2*3.194, 2*3.814, 2*4.484, 2*5.224, 2*6.004, 2*6.844,
     & 2*7.744, .048, .044, .049, .145, .141, .145/
C      WIDTHS IN MHZ/MB
      DATA WB300/.56/, X/.8/
      DATA W300/1.67, 1.646, 1.468, 1.449, 1.382, 1.360,
     & 1.319, 1.297, 1.266, 1.248, 1.221, 1.207, 1.181, 1.171,
     & 1.144, 1.139, 1.110, 1.108, 1.079, 1.078, 2*1.05,
     & 2*1.02,2*1.00,2*.97,2*.94,2*.92,2*.89, 3*1.64, 3*1.81/
      DATA Y300/  -0.036,  0.2408, -0.3486,  0.5227,
     & -0.5430,  0.5877, -0.3970,  0.3237, -0.1348,  0.0311,
     &  0.0725, -0.1663,  0.2832, -0.3629,  0.3970, -0.4599,
     &  0.4695, -0.5199,  0.5187, -0.5597,  0.5903, -0.6246,
     &  0.6656, -0.6942,  0.7086, -0.7325,  0.7348, -0.7546,
     &  0.7702, -0.7864,  0.8083, -0.8210,  0.8439, -0.8529, 6*0./
      DATA V/  0.0079, -0.0978,  0.0844, -0.1273,
     &  0.0699, -0.0776,  0.2309, -0.2825,  0.0436, -0.0584,
     &  0.6056, -0.6619,  0.6451, -0.6759,  0.6547, -0.6675,
     &  0.6135, -0.6139,  0.2952, -0.2895,  0.2654, -0.2590,
     &  0.3750, -0.3680,  0.5085, -0.5002,  0.6206, -0.6091,
     &  0.6526, -0.6393,  0.6640, -0.6475,  0.6729, -0.6545, 6*0./
C
      TH = 300./TEMP
      TH1 = TH-1.
      B = TH**X
      PRESWV = VAPDEN*TEMP/217.
      PRESDA = PRES -PRESWV
      DEN = .001*(PRESDA*B + 1.1*PRESWV*TH)
      DENS = .001*(PRESDA*TH**.9 + 1.1*PRESWV*TH)
      DFNR = WB300*DEN
      SUM = 1.6E-17*FREQ*FREQ*DFNR/(TH*(FREQ*FREQ + DFNR*DFNR))
      DO K=1,40
      IF(K.EQ.1) THEN !exception for 1- line
        DF = W300(1)*DENS
      ELSE
        DF = W300(K)*DEN
      ENDIF
      FCEN = F(K)
      Y = .001*PRES*B*(Y300(K)+V(K)*TH1)
      STR = S300(K)*EXP(-BE(K)*TH1)
      SF1 = (DF + (FREQ-FCEN)*Y)/((FREQ-FCEN)**2 + DF*DF)
      SF2 = (DF - (FREQ+FCEN)*Y)/((FREQ+FCEN)**2 + DF*DF)
      SUM = SUM + STR*(SF1+SF2)*(FREQ/F(K))**2
      END DO
      O2ABS = .5034E12*SUM*PRESDA*TH**3/3.14159
      O2ABS = AMAX1(O2ABS,0.)
      RETURN
      END

      Subroutine ABS_N2(T,P,F,ABSN2)

C
C     C. Kummerow, 8/2003.  Changed function to subroutine      
C  Copyright (c) 2002 Massachusetts Institute of Technology
C     ABSN2 = COLLISION-INDUCED ABSORPTION COEFFICIENT (NEPER/KM)

C     IN AIR

C     T = TEMPERATURE (K)

C     P = PRESSURE (MB)

C     F = FREQUENCY (GHZ)(valid 0-1000 GHz)

C

c     5/22/02 P.Rosenkranz

c

C     Equations based on:

C      Borysow, A, and L. Frommhold, 

C      Astrophysical Journal, v.311, pp.1043-1057 (1986)

C     with modification of 1.29 to account for O2-O2 and O2-N2

c     collisions, as suggested by

C      J.R. Pardo, E.Serabyn, J.Cernicharo, J. Quant. Spectros.

c      Radiat. Trans. v.68, pp.419-433 (2001).

c

      TH = 300./T

      FDEPEN = .5 + .5/(1.+(F/450.)**2)

      BF = 6.5E-14*FDEPEN*P*P*F*F*TH**3.6

      ABSN2 = 1.29*BF

      RETURN

      END


      
      subroutine absorb_clw(freq,temp,water,kabs_clw)
C     COMPUTES ABSORPTION IN NEPERS/KM BY SUSPENDED WATER DROPLETS
c     ARGUMENTS (INPUT):
C     WATER IN G/M**3
C     FREQ IN GHZ     (VALID FROM 0 TO 1000 GHZ)
C     TEMP IN KELVIN
      
      real  kabs_clw
C
C     REFERENCES:
C     LIEBE, HUFFORD AND MANABE, INT. J. IR & MM WAVES V.12, pp.659-675
C      (1991);  Liebe et al, AGARD Conf. Proc. 542, May 1993.
c
C     REVISION HISTORY:
C        PWR 8/3/92   original version
c        PWR 12/14/98 temp. dependence of EPS2 eliminated to agree 
c                     with MPM93 
c        pwr 2/27/02  use exponential dep. on T, eq. 2b instead of eq. 4a 
C
      COMPLEX EPS,RE
      IF(WATER.LE.0.) THEN
       kabs_clw = 0.
       RETURN
      ENDIF
      THETA1 = 1.-300./TEMP
      EPS0 = 77.66 - 103.3*THETA1
      EPS1 = .0671*EPS0
      EPS2 = 3.52                 ! from MPM93
      FP = 20.1*EXP(7.88*THETA1)  ! from eq. 2b
      FS = 39.8*FP
      EPS = (EPS0-EPS1)/CMPLX(1.,FREQ/FP) +
     & (EPS1-EPS2)/CMPLX(1.,FREQ/FS) +EPS2
      RE = (EPS-1.)/(EPS+2.)
      kabs_clw = -.06286*AIMAG(RE)*FREQ*WATER
      RETURN
      END

