 
      SUBROUTINE  GRQPTP(MODE,E,PT,DSG)
 
      IMPLICIT REAL*8  ( A-H, O-Z )
C--
C  Dummy arguments.
C     MODE < 4 : some parameters changed.
C          = 4 : nothing changed.
C--
      INTEGER*4       MODE
      REAL   *8       E, PT, DSG
C--
C  Parameters.
C--
      COMMON /TTPARM/ MT, ALFS, ALF, S2W, YCUT, RWT, GAMMAT, ALAMB,
     .                CF, MB, MZ, MW, EFFMT, TLETA1, TLETA2, KAPPA,
     .                PI, IMAG
      REAL   *8       MT, ALFS, ALF, S2W, YCUT, RWT, GAMMAT, ALAMB,
     .                CF, MB, MZ, MW, TLETA1, TLETA2, KAPPA,
     .                PI
      COMPLEX*16      EFFMT, IMAG
C--
      COMMON /QCDPAR/ ALFST, ALFSP, RQCD
      REAL*8          ALFST, ALFSP, RQCD
C--
C  Tabulated Green's functions.
C--
      COMMON /STHIRD/ PS, TILG
      REAL   *8       PS(0:100)
      COMPLEX*16      TILG(0:100)
C>>>
C     COMMON /PTHIRD/ PP, TILF
C     REAL   *8       PP(0:100)
C     COMPLEX*16      TILF(0:100)
C>>>
C--
C  Coefficients for flat and cos(theta) terms.
C--
      COMMON /COEFF/  TEVEN, TODD1, TODD2, TVTXC, TTOT
      REAL   *8       TEVEN, TODD1, TODD2, TVTXC, TTOT
C--
C  Other variables.
C--
      PARAMETER      ( NP = 100 )
      REAL   *8       DSDAT(0:NP)
C--
      REAL   *8       LNT0FI
      SAVE
C
C===========< Entry Point >=============================================
C
C--
C  Check if new ds/dp table is necessary.
C--
      IF ( MODE.LE.3 ) THEN
C--
C  Calculate S-wave contributions.
C--
         GV2PB = 3.8937966D8
         FACT  = 3*ALF**2/2/MT**4 * GV2PB
C--
         DP       = PS(1) - PS(0)
         DSDAT(0) = 0
C>>>
C        WRITE(*,'(A,11F6.3)')  'PS = ', ( PS(K), K = 0, 10 )
C        WRITE(*,'(A,11F6.3)')  'PP = ', ( PP(K), K = 0, 10 )
C        WRITE(*,'(A, 1F6.3)')  'DP = ', DP
C>>>
C--
         DO 100 IP = 1, 100
            P         = PS(IP)
            T0VC      = TEVEN * TVTXC * ABS(TILG(IP))**2
            T0FI      = TEVEN * CF*ALFSP/P * LNT0FI(IP)
            DSDAT(IP) = FACT * P**2*(T0VC+T0FI)*RUNWID(E,P)
100      CONTINUE
      ENDIF
C--
C  Calculate ds/dp, using the lookup table.
C--
      IP = PT/DP
      IF ( IP.GE.NP-1 ) THEN
         DSG = 0
      ELSE IF ( PT.LE.0.D0 ) THEN
         DSG = DSDAT(0)
      ELSE
         A1  = PT - PS(IP)
         A2  = PS(IP+1) - PT
         DSG = (A2*DSDAT(IP)+A1*DSDAT(IP+1))/DP
      END IF
C--
C  That's it.
C--
      RETURN
      END
