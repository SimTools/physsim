      SUBROUTINE  SETPRM
 
      IMPLICIT REAL*8 ( A-H, O-Z )
C--
      COMMON /PARAM1/ E, GAMMA0, ALF, ALFS, S2W, AMT, AMW, AMZ, GMZ,
     .                AMB, VTB, AMH, BTH, ALFST
      REAL   *8       E, GAMMA0, ALF, ALFS, S2W, AMT, AMW, AMZ, GMZ,
     .                AMB, VTB, AMH, BTH, ALFST
      COMMON /EFFMAS/ ZMU, BETA
      COMPLEX*16      ZMU, BETA
C--
      COMMON /KFFLAG/ IPOT, IRGM
      INTEGER*4       IPOT, IRGM
C--
      REAL   *8       CAPPA, C1, C2, GF, W, PI
      COMPLEX*16      ZI
      PARAMETER ( ZI = (0.,1.) )
      DATA NCALL /0/
C
C========< Entry Point >================================================
C
C--
C  Initialize constants.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL = 1
         PI = ACOS(-1.D0)
         GF = (PI/SQRT(2.D0))*(ALF/S2W)/AMW**2
      ENDIF
C--
C  Calculate on-shell width and effective mass.
C--
      W   = AMW/AMT
      C10 = 1 - 3*W**4 + 2*W**6
      IF ( IRGM.EQ.1 ) THEN
C>>>
C 94/08/01 Modified according to TUW-94-06.
C
CO       C1E = ( 4 - 4*W**6 )*(E/AMT)
C1       C1E = ( 3 + 3*W**4 - 6*W**6 )*(E/AMT)
         C1E = 0
C>>>
         C1  = C10 + C1E
C>>>
CVO      C2  = 23.D0/6 + 0.5D0*W**4 - 13.D0/3*W**6
CO       C2  = 5.D0/3 + 7*W**4 - 26.D0/3*W**6
C1       C2  = 13.D0/6 + 11.D0/2*W**4 - 23.D0/3*W**6
         C2  = C10/2
C>>>
         CAPPA  = VTB**2*GF/SQRT(2.D0) * AMT**2/8/PI * C2
      ELSE
         C1     = C10
         CAPPA  = 0
      END IF
      GAMMA0 = VTB**2*GF/SQRT(2.D0) * AMT**3/4/PI * C1
      ZMU    = AMT/( 1 + ZI*CAPPA )
C--
      IF ( IPOT.EQ.1 .OR. IPOT.EQ. 2 ) THEN
         BETA = (8*PI/23) * ZMU
      ELSE IF ( IPOT.EQ.3 ) THEN
         BETA = 4.D0/3 * ZMU * ALFS
      ELSE
         BETA = 0
      END IF
C--
C  That's it.
C--
      RETURN
      END
