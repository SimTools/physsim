C*  95/06/23  K.Fujii		Original version for BASES 5.1.
C*
      REAL FUNCTION FUNC*8(ZZ)

      IMPLICIT REAL*8 ( A-H, O-Z )
      REAL*8   ZZ(50), Z(50)
C--
C  Generator COMMON's.
C--
      INCLUDE 'BEMCNF.inc'
#include "BSHUFL.inc"
      DATA NCALL / 0 /
C
C========< Entry Point >================================================
C
C--
C  Set some variables.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL = 1
C--
C  Machine selection.
C        ITYP  = 31 : JLC-I S 300
C              = 32 : JLC-I C 300
C              = 33 : JLC-I X 300
C              = 34 : JLC-I S 500
C              = 35 : JLC-I C 500
C              = 36 : JLC-I X 500
C        SGEB  = 0.007 : JLC-I design
C              = 0.005 : tt nominal
C	       = 0.002 : controlled
C--
         ISRB  = ISRBM
         ITYP  = IBMTP
         SGEB  = SGEBM
         SQRTS = ROOTS
C--
         xALF  = 1/137.D0
         xPI   = ACOS(-1.D0)
         AME   = 0.511D-3
      ENDIF
C--
C  Shuffle integration variables.
C--
      DO 10 I = 1, NZZ
         Z(ISHUFL(I)) = ZZ(I)
10    CONTINUE 
C--
C  Set independent variables.
C     Z( 1) : Energy spread of e- beam.
C      ( 2) : Energy spread of e+ beam.
C      ( 3) : Beamstrahlung of e- beam.
C      ( 4) : Beamstrahlung of e+ beam.
C      ( 5) : Bremsstralung.
C--
      IF ( ISRB/2.EQ.1 ) THEN
         CALL EBDIST(ITYP,Z(1),Z(3),SGEB,FM)
         CALL EBDIST(ITYP,Z(2),Z(4),SGEB,FP)
         FRAC = SQRT(FM*FP)
      ELSE
         FRAC = 1
         FM   = 1
         FP   = 1
      ENDIF
C--
C  Then decide reduced sqrt(s) after bremsstrahlung.
C--
      FUNC = 1
      IF ( MOD(ISRB,2).EQ.1 ) THEN
         RS   = SQRTS*FRAC
         BTE  = (2*xALF/xPI)*(2*LOG(MAX(RS,2*AME)/AME)-1)
         ZG   = Z(5)**(1/BTE)
         FRAC = SQRT(1-ZG)*FRAC
         IF ( ZG .NE. 0.0D0 ) THEN
            FUNC = ( BTE*ZG**(BTE-1)*(1+3*BTE/4) - BTE*(1-ZG/2) )
            FUNC = FUNC/(BTE*ZG**(BTE-1)*(1+3*BTE/4))
         ENDIF
      ENDIF
C--
      XM = FM
      XP = FP
      FS = FRAC
C--
C  Fill hists. and plots.
C--
      FRAC = MAX(1.D-6,1+5*SGEB-FRAC)
      FRAC = -LOG(FRAC)
      IF ( FRAC.LT.3.D0 ) THEN
         CALL XHFILL(1,FRAC,FUNC)
      ELSE IF ( FRAC.LT.4.5D0 ) THEN
         CALL XHFILL(2,FRAC,FUNC)
      ELSE IF ( FRAC.LT.6.0D0 ) THEN
         CALL XHFILL(3,FRAC,FUNC)
      ELSE
         CALL XHFILL(4,FRAC,FUNC)
      ENDIF
      CALL XHFILL(5,DBLE(FS),FUNC)
C--
C  That's it.
C--
      RETURN
      END
