 
      SUBROUTINE PROBEB(X1,X2,E0,EBM,EBM0,PROB)
      IMPLICIT REAL*8 (Z)
#include "include/BEMCNS.inc"
      REAL*8     X1, X2, E0, EBM, EBM0, PROB
C>>> 93/05/08
C S300
C     DATA ITYP  / 31 /
C X300
      DATA ITYP  / 33 /
C>>>
      DATA NCALL /0/
C
C========< Entry Point >================================================
C
C--
C  Initialization.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL = 1
         FACT  = 1/SQRT(2*ACOS(-1.D0))
      ENDIF
C--
C  Decide flatly smeared beam energy.
C--
      EBM0 = E0*( 1 + SGEBM*2*(X1-0.5D0) )
      PROB = 1
C--
C  Convolute beamstrahlung.
C--
      Z    = X2
      CALL EPROBX(Z,ITYP,ZB)
      EBM  = EBM0*ZB
C--
C  Calculate corresponding probability density.
C--
C     PROB = PROB*EPROB(ZB,ZE)
C--
C  That's it.
C--
      RETURN
      END
