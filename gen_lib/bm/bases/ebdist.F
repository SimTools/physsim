 
      SUBROUTINE EBDIST(ITYP,X1,X2,SG,FRAC)
      IMPLICIT   REAL*8 ( Z )
      INTEGER*4  ITYP
      REAL*8     X1, X2, SG, FRAC
C
C========< Entry Point >================================================
C
C--
C  Decide Gaussian smeared beam energy.
C--
C     CHI  = UGAUIN(X1)
C     FRAC = 1 + SG*CHI
C--
C  Decide flat smeared beam energy.
C--
      FRAC = 1 + SG*2*(X1-0.5D0)
C--
C  Convolute beamstrahlung.
C--
      CALL EPROBX(X2,ITYP,ZB)
      FRAC = FRAC*ZB
C--
C  That's it.
C--
      RETURN
      END
