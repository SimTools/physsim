      SUBROUTINE AMPSSM (GAL,GAX,GZL,GZX,AMZ,GMZ,
     .                   EIN,EOT,XIN,XOT,AMP)
 
      IMPLICIT   REAL*4  ( A-H, O-Z )
      REAL   *4  GAL(2), GAX, GZL(2), GZX, AMZ, GMZ
      COMPLEX*8  EIN(6), EOT(6), XIN(3), XOT(3), AMP(0:2)
C--
      COMPLEX*8  CURR(6)
      DATA AMG, GMG   / 0., 0. /
C
C========< Entry Point >================================================
C
C--
C  Compute photon s-channel exchange diagram.
C--
      CALL JIOXXX(EIN,EOT,GAL,AMG,GMG,CURR)
      CALL VSSXXX(CURR,XOT,XIN,GAX,AMP(1))
C--
C  Compute Z s-channel exchange diagram.
C--
      CALL JIOXXX(EIN,EOT,GZL,AMZ,GMZ,CURR)
      CALL VSSXXX(CURR,XOT,XIN,GZX,AMP(2))
C--
C  Sum three amplitudes.
C     AMP(0) = sum.
C        (1) = photon-exchange.
C        (2) = Z-exchange.
C--
      AMP(0) = AMP(1) + AMP(2)
C--
C  That's it.
C--
      RETURN
      END
