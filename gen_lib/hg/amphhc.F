C*
C*  GAH and GZH are defined by out-going particle charges.
C*
 
      SUBROUTINE AMPHHC(GAL,GZL,GAH,GZH,AMZ,GMZ,HIN,HOT,EIN,EOT,AMP)
 
      IMPLICIT   REAL*4  ( A-H, O-Z )
      REAL   *4  GAL(2), GZL(2), GAH, GZH, AMZ, GMZ
      COMPLEX*8  EIN(6), EOT(6), HIN(3), HOT(3), AMP(0:2)
C--
      COMPLEX*8  CURR(6)
C
C========< Entry Point >================================================
C
C--
C  Compute photon s-channel exchange diagram.
C--
      CALL JIOXXX(EIN,EOT,GAL,0.,0., CURR)
      CALL VSSXXX(CURR,HOT,HIN,GAH,AMP(1))
C--
C  Compute Z s-channel exchange diagram.
C--
      CALL JIOXXX(EIN,EOT,GZL,AMZ,GMZ, CURR)
      CALL VSSXXX(CURR,HOT,HIN,GZH,AMP(2))
C--
C  Sum two amplitudes.
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
