C*(Update Record)
C*  95/03/16  K.Fujii		New version with HELAS V204.
C*

      SUBROUTINE AMPTT (GAL,GAU,GZL,GZU,AMZ,GMZ,
     .                  EIN,EOT,TIN,TOT,AMP)
 
      IMPLICIT   REAL*4  ( A-H, O-Z )
      REAL   *4  GAL(2), GAU(2), GZL(2), GZU(2), AMZ, GMZ
      COMPLEX*8  EIN(6), EOT(6), TIN(6), TOT(6), AMP(0:2)
C--
      COMPLEX*8  CURR(6)
C
C========< Entry Point >================================================
C
C--
C  Compute photon s-channel exchange diagram.
C--
      CALL JIOXXX(EIN,EOT,GAL,0.,0., CURR)
      CALL IOVXXX(TIN,TOT,CURR,GAU,AMP(1))
C--
C  Compute Z s-channel exchange diagram.
C--
      CALL JIOXXX(EIN,EOT,GZL,AMZ,GMZ, CURR)
      CALL IOVXXX(TIN,TOT,CURR,GZU,AMP(2))
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
