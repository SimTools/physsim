C*(Update Record)
C*  2010/11/29  K.Fujii		original version.
C*

      SUBROUTINE AMPTT(GAL,GAU,GZL,GZU,AMZ,GMZ,
     .                 EIN,EOT,TIN,TOT,AMP)
 
      IMPLICIT   REAL*8  ( A-H, O-Z )
      REAL   *8  GAL(2), GAU(2), GZL(2), GZU(2), AMZ, GMZ
      COMPLEX*16 EIN(6), EOT(6), TIN(6), TOT(6), AMP(0:2)
C--
      REAL   *8  VZT(2), AZT(2)
      COMPLEX*16 CURR(6)
C
C========< Entry Point >================================================
C
C--
C  Compute photon s-channel exchange diagram.
C--
      CALL JIOXXX(EIN,EOT,GAL,0.D0,0.D0, CURR)
      CALL IOVXXX(TIN,TOT,CURR,GAU,AMP(1))
C--
C  Compute Z s-channel exchange diagram.
C--
      CALL JIOXXX(EIN,EOT,GZL,AMZ,GMZ, CURR)
      CALL IOVXXX(TIN,TOT,CURR,GZU,AMP(2))
C--
C  Sum two amplitudes.
C     AMP(0) = sum.
C        (1) = phton.
C        (2) = Z.
C--
      AMP(0) = AMP(1) + AMP(2)
C--
C  That's it.
C--
      RETURN
      END
