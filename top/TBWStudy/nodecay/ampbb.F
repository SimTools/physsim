C*(Update Record)
C*  2010/11/29  K.Fujii		original version.
C*

      SUBROUTINE AMPBB(GAL,GAD,GZL,GZD,AMZ,GMZ,
     .                 EIN,EOT,BIN,BOT,AMP)
 
      IMPLICIT   REAL*8  ( A-H, O-Z )
      REAL   *8  GAL(2), GAD(2), GZL(2), GZD(2), AMZ, GMZ
      COMPLEX*16 EIN(6), EOT(6), BIN(6), BOT(6), AMP(0:2)
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
      CALL IOVXXX(BIN,BOT,CURR,GAD,AMP(1))
C--
C  Compute Z s-channel exchange diagram.
C--
      CALL JIOXXX(EIN,EOT,GZL,AMZ,GMZ, CURR)
      CALL IOVXXX(BIN,BOT,CURR,GZD,AMP(2))
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
