      SUBROUTINE AMPSSE (GAL,GAX,GZL,GZX,AMZ,GMZ,AMN,GMN,GNESE,
     .                   EIN,EOT,XIN,XOT,AMP)
 
      IMPLICIT   REAL*4  ( A-H, O-Z )
      REAL   *4  GAL(2), GAX, GZL(2), GZX, AMZ, GMZ, AMN(4), GMN(4)
      COMPLEX*8  GNESE(2,2,4)
      COMPLEX*8  EIN(6), EOT(6), XIN(3), XOT(3), AMP(0:6)
C--
      COMPLEX*8  CURR(6), FSI(6)
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
C  Compute t-channel neutralino exchange diagram.
C--
      DO 10 IX = 1, 4
         CALL FSIXXX(EIN,XOT,GNESE(1,1,IX),AMN(IX),GMN(IX),FSI)
         CALL IOSXXX(FSI,EOT,XIN,GNESE(1,2,IX),AMP(2+IX))
10    CONTINUE
C--
C  Sum three amplitudes.
C     AMP(0  ) = sum.
C        (1  ) = photon-exchange.
C        (2  ) = Z-exchange.
C        (3-6) = t-channel neutralino exchanges.
C--
      AMP(0) = AMP(1) + AMP(2) + AMP(3) + AMP(4) + AMP(5) + AMP(6)
C--
C  That's it.
C--
      RETURN
      END
