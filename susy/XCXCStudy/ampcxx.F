      SUBROUTINE AMPCXX (GAL,GAX,GZL,GZX,GSNX,AMSN,GMSN,AMZ,GMZ,
     .                   EIN,EOT,XIN,XOT,AMP)
 
      IMPLICIT   REAL*4  ( A-H, O-Z )
      REAL   *4  GAL(2), GAX(2), GZL(2), GZX(2),
     .           AMSN, GMSN, AMZ, GMZ
      COMPLEX*8  GSNX(2,2)
      COMPLEX*8  EIN(6), EOT(6), XIN(6), XOT(6), AMP(0:3)
C--
      COMPLEX*8  CURR(6), SCAL(3)
C
C========< Entry Point >================================================
C
C--
C  Compute photon s-channel exchange diagram.
C--
      CALL JIOXXX(EIN,EOT,GAL,0.,0.,CURR)
      CALL IOVXXX(XIN,XOT,CURR,GAX,AMP(1))
C--
C  Compute Z s-channel exchange diagram.
C--
      CALL JIOXXX(EIN,EOT,GZL,AMZ,GMZ,CURR)
      CALL IOVXXX(XIN,XOT,CURR,GZX,AMP(2))
C--
C  Compute t-channel snu exchange diagram.
C--
      CALL HIOXXX(EIN,XOT,GSNX(1,1),AMSN,GMSN,SCAL)
      CALL IOSXXX(XIN,EOT,SCAL,GSNX(1,2),AMP(3))
      AMP(3) = - AMP(3)
C--
C  Sum three amplitudes.
C     AMP(0) = sum.
C        (1) = photon-exchange.
C        (2) = Z-exchange.
C        (3) = snu-exchange.
C--
      AMP(0) = AMP(1) + AMP(2) + AMP(3)
C--
C  That's it.
C--
      RETURN
      END
