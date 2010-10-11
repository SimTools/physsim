      SUBROUTINE AMPNXX (GZL,GZX,GSEX,AMSE,GMSE,AMZ,GMZ,
     .                   EIN,EOT, XIIN,XJIN, XIOT,XJOT, AMPS)
 
      IMPLICIT   REAL*8  (A-H, O-Z)
      REAL   *8  GZL(2),
     .           AMSE(2), GMSE(2), AMZ, GMZ
C--
C  Complex couplings
C     GZX (*)           : X0_i-X0_j-Z
C     GSEX(*,IJ,LR,IOX) : IOX = (1,2) 
C                       :     = (se_LR-e-X0_[IJ], X0_[IJ]-e-se_LR)
C                       : IJ  = (1,2) = (i,j)
C--
      COMPLEX*16 GZX(2), GSEX(2,2,2,2)
      COMPLEX*16 EIN(6), EOT(6), XIIN(6), XJIN(6), XIOT(6), XJOT(6),
     .           AMPS(0:3)
C--
      COMPLEX*16 AMP(5), CURR(6), SCAL(3)
C
C========< Entry Point >================================================
C
C--
C (1) s-channel Z exchange.
C--
      CALL JIOXXX(EIN,EOT, GZL, AMZ,GMZ, CURR)
      CALL IOVXXC(XIIN,XJOT,CURR, GZX, AMP(1))
C--
C (2) t-channel seL exchange.
C--
      LR = 1
      CALL HIOXXX(EIN,XJOT, GSEX(1,2,LR,2), AMSE(LR),GMSE(LR), SCAL)
      CALL IOSXXX(XIIN,EOT,SCAL, GSEX(1,1,LR,1), AMP(2))
      AMP(2) = -AMP(2)
C--
C (3) u-channel seL exchange.
C--
      CALL HIOXXX(EIN,XIOT, GSEX(1,1,LR,2), AMSE(LR),GMSE(LR), SCAL)
      CALL IOSXXX(XJIN,EOT,SCAL, GSEX(1,2,LR,1), AMP(3))
C--
C (4) t-channel seR exchange.
C--
      LR = 2
      CALL HIOXXX(EIN,XJOT, GSEX(1,2,LR,2), AMSE(LR),GMSE(LR), SCAL)
      CALL IOSXXX(XIIN,EOT,SCAL, GSEX(1,1,LR,1), AMP(4))
      AMP(4) = -AMP(4)
C--
C (5) u-channel seR exchange.
C--
      CALL HIOXXX(EIN,XIOT, GSEX(1,1,LR,2), AMSE(LR),GMSE(LR), SCAL)
      CALL IOSXXX(XJIN,EOT,SCAL, GSEX(1,2,LR,1), AMP(5))
C--
C  Sum three amplitudes.
C     AMPS(0) = sum.
C         (1) = s-channel Z.
C         (2) = t/u-channel seL.
C         (3) = t/u-channel seR.
C--
      AMPS(0) = AMP(1) + AMP(2) + AMP(3) + AMP(4) + AMP(5)
      AMPS(1) = AMP(1)
      AMPS(2) = AMP(2) + AMP(3)
      AMPS(3) = AMP(4) + AMP(5)
C--
C  That's it.
C--
      RETURN
      END
