C **********************************************************************
C
C   SUBROUTINE EWINIT
C(FUNCTION)
C   INITIALIZE CONSTANTS FOR EW THEORY
C(OUTPUT)
C   TO COMMON /EWCONS/
C(UPDATE RECORD)
C   95/03/26  K.Fujii  	Modified to get constants from dummy
C			arguments.
C
C **********************************************************************

      SUBROUTINE EWINIT(GV2PB,ALF0,ALF,S2W,AMZ,GMZ,GMW,AMH,AMT)

      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /EWCONS/ EWPARA
      REAL*8          EWPARA(50)
C
C========< Entry Point >================================================
C
C-- 
C  Basic Constants.
C--
      EWPARA( 1) = DACOS(-1.D0)
      EWPARA( 2) = 2*EWPARA(1)
      EWPARA( 3) = GV2PB
      EWPARA( 4) = ALF0
      EWPARA( 5) = ALF
      EWPARA( 6) = S2W
      EWPARA( 7) = 1.D0 - EWPARA(6)
      EWPARA( 8) = DSQRT(EWPARA(6))
      EWPARA( 9) = DSQRT(EWPARA(7))
C--
C  Masses: order of 15 and 16 for AMH and AMT are KF's guess.
C  Not used in MHIGGS and does not matter.
C-- 
      EWPARA(11) = AMZ
      EWPARA(12) = GMZ
      EWPARA(13) = DSQRT(EWPARA(7)*EWPARA(11)**2)
      EWPARA(14) = GMW
      EWPARA(15) = AMH
      EWPARA(16) = AMT
C--
C  Couplings.
C-- 
      EWPARA(21) = 4.D0*EWPARA(1)*EWPARA( 5)
      EWPARA(22) = DSQRT(EWPARA(21))
      EWPARA(23) = EWPARA(22)/EWPARA(8)
      EWPARA(24) = EWPARA(23)/EWPARA(9)
C--
      EWPARA(31) =  0.25D0
      EWPARA(32) = -0.25D0
      EWPARA(33) = -0.25D0 + EWPARA( 6)
      EWPARA(34) =  0.25D0
      EWPARA(35) =  0.25 - 2./3.*EWPARA( 6)
      EWPARA(36) = -0.25D0
      EWPARA(37) = -0.25D0 + 1./3.*EWPARA( 6)
      EWPARA(38) = 0.25D0
C--
C  That's it.
C-- 
      RETURN
      END
