 
      SUBROUTINE DSGNXX(GZL,GZX,GSEX,AMSE,GMSE,AMZ,GMZ,AMX,GMX,
     .                  RS,POL,CSTH,PHI,SG)
 
      IMPLICIT   REAL*8  ( A-H, O-Y )
      IMPLICIT   COMPLEX*16 ( Z )
      REAL   *8  GZL(2)
      COMPLEX*16 GZX(2), GSEX(2,2)
      REAL   *8  AMSE(2), GMSE(2), AMZ, GMZ, AMX(2), GMX(2)
      REAL   *8  RS, POL, SG(0:6)
      REAL   *8  RPOL(2)
C--
      DATA NCALL / 0 /
C--
C  Statement function.
C--
      BETA(X1,X2) = SQRT( 1 - 2*(X1+X2) + (X1-X2)**2 )
C
C========< Entry Point >================================================
C
C--
C  Constants.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL  = 1
         NB     = 1
         xGV2PB = 3.893796623D8
         x2PI   = 2*ACOS(-1.D0)
         x4PI   = 2*x2PI
         FACT   = xGV2PB*x2PI/(x4PI**(3*NB))
         AME    = 0.511E-3
      ENDIF
C--
C  Set variables.
C--
      S    = RS*RS
      EBM  = RS/2
      PBM  = SQRT((EBM-AME)*(EBM+AME))
      BT0  = PBM/EBM
      BT1  = BETA(AMX(1)**2/S,AMX(2)**2/S)
      P1   = BT1*EBM
      E1   = SQRT(AMX(1)**2+P1**2)
      E2   = SQRT(AMX(2)**2+P1**2)
      EE1  = E1/EBM
      EE2  = E2/EBM
      T    = AMX(1)**2 - 2*EBM*E1 + 2*PBM*P1*CSTH
      U    = AMX(1)**2 - 2*EBM*E1 - 2*PBM*P1*CSTH
      AMM  = 4*AMX(1)*AMX(2)/S
C--
C  Calculate coupling constants.
C--
      VE   = (GZL(1)+GZL(2))/2
      AE   = (GZL(2)-GZL(1))/2
      ZV   = (GZX(1)+GZX(2))/2
      ZA   = (GZX(2)-GZX(1))/2
C--
C  Calculate amplitude squared.
C--
      RPOL(1) = (1-POL)/2
      RPOL(2) = (1+POL)/2
      DZ      = S/(S-AMZ**2)
C--
      TTASS = 0
      TTAST = 0
      TTASU = 0
      TTATT = 0
      TTAUU = 0
      TTATU = 0
C--
      DO 100 LR = 1, 2
         ZG1 = GSEX(1,LR)/2
         ZG2 = GSEX(2,LR)/2
         DT = S/(T-AMSE(LR)**2)
         DU = S/(U-AMSE(LR)**2)
C--
         TTASS = TTASS + RPOL(LR)*(VE-AE)**2*DZ**2
     .                  *( ABS(ZV)**2*( AMM+EE1*EE2+BT1**2*CSTH**2)
     .                   + ABS(ZA)**2*(-AMM+EE1*EE2+BT1**2*CSTH**2) )
         TTAST = TTAST + RPOL(LR)*2
     .                   *( DBLE(CONJG(ZG1)*ZG2*(ZV+ZA))*(VE-AE)*DZ*DT
     .                   *(EE1-BT1*CSTH)*(EE2-BT1*CSTH)
     .              + DBLE(CONJG(ZG1)*ZG2*(ZV-ZA))*(VE-AE)*DZ*DT*AMM )
         TTASU = TTASU - RPOL(LR)*2
     .                   *( DBLE(CONJG(ZG2)*ZG1*(ZV-ZA))*(VE-AE)*DZ*DU
     .                   *(EE1+BT1*CSTH)*(EE2+BT1*CSTH)
     .              + DBLE(CONJG(ZG2)*ZG1*(ZV+ZA))*(VE-AE)*DZ*DU*AMM )
         TTATU = TTATU - RPOL(LR)*4*DBLE((CONJG(ZG1)*ZG2)**2)*DT*DU*AMM
         TTATT = TTATT + RPOL(LR)*2*ABS(ZG1*ZG2)**2*DT**2
     .                   *(EE1-BT1*CSTH)*(EE2-BT1*CSTH)
         TTAUU = TTAUU + RPOL(LR)*2*ABS(ZG1*ZG2)**2*DU**2
     .                   *(EE1+BT1*CSTH)*(EE2+BT1*CSTH)
         AE = - AE
         ZA = - ZA
100   CONTINUE
C--
      TTASM = TTASS + TTAST + TTASU + TTATU + TTATT + TTAUU
C--
C  Calculate phase space weight.
C--
      WAT = FACT*BT1/(2*S*BT0)
C--
C  Differential cross section.
C--
      SG(0) = TTASM*WAT
      SG(1) = TTASS*WAT
      SG(2) = TTAST*WAT
      SG(3) = TTASU*WAT
      SG(4) = TTATU*WAT
      SG(5) = TTATT*WAT
      SG(6) = TTAUU*WAT
C--
C  That's it.
C--
      RETURN
      END
