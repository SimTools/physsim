C***********************************************************************
C*
C*==========================================================
C* Subrouine SGCXXA(GAL,GAX,GZL,GZX,GSNX,AMSN,GMSN,AMZ,GMZ,
C*                                AMX, GMX ,RS, POL, SGEFF)
C*===================================================-----==
C*
C* (Purpose)
C*   Calculates the tree level total cross section for
C*   chargino pair production.
C* (Inputs)
C*   The unit for masses and energies is GeV.
C*      GAL(2)    : (R*8) : L-A-L chiral couplings.
C*      GAX(2)    : (R*8) : XC-A-XC chiral couplings.
C*      GZL(2)    : (R*8) : L-Z-L chiral couplings.
C*      GZX(2)    : (R*8) : XC-Z-XC chiral couplings.
C*      GSNX(2,2) : (R*8) : E-SNUE-XC couplings.
C*      AMSN      : (R*8) : sneutrino mass.
C*      GMSN      : (R*8) : sneutrino width (not used).
C*      AMZ       : (R*8) : Z mass.
C*      GMZ       : (R*8) : Z width.
C*      AMX(2)    : (R*8) : XC masses.
C*      GMX(2)    : (R*8) : XC widths (not used).
C*      RS        : (R*8) : sqrt(s).
C*      POL       : (R*8) : electron beam polarization.
C* (Output)
C*      SGEFF     : (R*8) : effective cross section [pb].
C* (Relation)
C*   Invokes 
C*		none.
C* (Update Record)
C*   91/05/16  K.Fujii       Original version. 
C*
C***********************************************************************
 
      SUBROUTINE SGCXXA(GAL,GAX,GZL,GZX,GSNX,AMSN,GMSN,AMZ,GMZ,
     .                 AMX, GMX ,RS, POL, SG)
 
      IMPLICIT   REAL*8  ( A-H, O-Z )
      REAL   *8  GAL(2), GAX(2), GZL(2), GZX(2), GSNX(2,2),
     .           AMSN, GMSN, AMZ, GMZ, AMX(2), GMX(2)
      REAL   *8  RS, POL, SG(0:3)
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
C  Check if threshold crossed.
C--
      IF ( RS.LE.AMX(1)+AMX(2) ) THEN
         SG(0) = 0
         SG(1) = 0
         SG(2) = 0
         SG(3) = 0
         RETURN
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
      AMM  = 4*AMX(1)*AMX(2)/S
C--
C  Calculate coupling constants.
C--
      VE = (GZL(1)+GZL(2))/2
      AE = (GZL(2)-GZL(1))/2
      V  = (GZX(1)+GZX(2))/2
      A  = (GZX(2)-GZX(1))/2
      G1 = GSNX(1,1)/2
      G2 = GSNX(2,2)/2
      QE = (GAL(1)+GAL(2))/2
      Q  = (GAX(1)+GAX(2))/2
C--
C  Calculate amplitude squared.
C--
      DZ = S/(S-AMZ**2)
      AC = 2*(AMX(1)-AMSN)*(AMX(1)+AMSN)/S - E1/EBM + 1
      AL = LOG((AC-1+BT1)/(AC-1-BT1))
C--
      A0 = x4PI
      A2 = x4PI/3
      B0 = (x4PI/BT1)*AL
      B1 = (x4PI/BT1)*( 2*BT1*(1-AC-EE1-EE2)+(AC-1+EE1)*(AC-1+EE2)*AL )
      B2 = (x4PI/BT1)*4*( BT1
     .         +(AC-1+EE1)*(AC-1+EE2)*BT1/((AC-1-BT1)*(AC-1+BT1))
     .         +(1-AC-(EE1+EE2)/2)*AL )
C--
      RL = (1-POL)/2
      RR = (1+POL)/2
C--
      TTASSL = (QE*Q+(VE-AE)*V*DZ)**2*((AMM+EE1*EE2)*A0+BT1**2*A2)
     .       + (VE-AE)**2*A**2*DZ**2*((-AMM+EE1*EE2)*A0+BT1**2*A2)
      TTASSR = (QE*Q+(VE+AE)*V*DZ)**2*((AMM+EE1*EE2)*A0+BT1**2*A2)
     .       + (VE+AE)**2*A**2*DZ**2*((-AMM+EE1*EE2)*A0+BT1**2*A2)
      TTAST = 2*G1*G2*( (QE*Q+(VE-AE)*(V+A)*DZ)*B1
     .                + (QE*Q+(VE-AE)*(V-A)*DZ)*AMM*B0 )
      TTATT = 2*(G1*G2)**2*B2
C--
      TTASM = RL*( TTASSL + TTAST + TTATT ) + RR*TTASSR
C--
C  Calculate phase space weight.
C--
      WAT = FACT*BT1/(2*S*BT0)
C--
C  Total cross section.
C--
      SG(0) = TTASM*WAT
      SG(1) = (TTASSL*RL+TTASSR*RR)*WAT
      SG(2) = TTAST*RL*WAT
      SG(3) = TTATT*RL*WAT
C--
C  That's it.
C--
      RETURN
      END
