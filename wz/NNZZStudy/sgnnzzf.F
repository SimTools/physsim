CC**********************************************************************
C* 
C*==============================-----===
C* Subroutine SGNNZZF(IDP,IHEL,X,SG,PV)
C*==============================-----===
C*
C* (Purpose)
C*    Calculate differential cross section for e+ e- ---> nb n Z Z.
C*    with the Z's subsequent decays.
C* (Inputs)
C*       IDP(1,i) : (I*4) : generation number.
C*          (2,i) : (I*4) : T3 + 1.5.
C*          (3,i) : (I*4) : l/q flag; (1,2) = (l,q).
C*       IHEL( i) : (I*4) : helicity combination.
C*                        : i = (1,2,3,4,5,6,7,8)
C*                        :   = (e-,e+, nu,nub, fd,fub, fdb,fu)
C*                        :                   <- Z1 -><- Z2 ->
C*       X(1,1)   : (R*8) : s.
C*        (2,1)   : (R*8) : eta_1.
C*        (3,1)   : (R*8) : phi_1.
C*       X(1,2)   : (R*8) : log(1-E_i/E_bm).
C*        (2,2)   : (R*8) : eta_2
C*        (3,2)   : (R*8) : phi_2.
C*       X(1,3)   : (R*8) : invariant mass squared for Z1Z2.
C*        (2,3)   : (R*8) : cos_Z1 in (Z1Z2) rest frame.
C*        (3,3)   : (R*8) : phi_Z1 in (Z1Z2) rest frame.
C*       X(1,4)   : (R*8) : invariant mass squared for Z1.
C*        (2,4)   : (R*8) : cos_fd  in Z1 rest frame.
C*        (3,4)   : (R*8) : phi_fub in Z1 rest frame.
C*       X(1,5)   : (R*8) : invariant mass squared for Z2.
C*        (2,5)   : (R*8) : cos_fdb in Z2 rest frame.
C*        (3,5)   : (R*8) : phi_fu  in Z2 rest frame.
C* (Output)
C*       PV(*,i)  : (R*8) : 4-momentum of i-th particle, where
C*                        : numbering convention is that of IDP.
C*       SG       : (R*8) : d(sigma)/dX1dX2....dXn.
C* (Relation)
C*    Invokes UVZERO, UHSETF, UBTRAN, FULNNZZ, and routines
C*    in the HELAS library.
C* (Update Record)
C*    2005/03/06  K.Fujii	Original JSF version.
C*
CC**********************************************************************
 
      SUBROUTINE SGNNZZF(IDP,IHEL,X,SG,PV)
 
      IMPLICIT     REAL*8  ( A-H, O-Z )
#include "SMCONS.inc"
#include "SMPTAB.inc"
      PARAMETER    ( NP = 8, NB = NP - 3 )
      INTEGER*4    IDP(3,NP), IHEL(NP)
      REAL   *8    X(3,NB), SG, PV(0:3,NP)
C--
      REAL   *8    QV(0:3,NP), EC(3,3), EB(3,3), EA(3,3), 
     .             XX(3,2), AMP2(0:3)
      DATA NCALL /0/
C
C========< Entry Point >================================================
C
C--
C  Constants.
C  Spin average is taken only for positron.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL  = 1
         SPIN   = 2
         FACT   = xGV2FB*x2PI/(x4PI**(3*NB))/SPIN
C        FACT   = FACT/(xALF*137.0359895D0)**2
         IDNT   = 2
         FACT   = FACT/IDNT
         AME    = AMF(1,2,1)
         CALL UVZERO(18,EC)
         EC(1,1) = 1
         EC(2,2) = 1
         EC(3,3) = 1
      ENDIF
C--
C  Reset phase space weight.
C--
      WAT = 1
C--
C  Define some kinematical variables.
C--
      S       = X(1,1)
      RS      = SQRT(S)
      EBM     = RS/2
C--
      AM3     = AMF(IDP(1,3),IDP(2,3),IDP(3,3))
      AM4     = AMF(IDP(1,4),IDP(2,4),IDP(3,4))
      Q3      = SQRT(X(1,3))
      Q4      = SQRT(X(1,4))
      Q5      = SQRT(X(1,5))
C--      
      RXI     = EXP(X(1,2))
      WAT     = WAT*RXI
C-- 
C For gamma gamma     
C      DM1     = (AME*AME/S)*RXI*RXI/(1-RXI)
C      DP1     = (AME*AME/S)
C For W W fusion
C
      DM1     = AMW*AMW/S
      DP1     = 1
      DM2     = DP1
      DP2     = DM1
C--
      ET1     = X(2,1)
      EXP2ET  = EXP(2*ET1)
      SH1     = MIN( SQRT( (1+DP1+DM1)/(1+EXP2ET)   - DM1 ), 1.D0 )
      CH1     = MIN( SQRT( (1+DP1+DM1)/(1+1/EXP2ET) - DP1 ), 1.D0 )
      SN1     = 2*SH1*CH1
      CS1     = (1+DP1+DM1)*TANH(ET1) - DP1 + DM1
      FI1     = X(3,1)
      WAT     = WAT*4*(1+DP1+DM1)/((1+EXP2ET)*(1+1/EXP2ET))
C--
      ET2     = X(2,2)
      EXP2ET  = EXP(2*ET2)
      SH2     = MIN( SQRT( (1+DP2+DM2)/(1+EXP2ET)   - DM2 ), 1.D0 )
      CH2     = MIN( SQRT( (1+DP2+DM2)/(1+1/EXP2ET) - DP2 ), 1.D0 )
      SN2     = 2*SH2*CH2
      CS2     = (1+DP2+DM2)*TANH(ET2) - DP2 + DM2
      FI2     = X(3,2)
      WAT     = WAT*4*(1+DP2+DM2)/((1+EXP2ET)*(1+1/EXP2ET))
C--
      CS12    = CS1*CS2 + SN1*SN2*COS(FI2-FI1)
C--
      IF ( ET1.GT.-ET2 ) THEN
         RXIMN   = (Q3-AM3+AM4)*(Q3+AM3+AM4)/S
         RXIMX   = 1 - 2*AM3/RS
         X1      = 1 - RXI
         X2      = ( RXI - X(1,3)/S )/( 1 - X1*(1-CS12)/2 )
         IF ( 1-X2.LT.RXIMN .OR. 1-X2.GT.RXIMX )	GO TO 9999
         WAT     = WAT*S*X1*X2*X2/( RXI - X(1,3)/S )
      ELSE
         RXIMN   = (Q3+AM3-AM4)*(Q3+AM3+AM4)/S
         RXIMX   = 1 - 2*AM4/RS
         X2      = 1 - RXI
         X1      = ( RXI - X(1,3)/S )/( 1 - X2*(1-CS12)/2 )
         IF ( 1-X1.LT.RXIMN .OR. 1-X1.GT.RXIMX )	GO TO 9999
         WAT     = WAT*S*X2*X1*X1/( RXI - X(1,3)/S )
      ENDIF
C--
C  Set 4-momenta.
C--
      CALL UVZERO(8*NP,PV)
      CALL UVZERO(8*NP,QV)
C
C-- e+ + e- --> CM.
C
      PBM     = SQRT((EBM-AME)*(EBM+AME))
      BT0     = PBM/EBM
      PV(0,1) = EBM
      PV(3,1) = PBM
      CALL PMIRRx(PV(0,1),PV(0,2))
      QV(0,1) = RS
C
C-- CM --> a + b + X.
C
C-(a)
      E1      = EBM*X1
      IF ( E1.LT.AM3 )				 GO TO 9999
      APV     = SQRT((E1-AM3)*(E1+AM3))
      PV(0,3) = E1
      PV(1,3) = APV*SN1*COS(FI1)
      PV(2,3) = APV*SN1*SIN(FI1)
      PV(3,3) = APV*CS1
C-(b)
      E2      = EBM*X2
      IF ( E2.LT.AM4 ) 				 GO TO 9999
      APV     = SQRT((E2-AM4)*(E2+AM4))
      PV(0,4) = E2
      PV(1,4) = APV*SN2*COS(FI2)
      PV(2,4) = APV*SN2*SIN(FI2)
      PV(3,4) = APV*CS2
C-(X)
      CALL PSUBxx(QV(0,1),PV(0,3),QV(0,3))
      CALL PSUBxx(QV(0,3),PV(0,4),QV(0,3))
C
C-- X --> C + D.
C
      AM12 = X(1,4)
      AM22 = X(1,5)
      COS1 = X(2,3)
      PHI1 = X(3,3)
      CALL UHPHS2(2,QV(0,3),AM12,AM22,EC,COS1,PHI1,
     .            QV(0,4),QV(0,5),BT3,EB)
      IF ( BT3.EQ.0.D0 ) 				GO TO 9999
C
C-- C --> f + fbar.
C
      AM12 = AMF(IDP(1,5),IDP(2,5),IDP(3,5))**2
      AM22 = AMF(IDP(1,6),IDP(2,6),IDP(3,6))**2
      COS1 = X(2,4)
      PHI1 = X(3,4)
      CALL UHPHS2(2,QV(0,4),AM12,AM22,EB,COS1,PHI1,
     .            PV(0,5),PV(0,6),BT4,EA)
      IF ( BT4.EQ.0.D0 )                           	GO TO 9999
C
C-- D --> f' + f'bar.
C
      AM12 = AMF(IDP(1,7),IDP(2,7),IDP(3,7))**2
      AM22 = AMF(IDP(1,8),IDP(2,8),IDP(3,8))**2
      COS1 = X(2,5)
      PHI1 = X(3,5)
      CALL UHPHS2(2,QV(0,5),AM12,AM22,EB,COS1,PHI1,
     .            PV(0,7),PV(0,8),BT5,EA)
      IF ( BT5.EQ.0.D0 )                           	GO TO 9999
C--
      XX(1,1) = SH1
      XX(2,1) = CH1
      XX(3,1) = FI1
      XX(1,2) = SH2
      XX(2,2) = CH2
      XX(3,2) = FI2
C--
C  Calculate phase space weight.
C--
      WAT = WAT*FACT*BT3*BT4*BT5/(2*S*BT0)
C--
C  Calculate amplitude squared.
C--
C>>>
CDBG       CALL PSUMxx(PV(0,3),PV(0,4),QV(0,6))
CDBG       CALL PSUMxx(QV(0,6),PV(0,5),QV(0,6))
CDBG       CALL PSUMxx(QV(0,6),PV(0,6),QV(0,6))
CDBG       CALL PSUMxx(QV(0,6),PV(0,7),QV(0,6))
CDBG       CALL PSUMxx(QV(0,6),PV(0,8),QV(0,6))
CDBG       PRINT *, ' '
CDBG       PRINT *, ' --- '
CDBG       PRINT *, ' ET1, SN1, CS1 = ', ET1, SN1, CS1
CDBG       PRINT *, ' ET2, SN2, CS2 = ', ET2, SN2, CS2
CDBG       PRINT *, ' SH1, CH1, FI1 = ', (XX(K,1),K=1,3)
CDBG       PRINT *, ' SH2, CH2, FI2 = ', (XX(K,2),K=1,3)
CDBG       PRINT *, ' '
CDBG       PRINT *, ' IDP  = ', IDP
CDBG       PRINT *, ' IHEL = ', (IHEL(K),K=1,8)
CDBG       PRINT *, ' PV(*,1) = ', (PV(K,1),K=0,3)
CDBG       PRINT *, '   (*,2) = ', (PV(K,2),K=0,3)
CDBG       PRINT *, '   (*,3) = ', (PV(K,3),K=0,3)
CDBG       PRINT *, '   (*,4) = ', (PV(K,4),K=0,3)
CDBG       PRINT *, '   (*,5) = ', (PV(K,5),K=0,3)
CDBG       PRINT *, '   (*,6) = ', (PV(K,6),K=0,3)
CDBG       PRINT *, '   (*,7) = ', (PV(K,7),K=0,3)
CDBG       PRINT *, '   (*,8) = ', (PV(K,8),K=0,3)
CDBG       PRINT *, ' '
CDBG       PRINT *, ' PV_SUM  = ', (QV(K,6),K=0,3)
CDBG       PRINT *, ' '
C>>>
      CALL FULNNZZ(IDP,IHEL,XX,PV,AMP2)
C--
C  Differenctial cross section.
C--
      SG  = AMP2(0)*WAT
CNF       SG = AMP2(1)*WAT
CVV       SG = AMP2(2)*WAT
CHG       SG = AMP2(3)*WAT
C>>>
CDBG       PRINT *, ' '
CDBG       PRINT *, ' SG, WAT = ', SG, WAT
CDBG       PRINT *, ' AMP2    = ', (AMP2(K),K=0,3)
CDBG       PRINT *, ' '
C>>>
C--
C  That's it.
C--
      RETURN
C--
C  Kinematically forbidden.
C--
9999  SG = 0
      RETURN
      END
