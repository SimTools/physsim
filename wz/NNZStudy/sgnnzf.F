CC**********************************************************************
C*
C*==============================-----===
C* Subroutine SGNNZF(IDP,IHEL,X,SG,PV)
C*==============================-----===
C*
C* (Purpose)
C*    Calculate differential cross section for e+ + e- ---> nu nub + Z.
C*    with the Z's subsequent decays.
C* (Inputs)
C*       IDP(1,i) : (I*4) : generation number.
C*          (2,i) : (I*4) : T3 + 1.5.
C*          (3,i) : (I*4) : l/q flag; (1,2) = (l,q).
C*       IHEL( i) : (I*4) : helicity combination.
C*                        : i = (1,2,3,4,5,6,7,8)
C*                        :   = (e-,e+,nu,nub,fb,f)
C*                        :                   <-Z->
C*       X(1,1)   : (R*8) : s.
C*        (2,1)   : (R*8) : eta_nu.
C*        (3,1)   : (R*8) : phi_nu.
C*       X(1,2)   : (R*8) : log(1-E_nu/E_bm).
C*        (2,2)   : (R*8) : zeta_nub.
C*        (3,2)   : (R*8) : phi_nub.
C*       X(1,3)   : (R*8) : invariant mass squared for Z.
C*        (2,3)   : (R*8) : cos_fb in Z rest frame.
C*        (3,3)   : (R*8) : phi_f  in Z rest frame.
C* (Output)
C*       PV(*,i)  : (R*8) : 4-momentum of i-th particle, where
C*                        : numbering convention is that of IDP.
C*       SG       : (R*8) : d(sigma)/dX1dX2....dXn.
C* (Relation)
C*    Invokes UVZERO, UHSETF, UBTRAN, FULNNZ, and routines
C*    in the HELAS library.
C* (Update Record)
C*    95/04/26  K.Fujii		Original version.
C*
CC**********************************************************************
 
      SUBROUTINE SGNNZF(IDP,IHEL,X,SG,PV)
 
      IMPLICIT     REAL*8  ( A-H, O-Z )
#include "SMCONS.inc"
#include "SMPTAB.inc"
      INTEGER*4    IDP(3,*), IHEL(*)
      REAL   *8    X(3,*), SG, PV(0:3,*)
C--
      REAL   *8    QV(0:3,6), EC(3,3), EB(3,3), EA(3,3), AMP2(0:3)
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
         NB     = 3
         FACT   = xGV2FB*x2PI/(x4PI**(3*NB))/SPIN
         AME    = AMF(1,2,1)
         CALL UVZERO(2*9,EC)
         EC(1,1) = 1
         EC(2,2) = 1
         EC(3,3) = 1
      ENDIF
C--
C  Set 4-momenta.
C--
      CALL UVZERO(2*24,PV)
      CALL UVZERO(2*24,QV)
C--
      WAT = 1
C
C-- e+ + e- --> CM.
C
      S       = X(1,1)
      RS      = SQRT(S)
      EBM     = RS/2
      PBM     = SQRT((EBM-AME)*(EBM+AME))
      BT0     = PBM/EBM
      PV(0,1) = EBM
      PV(3,1) = PBM
      CALL PMIRRx(PV(0,1),PV(0,2))
      QV(0,1) = RS
C--
      AM3     = AMF(IDP(1,3),IDP(2,3),IDP(3,3))
C--
      EXPXI   = EXP(X(1,2))
      Q232    = S*EXPXI + AM3*AM3
      Q23     = SQRT(Q232)
      WAT     = WAT*S*EXPXI
C--
C>>>
      DM1     = AMW*AMW/S
      DP1     = 1
C>>>
C--
      EXP2ET  = EXP(2*X(2,1))
      SH1     = MIN( SQRT( (1+DP1+DM1)/(1+EXP2ET) - DM1 ), 1.D0 )
      CH1     = SQRT( (1-SH1)*(1+SH1) )
      CS1     = ( 1 - xSQRT2*SH1 )*( 1 + xSQRT2*SH1 )
      FI1     = X(3,1)
      WAT     = WAT*4*(1+DP1+DM1)/((1+EXP2ET)*(1+1/EXP2ET))
C--
      Q32     = X(1,3)
      Q3      = SQRT(Q32)
      IF ( Q3.GT.0.D0 .AND. Q23.GT.Q3 ) THEN
         DLT  = 2*Q32/((Q23-Q3)*(Q23+Q3))
      ELSE
         DLT  = 10
      ENDIF
C--      
      OPDZT   = ( 1 + 1/DLT )**X(2,2)
      SH2     = MIN( SQRT( DLT*( OPDZT - 1 ) ), 1.D0 )
      CS2     = ( 1 - xSQRT2*SH2 )*( 1 + xSQRT2*SH2 )
      FI2     = X(3,2)
      WAT     = WAT*2*DLT*OPDZT*LOG(1+1/DLT)
C
C-- CM --> a + (b C).
C
      AM12 = AM3*AM3
      AM22 = Q232
      CALL UHPHS2(1,QV(0,1),AM12,AM22,EC,CS1,FI1,
     .            PV(0,3),QV(0,2),BT1,EA)
      IF ( BT1.EQ.0.D0 )                           GO TO 9999
C--
      APV = UABV3D(PV(1,3))
      SN1 = 2*SH1*CH1
      PV(1,3) = APV*SN1*COS(FI1)
      PV(2,3) = APV*SN1*SIN(FI1)
C
C-- (b C) --> b + C.
C
      AM12 = AMF(IDP(1,4),IDP(2,4),IDP(3,4))**2
      AM22 = Q32
      CALL UHPHS2(2,QV(0,2),AM12,AM22,EC,CS2,FI2,
     .            PV(0,4),QV(0,3),BT2,EB)
      IF ( BT2.EQ.0.D0 )                           GO TO 9999
C
C-- C --> f-bar + f.
C
      AM12 = AMF(IDP(1,5),IDP(2,5),IDP(3,5))**2
      AM22 = AMF(IDP(1,6),IDP(2,6),IDP(3,6))**2
      CALL UHPHS2(2,QV(0,3),AM12,AM22,EB,X(2,3),X(3,3),
     .            PV(0,5),PV(0,6),BT3,EA)
      IF ( BT3.EQ.0.D0 )                           GO TO 9999
C--
C  Calculate phase space weight.
C--
      WAT = WAT*FACT*BT1*BT2*BT3/(2*S*BT0)
C--
C  Calculate amplitude squared.
C--
C>>>
CDBG       QV(0,5) = PV(0,3) + PV(0,4) + PV(0,5) + PV(0,6)
CDBG       QV(1,5) = PV(1,3) + PV(1,4) + PV(1,5) + PV(1,6)
CDBG       QV(2,5) = PV(2,3) + PV(2,4) + PV(2,5) + PV(2,6)
CDBG       QV(3,5) = PV(3,3) + PV(3,4) + PV(3,5) + PV(3,6)
CDBG       APVSUM  = UABV3D(QV(1,5))
CDBG       APV3    = UABV3D(PV(1,3))
CDBG       APV4    = UABV3D(PV(1,4))
CDBG       APV5    = UABV3D(PV(1,5))
CDBG       APV6    = UABV3D(PV(1,6))
CDBG       IF ( APVSUM.GT.5.E-3 .OR. 
CDBG      .     APV3.EQ.0. .OR. APV4.EQ.0. .OR. 
CDBG      .     APV5.EQ.0. .OR. APV6.EQ.0. ) THEN
CDBG       PRINT *, ' '
CDBG       PRINT *, ' --- '
CDBG       PRINT *, ' SH1, CH1, FI1 = ', SH1, CH1, FI1
CDBG       PRINT *, ' '
CDBG       PRINT *, ' IHEL = ', (IHEL(K),K=1,6)
CDBG       PRINT *, ' PV(*,1) = ', (PV(K,1),K=0,3)
CDBG       PRINT *, '   (*,2) = ', (PV(K,2),K=0,3)
CDBG       PRINT *, '   (*,3) = ', (PV(K,3),K=0,3)
CDBG       PRINT *, '   (*,4) = ', (PV(K,4),K=0,3)
CDBG       PRINT *, '   (*,5) = ', (PV(K,5),K=0,3)
CDBG       PRINT *, '   (*,6) = ', (PV(K,6),K=0,3)
CDBG       PRINT *, ' '
CDBG       PRINT *, ' PV_SUM  = ', (QV(K,5),K=0,3)
CDBG       PRINT *, ' '
CDBG       ENDIF
C>>>
      CALL FULNNZ(IDP,IHEL,PV,AMP2)
C>>>
CCDBG       IF ( AMP2(0).EQ.0. .OR. AMP2(1).EQ.0. .OR.
CCDBG      .     AMP2(2).EQ.0. .OR. AMP2(3).EQ.0. ) THEN
CCDBG          PRINT *, ' AMP2 = ',  AMP2
CCDBG       ENDIF
C>>>
C--
C  Differenctial cross section.
C--
      SG  = AMP2(0)*WAT
CE       SG = AMP2(1)*WAT
CN       SG = AMP2(2)*WAT
CW       SG = AMP2(3)*WAT
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
