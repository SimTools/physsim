CC**********************************************************************
C*
C*=====================================-----===
C* Subroutine SGNXXF(IDP,IHEL,AM,IPV,X,SG,PV)
C*=====================================-----===
C*
C* (Purpose)
C*    Calculate differential cross section for e+ + e- ---> X0i + X0j
C*    with their subsequent decays.
C* (Inputs)
C*       IDP(1,i) : (I*4) : generation number.
C*          (2,i) : (I*4) : T3 + 1.5.
C*          (3,i) : (I*4) : l/q flag; (1,2) = (l,q).
C*       IHEL( i) : (I*4) : helicity combination.
C*                        : i = (1,2,3,4,5,6)
C*                        :   = (e-,e+,tau-,tau+,X0,tau-,tau+,X0)
C*                        :            <-- X0i ---> <-- X0j --->
C*       AM(i)    : (R*8) : mass of i.
C*       IPV(1,1) : (I*4) : a from X0i
C*          (2,1) : (I*4) : b
C*          (3,1) : (I*4) : c
C*          (1,2) : (I*4) : a from X0j
C*          (2,2) : (I*4) : b
C*          (3,2) : (I*4) : c
C*       X(1,1)   : (R*8) : s.
C*        (2,1)   : (R*8) : cos(theta_X0i).
C*        (3,1)   : (R*8) : phi_X0i.
C*       X(1,2)   : (R*8) : invariant mass squared for X0i.
C*        (2,2)   : (R*8) : cos(theta_a) in X0i rest frame.
C*        (3,2)   : (R*8) : phi_a in X0i rest frame.
C*       X(1,3)   : (R*8) : invariant mass squared for ab.
C*        (2,3)   : (R*8) : invariant mass squared for ac.
C*        (3,3)   : (R*8) : phi_b in X0i rest frame.
C*       X(1,4)   : (R*8) : invariant mass squared for X0j.
C*        (2,4)   : (R*8) : cos(theta_a) in X0j rest frame.
C*        (3,4)   : (R*8) : phi_a in X0j rest frame.
C*       X(1,5)   : (R*8) : invariant mass squared for ab.
C*        (2,5)   : (R*8) : invariant mass squared for ac.
C*        (3,5)   : (R*8) : phi_b in X0j rest frame.
C* (Output)
C*       PV(*,i)  : (R*8) : 4-momentum of i-th particle, where
C*                        : numbering convention is that of IDP.
C*       SG       : (R*8) : d(sigma)/dX1dX2....dXn.
C* (Relation)
C*    Invokes UVZERO, UHSETF, UBTRAN, FULXNXN, and routines
C*    in the HELAS library.
C* (Update Record)
C*    03/08/28  K.Fujii		Original version.
C*
CC**********************************************************************
 
      SUBROUTINE SGNXXF(IDP,IHEL,AM,IPV,X,SG,PV)
 
      IMPLICIT     REAL*8  ( A-H, O-Z )
#include "SMCONS.inc"
#include "SMPTAB.inc"
#include "SSPTAB.inc"
      INTEGER*4    IDP(3,*), IHEL(*), IPV(3,*)
      REAL   *8    AM(*), X(3,*), SG, PV(0:3,*)
      REAL   *8    QV(0:3,8), EC(3,3), EB(3,3), EA(3,3), AMP2(0:3)
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
         FACT   = xGV2PB*2/(x4PI**4)/SPIN
         AME    = AMF(1,2,1)
         CALL UVZERO(2*9,EC)
         EC(1,1) = 1
         EC(2,2) = 1
         EC(3,3) = 1
      ENDIF
C--
C  Set 4-momenta.
C--
      CALL UVZERO(2*32,PV)
      CALL UVZERO(2*32,QV)
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
C
C-- CM --> X0i X0j
C
      AM12 = X(1,2)
      AM22 = X(1,4)
      CALL UHPHS2(1,QV(0,1),AM12,AM22,EC,X(2,1),X(3,1),
     .            QV(0,3),QV(0,4),BT1,EA)
      IF ( BT1.EQ.0.D0 )                           GO TO 9999
      WAT = WAT*BT1
      CALL UVCOPY(2*9,EC,EA)
      CALL UVCOPY(2*9,EC,EB)
C
C-- X0i --> a + b + c
C
      AM1 = AM(IPV(1,1))
      AM2 = AM(IPV(2,1))
      AM3 = AM(IPV(3,1))
      Q12 = SQRT(X(1,3))
      Q13 = SQRT(X(2,3))
C>>>
#ifdef __DEBUG__
      PRINT *, ' P_X0i = ', (QV(K,3),K=0,3)
      print *, ' IPV(*) = ', (IPV(K,1),K=1,3)
      PRINT *, ' AM1, AM2, AM3 = ', AM1, AM2, AM3
      PRINT *, ' CS1, FI1, FI2 = ', X(2,2), X(3,2), X(3,3)
      PRINT *, 'RS = ', RS, ' Q12 = ', Q12, ' Q13 = ', Q13
#endif
C>>>      
      CALL UHPHS3(QV(0,3),AM1,AM2,AM3,X(2,2),
     .            X(3,2),X(3,3),Q12,Q13,WT,QV(0,5))
      IF ( WT.EQ.0.D0 )                            GO TO 9999
      WAT = WAT*WT
C--
      CALL UVCOPY(2*4,QV(0,5),PV(0,IPV(1,1)))
      CALL UVCOPY(2*4,QV(0,6),PV(0,IPV(2,1)))
      CALL UVCOPY(2*4,QV(0,7),PV(0,IPV(3,1)))
C
C-- X0j --> a + b + c
C
      AM1 = AM(IPV(1,2))
      AM2 = AM(IPV(2,2))
      AM3 = AM(IPV(3,2))
      Q12 = SQRT(X(1,5))
      Q13 = SQRT(X(2,5))
C>>>
C>>>     
#ifdef __DEBUG__
      PRINT *, ' P_X0j = ', (QV(K,4),K=0,3)
      print *, ' IPV(*) = ', (IPV(K,2),K=1,3)
      PRINT *, ' AM1, AM2, AM3 = ', AM1, AM2, AM3
      PRINT *, ' CS1, FI1, FI2 = ', X(2,4), X(3,4), X(3,5)
      PRINT *, 'RS = ', RS, ' Q12 = ', Q12, ' Q13 = ', Q13
#endif
C>>>
      CALL UHPHS3(QV(0,4),AM1,AM2,AM3,X(2,4),
     .            X(3,4),X(3,5),Q12,Q13,WT,QV(0,5))
      IF ( WT.EQ.0.D0 )                            GO TO 9999
      WAT = WAT*WT
C--
      CALL UVCOPY(2*4,QV(0,5),PV(0,IPV(1,2)))
      CALL UVCOPY(2*4,QV(0,6),PV(0,IPV(2,2)))
      CALL UVCOPY(2*4,QV(0,7),PV(0,IPV(3,2)))
C--
C  Calculate phase space weight.
C--
      WAT = WAT*FACT/(2*S*BT0)
C--
C  Calculate amplitude squared.
C--
CPHASE                SG  = WAT
CPHASE                RETURN
C>>>
#ifdef __DEBUG__
      PRINT *, ' QXI = ', SQRT(AM12), ' QXJ = ', SQRT(AM22)
      PRINT *, ' PV_3 = ', (PV(K,3),K=0,3)
      PRINT *, ' PV_4 = ', (PV(K,4),K=0,3)
      PRINT *, ' PV_5 = ', (PV(K,5),K=0,3)
      PRINT *, ' PV_6 = ', (PV(K,6),K=0,3)
      PRINT *, ' PV_7 = ', (PV(K,7),K=0,3)
      PRINT *, ' PV_8 = ', (PV(K,8),K=0,3)
      DO 66 I = 3, 8
         DO 6 K = 0, 3
            QV(K,2) = QV(K,2) + PV(K,I)
6        CONTINUE
66    CONTINUE
      PRINT *, ' SUM = ', (QV(K,2),K=0,3)
#endif
C>>>
      CALL FULNXX(IDP,IHEL,PV,AMP2)
C>>>
#ifdef __DEBUG__
      PRINT *, ' WAT = ', WAT, ' AMP2 = ', AMP2
#endif
C>>>
C--
C  Differenctial cross section.
C--
      SG  = AMP2(0)*WAT
CSZ0ONLY                SG  = AMP2(1)*WAT
CSELONLY                SG  = AMP2(2)*WAT
CSERONLY                SG  = AMP2(3)*WAT
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
