CC**********************************************************************
C*
C*=====================================-----===
C* Subroutine SGNXXF(IDP,IHEL,AM,X,SG,PV)
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
C*                        :   = (e-,e+,tau-/+,stau+/-, tau+/-,stau-/+)
C*                        :            <--- X0i ---->  <--- X0j ---->
C*       AM(i)    : (R*8) : mass of i.
C*       X(1,1)   : (R*8) : s.
C*        (2,1)   : (R*8) : cos(theta_X0i).
C*        (3,1)   : (R*8) : phi_X0i.
C*       X(1,2)   : (R*8) : invariant mass squared for X0i.
C*        (2,2)   : (R*8) : cos(theta_tau-/+) in X0i rest frame.
C*        (3,2)   : (R*8) : phi_tau-/+ in X0i rest frame.
C*       X(1,3)   : (R*8) : invariant mass squared for X0j.
C*        (2,3)   : (R*8) : cos(theta_tau+/-) in X0j rest frame.
C*        (3,3)   : (R*8) : phi_tua+/- in X0j rest frame.
C* (Output)
C*       PV(*,i)  : (R*8) : 4-momentum of i-th particle, where
C*                        : numbering convention is that of IDP.
C*       SG       : (R*8) : d(sigma)/dX1dX2....dXn.
C* (Relation)
C*    Invokes UVZERO, UHSETF, UBTRAN, FULXNXN, and routines
C*    in the HELAS library.
C* (Update Record)
C*    2010/10/07  K.Fujii		Original version.
C*
CC**********************************************************************
 
      SUBROUTINE SGNXXF(IDP,IHEL,AM,X,SG,PV)
 
      IMPLICIT     REAL*8  ( A-H, O-Z )
#include "SMCONS.inc"
#include "SMPTAB.inc"
#include "SSPTAB.inc"
      INTEGER*4    IDP(3,*), IHEL(*)
      REAL   *8    AM(*), X(3,*), SG, PV(0:3,*)
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
C
C-- CM --> X0i X0j
C
      AM12 = X(1,2)
      AM22 = X(1,3)
      CALL UHPHS2(1,QV(0,1),AM12,AM22,EC,X(2,1),X(3,1),
     .            QV(0,2),QV(0,3),BT1,EA)
      IF ( BT1.EQ.0.D0 )                           GO TO 9999
      WAT = WAT*BT1
      CALL UVCOPY(2*9,EC,EA)
      CALL UVCOPY(2*9,EC,EB)
C
C-- X0i --> tau-/+ stau+/-
C
      AM12 = AM(3)**2
      AM22 = AM(4)**2
      CALL UHPHS2(2,QV(0,2),AM12,AM22,EC,X(2,2),X(3,2),
     .            PV(0,3),PV(0,4),BT2,EA)
      IF ( BT2.EQ.0. )                           GO TO 9999
      WAT = WAT*BT2
C
C-- X0j --> tau+/- stau-/+
C
      AM12 = AM(5)**2
      AM22 = AM(6)**2
      CALL UHPHS2(2,QV(0,3),AM12,AM22,EC,X(2,3),X(3,3),
     .            PV(0,5),PV(0,6),BT3,EA)
      IF ( BT3.EQ.0. )                           GO TO 9999
      WAT = WAT*BT3
C--
C  Calculate phase space weight.
C--
      WAT = WAT*FACT/(2*S*BT0)
C--
C  Calculate amplitude squared.
C--
C>>>
#ifdef __DEBUG__
      PRINT *, ' PV_3 = ', (PV(K,3),K=0,3)
      PRINT *, ' PV_4 = ', (PV(K,4),K=0,3)
      PRINT *, ' PV_5 = ', (PV(K,5),K=0,3)
      PRINT *, ' PV_6 = ', (PV(K,6),K=0,3)
#endif
#if 0
      CALL PSUMXX(PV(0,3),PV(0,4),QV(0,3))
      CALL PSUMXX(PV(0,5),PV(0,6),QV(0,4))
      CALL PSUMXX(QV(0,3),QV(0,4),QV(0,5))
      PRINT *, ' QV_ALL = ', (QV(K,5),K=0,3)
#endif
C>>>
      CALL FULNXX(IDP,IHEL,PV,AMP2)
C--
C  Differenctial cross section.
C--
#if 1
      SG  = AMP2(0)*WAT
#else
CS      SG  = AMP2(1)*WAT
CTEL      SG  = AMP2(2)*WAT
CTER      SG  = AMP2(3)*WAT
#endif
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
