CC**********************************************************************
C*
C*=============================-----===
C* Subroutine SGTBWF(IDP,IHEL,X,SG,PV)
C*=============================-----===
C*
C* (Purpose)
C*    Calculate amplitude for e+ + e- ---> t + t-bar.
C* (Inputs)
C*       IDP(1,i) : (I*4) : generation number.
C*          (2,i) : (I*4) : T3 + 1.5.
C*          (3,i) : (I*4) : l/q flag; (1,2) = (l,q).
C*       IHEL( i) : (I*4) : helicity combination.
C*                        : i = (1,2,3,4,5)
C*                        :   = (e-,e+,bb,t,W-)
C*       X(1,1)   : (R*8) : s.
C*        (2,1)   : (R*8) : cos(theta_t).
C*        (3,1)   : (R*8) : phi_t.
C*       X(1,2)   : (R*8) : invariant mass squared for t-bar.
C*        (2,2)   : (R*8) : cos(theta_q_bar) in t-bar rest frame.
C*        (3,2)   : (R*8) : phi_q_bar in t-bar rest frame.
C* (Output)
C*       PV(*,i)  : (R*8) : 4-momentum of i-th particle, where
C*                        : nubering convention is that of IDP.
C*       SG       : (R*8) : d(sigma)/dX1dX2....dXn.
C* (Relation)
C*    Invokes UVZERO, UHSETF, UBTRAN, FULTBW, and library routines
C*    in 'T#FP.HELAS.FORT'.
C* (Update Record)
C*   2010/11/29  K.Fujii          Original version.
C*
CC**********************************************************************
 
      SUBROUTINE SGTBWF(IDP,IHEL,X,SG,PV)
 
      IMPLICIT     REAL*8  ( A-H, O-Z )
#include "SMCONS.inc"
#include "SMPTAB.inc"
      PARAMETER    ( IxPRC = 2 )
      INTEGER*4    IDP(3,*), IHEL(*)
      REAL   *8    X(3,*), SG, PV(0:3,*)
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
         NB     = 2
         FACT   = xGV2FB*x2PI/(x4PI**(3*NB))/SPIN
         AME    = AMF(1,2,1)
         CALL UVZERO(IxPRC*9,EC)
         EC(1,1) = 1
         EC(2,2) = 1
         EC(3,3) = 1
      ENDIF
C--
C  Set 4-momenta.
C--
      CALL UVZERO(20*IxPRC,PV)
      CALL UVZERO(32*IxPRC,QV)
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
C-- CM --> t t-bar
C
      AM12 = AMF(IDP(1, 4),IDP(2, 4),IDP(3, 4))**2
      AM22 = X(1,2)
      CALL UHPHS2(1,QV(0,1),AM12,AM22,EC,X(2,1),X(3,1),
     .            PV(0,4),QV(0,3),BT1,EA)
      IF ( BT1.EQ.0.D0 )                           GO TO 9999
      CALL UVCOPY(9*IxPRC,EC,EA)
      CALL UVCOPY(9*IxPRC,EC,EB)
C
C-- t-bar --> q(d)-bar + W-
C
      AM12 = AMF(IDP(1, 3),IDP(2, 3),IDP(3, 3))**2
      AM22 = AMW*AMW
      CALL UHPHS2(2,QV(0,3),AM12,AM22,EB,X(2,2),X(3,2),
     .            PV(0,3),PV(0,5),BT2,EB)
      IF ( BT2.EQ.0.D0 )                           GO TO 9999
C--
C  Calculate phase space weight.
C--
      WAT = FACT*BT1*BT2/(2*S*BT0)
C--
C  Calculate amplitude squared.
C--
      CALL FULTBW(IDP,IHEL,PV,AMP2)
C--
C  Differenctial cross section.
C--
C     SG  = AMP2(0)*WAT
C     SG  = AMP2(1)*WAT
C     SG  = AMP2(2)*WAT
      SG  = AMP2(3)*WAT
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
