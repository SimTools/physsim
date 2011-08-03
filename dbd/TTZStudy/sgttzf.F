CC**********************************************************************
C*
C*==============================-----===
C* Subroutine SGTTZF(IDP,IHEL,X,SG,PV)
C*==============================-----===
C*
C* (Purpose)
C*    Calculate amplitude for e+ + e- ---> t + t-bar + Z.
C* (Inputs)
C*       IDP(1,i) : (I*4) : generation number.
C*          (2,i) : (I*4) : T3 + 1.5.
C*          (3,i) : (I*4) : l/q flag; (1,2) = (l,q).
C*       IHEL( i) : (I*4) : helicity combination.
C*                        : i = (1,2,3,4,5,6,7,8,9,10)
C*                        :   = (e-,e+,qdb,fub,fd,qd,fu,fdb,fb,f)
C*                        :            <--tbar--> <---t---> <-Z->
C*       X(1,1)   : (R*4) : s.
C*        (2,1)   : (R*4) : cos(theta_Z).
C*        (3,1)   : (R*4) : phi_Z.
C*       X(1,2)   : (R*4) : invariant mass squared for t t-bar system.
C*        (2,2)   : (R*4) : cos(theta_t) in t t-bar rest frame.
C*        (3,2)   : (R*4) : phi_t in t t-bar rest frame.
C*       X(1,3)   : (R*4) : invariant mass squared for t-bar.
C*        (2,3)   : (R*4) : cos(theta_q_bar) in t-bar rest frame.
C*        (3,3)   : (R*4) : phi_q_bar in t-bar rest frame.
C*       X(1,4)   : (R*4) : invariant mass squared for W-.
C*        (2,4)   : (R*4) : cos(theta_f) in W- rest frame.
C*        (3,4)   : (R*4) : phi_f in W- rest frame.
C*       X(1,5)   : (R*4) : invariant mass squared for t.
C*        (2,5)   : (R*4) : cos(theta_q) in t rest frame.
C*        (3,5)   : (R*4) : phi_q in t rest frame.
C*       X(1,6)   : (R*4) : invariant mass squared for W+.
C*        (2,6)   : (R*4) : cos(theta_f_bar) in W+ rest frame.
C*        (3,6)   : (R*4) : phi_f_bar in W+ rest frame.
C*       X(1,7)   : (R*4) : invariant mass squared for Z.
C*        (2,7)   : (R*4) : cos(theta_f_bar) in Z rest frame.
C*        (3,7)   : (R*4) : phi_f_bar in Z rest frame.
C* (Output)
C*       PV(*,i)  : (R*4) : 4-momentum of i-th particle, where
C*                        : nubering convention is that of IDP.
C*       SG       : (R*4) : d(sigma)/dX1dX2....dXn.
C* (Relation)
C*    Invokes UVZERO, UHSETF, UBTRAN, TTATTZ, and library routines
C*    in 'T#FP.HELAS.FORT'.
C* (Update Record)
C*    90/06/23  K.Fujii          Original version.
C*
CC**********************************************************************
 
      SUBROUTINE SGTTZF(IDP,IHEL,X,SG,PV)
 
      IMPLICIT     REAL*4  ( A-H, O-Z )
#include "SMCONS.inc"
#include "SMPTAB.inc"
      PARAMETER    ( IxPRC = 1 )
      INTEGER*4    IDP(3,*), IHEL(*)
      REAL   *4    X(3,*), SG, PV(0:3,*)
      REAL   *4    QV(0:3,8), EC(3,3), EB(3,3), EA(3,3), AMP2(0:2)
      DATA NCALL /0/
C
C========< Entry Point >================================================
C
C--
C  Constants.
C  Spin average is taken only for positron, since electron is
C  assumed to be polarized and its spin average is taken in TTATTZ.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL  = 1
         SPIN   = 1
         NB     = 7
         FACT   = xGV2FB*x2PI/(x4PI**(3*NB))/SPIN
         AME    = AMF(1,2,1)
         CALL UVZERO(9,EC)
         EC(1,1) = 1
         EC(2,2) = 1
         EC(3,3) = 1
      ENDIF
C--
C  Set 4-momenta.
C--
      CALL UVZERO(40*IxPRC,PV)
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
C-- CM --> Z+(tt-bar)
C
      AM12 = X(1,7)
      AM22 = X(1,2)
      CALL UHPHS2(1,QV(0,1),AM12,AM22,EC,X(2,1),X(3,1),
     .            QV(0,5),QV(0,6),BT1,EA)
      IF ( BT1.EQ.0. )                           GO TO 9999
C
C-- Z --> f+f-bar
C
      AM12 = AMF(IDP(1, 9),IDP(2, 9),IDP(3, 9))**2
      AM22 = AMF(IDP(1,10),IDP(2,10),IDP(3,10))**2
      CALL UHPHS2(2,QV(0,5),AM12,AM22,EC,X(2,7),X(3,7),
     .            PV(0,9),PV(0,10),BT7,EB)
      IF ( BT7.EQ.0. )                           GO TO 9999
C
C-- (tt-bar) --> t+t-bar
C
      AM12 = X(1,5)
      AM22 = X(1,3)
      CALL UHPHS2(2,QV(0,6),AM12,AM22,EC,X(2,2),X(3,2),
     .            QV(0,4),QV(0,3),BT2,EB)
      IF ( BT2.EQ.0. )                           GO TO 9999
      CALL UVCOPY(9,EB,EA)
C
C-- t-bar --> q(d)-bar + W-
C
      AM12 = AMF(IDP(1, 3),IDP(2, 3),IDP(3, 3))**2
      AM22 = X(1,4)
      CALL UHPHS2(2,QV(0,3),AM12,AM22,EB,X(2,3),X(3,3),
     .            PV(0,3),QV(0,7),BT3,EB)
      IF ( BT3.EQ.0. )                           GO TO 9999
C
C-- W- --> f(u)-bar + f(d)
C
      AM12 = AMF(IDP(1, 5),IDP(2, 5),IDP(3, 5))**2
      AM22 = AMF(IDP(1, 4),IDP(2, 4),IDP(3, 4))**2
      CALL UHPHS2(2,QV(0,7),AM12,AM22,EB,X(2,4),X(3,4),
     .            PV(0,5),PV(0,4),BT4,EB)
      IF ( BT4.EQ.0. )                           GO TO 9999
C
C-- t --> q(d) + W+
C
      AM12 = AMF(IDP(1, 6),IDP(2, 6),IDP(3, 6))**2
      AM22 = X(1,6)
      CALL UHPHS2(2,QV(0,4),AM12,AM22,EA,X(2,5),X(3,5),
     .            PV(0,6),QV(0,8),BT5,EA)
      IF ( BT5.EQ.0. )                           GO TO 9999
C
C-- W+ --> f(u) + f(d)-bar
C
      AM12 = AMF(IDP(1, 8),IDP(2, 8),IDP(3, 8))**2
      AM22 = AMF(IDP(1, 7),IDP(2, 7),IDP(3, 7))**2
      CALL UHPHS2(2,QV(0,8),AM12,AM22,EA,X(2,6),X(3,6),
     .            PV(0,8),PV(0,7),BT6,EA)
      IF ( BT6.EQ.0. )                           GO TO 9999
C--
C  Calculate phase space weight.
C--
      WAT = FACT*BT1*BT2*BT3*BT4*BT5*BT6*BT7/(2*S*BT0)
C--
C  Calculate amplitude squared.
C--
      CALL FULTTZ(IDP,IHEL,PV,AMP2)
C--
C  Differenctial cross section.
C--
      SG  = AMP2(0)*WAT
#ifdef NOHIGGS
      SG = AMP2(1)*WAT
#endif
#ifdef HIGGSONLY
      SG = AMP2(2)*WAT
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
