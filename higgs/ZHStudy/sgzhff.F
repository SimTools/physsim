CC**********************************************************************
C*
C*==============================-----===
C* Subroutine SGZHFF(IDP,IHEL,X,SG,PV)
C*==============================-----===
C*
C* (Purpose)
C*    Calculate amplitude for e+ + e- ---> Z + H.
C* (Inputs)
C*       IDP(1,i) : (I*4) : generation number.
C*          (2,i) : (I*4) : T3 + 1.5.
C*          (3,i) : (I*4) : l/q flag; (1,2) = (l,q).
C*       IHEL( i) : (I*4) : helicity combination.
C*                        : i = (1, 2,  3, 4, 5, 6)
C*                        :   = (e-,e+, bb,b, mb,m)
C*                        :            <-H->  <-Z->
C*       X(1,1)   : (R*4) : s.
C*        (2,1)   : (R*4) : cos(theta_H).
C*        (3,1)   : (R*4) : phi_H.
C*       X(1,2)   : (R*4) : invariant mass squared for b bb system.
C*        (2,2)   : (R*4) : cos(theta_bb) in b bb rest frame.
C*        (3,2)   : (R*4) : phi_bb in b bb rest frame.
C*       X(1,3)   : (R*4) : invariant mass squared for m mb system.
C*        (2,3)   : (R*4) : cos(theta_mb) in Z rest frame.
C*        (3,3)   : (R*4) : phi_mb in Z rest frame.
C* (Output)
C*       PV(*,i)  : (R*4) : 4-momentum of i-th particle, where
C*                        : nubering convention is that of IDP.
C*       SG       : (R*4) : d(sigma)/dX1dX2....dXn.
C* (Relation)
C*    Invokes UVZERO, UHSETF, UBTRAN, FULZHF, and library routines
C*    in 'T#FP.HELAS.FORT'.
C* (Update Record)
C*    92/02/24  K.Fujii          Original version.
C*    95/03/25  K.Fujii          Use INSMCN instead of CONINT.
C*    00/06/24  K.Fujii          Use xGV2FB instead of xGV2PB.
C*
CC**********************************************************************
 
      SUBROUTINE SGZHFF(IDP,IHEL,X,SG,PV)
 
      IMPLICIT     REAL*4  ( A-H, O-Z )
#include "SMCONS.inc"
#include "SMPTAB.inc"
      PARAMETER    ( IxPRC = 1 )
      INTEGER*4    IDP(3,*), IHEL(*)
      REAL   *4    X(3,*), SG, PV(0:3,*)
      REAL   *4    QV(0:3,6), EC(3,3), EA(3,3), AMP2
      DATA NCALL /0/
C
C========< Entry Point >================================================
C
C--
C  Constants.
C  Spin average is taken only for positron, since electron is
C  assumed to be polarized and its spin average is taken in FUNCZHFF.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL  = 1
         SPIN   = 2
         NB     = 3
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
      CALL UVZERO(24*IxPRC,PV)
      CALL UVZERO(24*IxPRC,QV)
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
C-- CM --> H + Z.
C
      AM12 = X(1,2)
      AM22 = X(1,3)
      CALL UHPHS2(1,QV(0,1),AM12,AM22,EC,X(2,1),X(3,1),
     .            QV(0,2),QV(0,3),BT1,EA)
      IF ( BT1.EQ.0. )                           GO TO 9999
C
C-- H --> b + bb.
C
      AM12 = AMF(IDP(1,3),IDP(2,3),IDP(3,3))**2
      AM22 = AMF(IDP(1,4),IDP(2,4),IDP(3,4))**2
      CALL UHPHS2(2,QV(0,2),AM12,AM22,EC,X(2,2),X(3,2),
     .            PV(0,3),PV(0,4),BT2,EA)
      IF ( BT2.EQ.0. )                           GO TO 9999
C
C-- Z --> m + mb.
C
      AM12 = AMF(IDP(1,5),IDP(2,5),IDP(3,5))**2
      AM22 = AMF(IDP(1,6),IDP(2,6),IDP(3,6))**2
      CALL UHPHS2(2,QV(0,3),AM12,AM22,EC,X(2,3),X(3,3),
     .            PV(0,5),PV(0,6),BT3,EA)
      IF ( BT3.EQ.0. )                           GO TO 9999
C--
C  Calculate phase space weight.
C--
      WAT = FACT*BT1*BT2*BT3/(2*S*BT0)
C--
C  Calculate amplitude squared.
C--
      CALL FULZHF(IDP,IHEL,PV,AMP2)
C--
C  Differenctial cross section.
C--
      SG  = AMP2*WAT
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
