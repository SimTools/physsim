CC**********************************************************************
C*
C*=====================================-----===
C* Subroutine SGCXXF(IDP,IHEL,AM,IPV,X,SG,PV)
C*=====================================-----===
C*
C* (Purpose)
C*    Calculate differential cross section for e+ + e- ---> X+ + X-
C*    with their subsequent decays.
C* (Inputs)
C*       IDP(1,i) : (I*4) : generation number.
C*          (2,i) : (I*4) : T3 + 1.5.
C*          (3,i) : (I*4) : l/q flag; (1,2) = (l,q).
C*       IHEL( i) : (I*4) : helicity combination.
C*                        : i = (1,2,3,4,5,6,7,8)
C*                        :   = (e-,e+,X0,fu,fdb,X0,fub,fd)
C*                        :            <-- X+ --><-- X- -->
C*       AM(i)    : (R*4) : mass of i.
C*       IPV(1,1) : (I*4) : a+
C*          (2,1) : (I*4) : b+
C*          (3,1) : (I*4) : c+
C*          (1,2) : (I*4) : a-
C*          (2,2) : (I*4) : b-
C*          (3,2) : (I*4) : c-
C*       X(1,1)   : (R*4) : s.
C*        (2,1)   : (R*4) : cos(theta_X-).
C*        (3,1)   : (R*4) : phi_X-.
C*       X(1,2)   : (R*4) : invariant mass squared for X+.
C*        (2,2)   : (R*4) : cos(theta_a+) in X+ rest frame.
C*        (3,2)   : (R*4) : phi_a+ in X+ rest frame.
C*       X(1,3)   : (R*4) : invariant mass squared for a+b+.
C*        (2,3)   : (R*4) : invariant mass squared for a+c+.
C*        (3,3)   : (R*4) : phi_b+ in X+ rest frame.
C*       X(1,4)   : (R*4) : invariant mass squared for X-.
C*        (2,4)   : (R*4) : cos(theta_a-) in X- rest frame.
C*        (3,4)   : (R*4) : phi_a- in X- rest frame.
C*       X(1,5)   : (R*4) : invariant mass squared for a-b-.
C*        (2,5)   : (R*4) : invariant mass squared for a-c-.
C*        (3,5)   : (R*4) : phi_b- in X- rest frame..
C* (Output)
C*       PV(*,i)  : (R*4) : 4-momentum of i-th particle, where
C*                        : numbering convention is that of IDP.
C*       SG       : (R*4) : d(sigma)/dX1dX2....dXn.
C* (Relation)
C*    Invokes UVZERO, UHSETF, UBTRAN, FULXCXC, and library routines
C*    in 'T#FP.HELAS.FORT'.
C* (Update Record)
C*    92/06/09  K.Fujii		Original version.
C*    95/04/06  K.Fujii		New kinematics.
C*    00/06/24  K.Fujii		Use fb instead of pb.
C*
CC**********************************************************************
 
      SUBROUTINE SGCXXF(IDP,IHEL,AM,IPV,X,SG,PV)
 
      IMPLICIT     REAL*4  ( A-H, O-Z )
#include "SMCONS.inc"
#include "SMPTAB.inc"
#include "SSPTAB.inc"
      INTEGER*4    IDP(3,*), IHEL(*), IPV(3,*)
      REAL   *4    AM(*), X(3,*), SG, PV(0:3,*)
      REAL   *4    QV(0:3,8), EC(3,3), EB(3,3), EA(3,3), AMP2(0:3)
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
         FACT   = xGV2FB*2/(x4PI**4)/SPIN
         AME    = AMF(1,2,1)
         CALL UVZERO(9,EC)
         EC(1,1) = 1
         EC(2,2) = 1
         EC(3,3) = 1
      ENDIF
C--
C  Set 4-momenta.
C--
      CALL UVZERO(32,PV)
      CALL UVZERO(32,QV)
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
C-- CM --> X+ X-
C
      AM12 = X(1,4)
      AM22 = X(1,2)
      CALL UHPHS2(1,QV(0,1),AM12,AM22,EC,X(2,1),X(3,1),
     .            QV(0,4),QV(0,3),BT1,EA)
      IF ( BT1.EQ.0. )                           GO TO 9999
      WAT = WAT*BT1
      CALL UVCOPY(9,EC,EA)
      CALL UVCOPY(9,EC,EB)
C
C-- X+ --> a + b + c
C
      AM1 = AM(IPV(1,1))
      AM2 = AM(IPV(2,1))
      AM3 = AM(IPV(3,1))
      Q12 = SQRT(X(1,3))
      Q13 = SQRT(X(2,3))
C>>>
CCDBG       PRINT *, ' P_X+ = ', (QV(K,3),K=0,3)
CCDBG       PRINT *, ' AM1, AM2, AM3 = ', AM1, AM2, AM3
CCDBG       PRINT *, ' CS1, FI1, FI2 = ', X(2,2), X(3,2), X(3,3)
CCDBG       PRINT *, ' Q12 = ', Q12, ' Q13 = ', Q13
C>>>      
      CALL UHPHS3(QV(0,3),AM1,AM2,AM3,X(2,2),
     .            X(3,2),X(3,3),Q12,Q13,WT,QV(0,5))
      IF ( WT.EQ.0. )                            GO TO 9999
      WAT = WAT*WT
C--
      CALL UVCOPY(4,QV(0,5),PV(0,IPV(1,1)))
      CALL UVCOPY(4,QV(0,6),PV(0,IPV(2,1)))
      CALL UVCOPY(4,QV(0,7),PV(0,IPV(3,1)))
C
C-- X- --> a + b + c
C
      AM1 = AM(IPV(1,2))
      AM2 = AM(IPV(2,2))
      AM3 = AM(IPV(3,2))
      Q12 = SQRT(X(1,5))
      Q13 = SQRT(X(2,5))
C>>>
CCDBG       PRINT *, ' P_X- = ', (QV(K,4),K=0,3)
CCDBG       PRINT *, ' AM1, AM2, AM3 = ', AM1, AM2, AM3
CCDBG       PRINT *, ' CS1, FI1, FI2 = ', X(2,4), X(3,4), X(3,5)
CCDBG       PRINT *, ' Q12 = ', Q12, ' Q13 = ', Q13
C>>>      
      CALL UHPHS3(QV(0,4),AM1,AM2,AM3,X(2,4),
     .            X(3,4),X(3,5),Q12,Q13,WT,QV(0,5))
      IF ( WT.EQ.0. )                            GO TO 9999
      WAT = WAT*WT
C--
      CALL UVCOPY(4,QV(0,5),PV(0,IPV(1,2)))
      CALL UVCOPY(4,QV(0,6),PV(0,IPV(2,2)))
      CALL UVCOPY(4,QV(0,7),PV(0,IPV(3,2)))
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
CCDBG        PRINT *, ' PV_3 = ', (PV(K,3),K=0,3)
CCDBG        PRINT *, ' PV_4 = ', (PV(K,4),K=0,3)
CCDBG        PRINT *, ' PV_5 = ', (PV(K,5),K=0,3)
CCDBG        PRINT *, ' PV_6 = ', (PV(K,6),K=0,3)
CCDBG        PRINT *, ' PV_7 = ', (PV(K,7),K=0,3)
CCDBG        PRINT *, ' PV_8 = ', (PV(K,8),K=0,3)
C>>>
      CALL FULCXX(IDP,IHEL,PV,AMP2)
C>>>
CCDBG        PRINT *, ' WAT = ', WAT, ' AMP2 = ', AMP2
C>>>
C--
C  Differenctial cross section.
C--
      SG  = AMP2(0)*WAT
CAONLY                SG  = AMP2(1)*WAT
CZONLY                SG  = AMP2(2)*WAT
CTONLY                SG  = AMP2(3)*WAT
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
