CC**********************************************************************
C*
C*=====================================-----===
C* Subroutine GMCXDF(IDP,IHEL,AM,IPV,X,SG,PV)
C*=====================================-----===
C*
C* (Purpose)
C*    Calculate decay width for X+ ---> X0 fu fdb.
C* (Inputs)
C*       IDP(1,i) : (I*4) : generation number.
C*          (2,i) : (I*4) : T3 + 1.5.
C*          (3,i) : (I*4) : l/q flag; (1,2) = (l,q).
C*       IHEL( i) : (I*4) : helicity combination.
C*                        : i = (1,2,3,4)
C*                        :   = (X+,X0,fu,fdb)
C*                        :         <-- X+ -->
C*       AM(i)    : (R*4) : mass of i.
C*       IPV(1)   : (I*4) : a+
C*          (2)   : (I*4) : b+
C*          (3)   : (I*4) : c+
C*       X(1,1)   : (R*4) : m_X+^2.
C*        (2,1)   : (R*4) : cos(theta_a+).
C*        (3,1)   : (R*4) : phi_a+.
C*       X(1,2)   : (R*4) : invariant mass squared for a+b+.
C*        (2,2)   : (R*4) : invariant mass squared for a+c+.
C*        (3,2)   : (R*4) : phi_b+.
C* (Output)
C*       PV(*,i)  : (R*4) : 4-momentum of i-th particle, where
C*                        : numbering convention is that of IDP.
C*       SG       : (R*4) : d(sigma)/dX1dX2....dXn.
C* (Relation)
C*    Invokes UHSETF, UBTRAN, FULCXD, and library routines in FORTLIB.
C* (Update Record)
C*    95/04/09  K.Fujii		Original version.
C*
CC**********************************************************************
 
      SUBROUTINE GMCXDF(IDP,IHEL,AM,IPV,X,SG,PV)
 
      IMPLICIT     REAL*4  ( A-H, O-Z )
      INTEGER*4    IDP(3,*), IHEL(*), IPV(3)
      REAL   *4    AM(*), X(3,*), SG, PV(0:3,*)
      REAL   *4    QV(0:3,4), AMP2(0:3)
C
C========< Entry Point >================================================
C
C--
C  Set 4-momenta.
C--
      CALL UVZERO(16,PV)
      CALL UVZERO(16,QV)
C--
      WAT = 1
C
C-- X+.
C
      AMX2    = X(1,1)
      AMX     = SQRT(AMX2)
      PV(0,1) = AMX
      QV(0,1) = AMX
C
C-- X+ --> a + b + c
C
      AM1 = AM(IPV(1))
      AM2 = AM(IPV(2))
      AM3 = AM(IPV(3))
      Q12 = SQRT(X(1,2))
      Q13 = SQRT(X(2,2))
C>>>
CCDBG       PRINT *, ' P_X+ = ', (QV(K,1),K=0,3)
CCDBG       PRINT *, ' AM1, AM2, AM3 = ', AM1, AM2, AM3
CCDBG       PRINT *, ' CS1, FI1, FI2 = ', X(2,1), X(3,1), X(3,2)
CCDBG       PRINT *, ' Q12 = ', Q12, ' Q13 = ', Q13
C>>>      
      CALL UHPHS3(QV(0,1),AM1,AM2,AM3,X(2,1),
     .            X(3,1),X(3,2),Q12,Q13,WT,QV(0,2))
      IF ( WT.EQ.0. )                            GO TO 9999
      WAT = WAT*WT
C--
      CALL UVCOPY(4,QV(0,2),PV(0,IPV(1)))
      CALL UVCOPY(4,QV(0,3),PV(0,IPV(2)))
      CALL UVCOPY(4,QV(0,4),PV(0,IPV(3)))
C--
C  Calculate phase space weight.
C--
      WAT = WAT/(2*AMX)
C--
C  Calculate amplitude squared.
C--
C>>>
CCDBG        PRINT *, ' PV_2 = ', (PV(K,2),K=0,3)
CCDBG        PRINT *, ' PV_3 = ', (PV(K,3),K=0,3)
CCDBG        PRINT *, ' PV_4 = ', (PV(K,4),K=0,3)
CCDBG        PRINT *, ' IDP  = ', ((IDP(K,L),K=1,3),L=1,4)
CCDBG        PRINT *, ' IHEL = ', (IHEL(K),K=1,4)
C>>>
      CALL FULCXD(IDP,IHEL,PV,AMP2)
C>>>
CCDBG        PRINT *, ' WAT = ', WAT, ' AMP2 = ', AMP2
C>>>
C--
C  Differenctial width.
C--
      SG  = AMP2(0)*WAT
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
