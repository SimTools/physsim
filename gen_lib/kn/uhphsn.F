C     NP = # external particles.
C     IB = 1, ..., NP-3 : branch sequence number.
C
C   Branch sequence is controlled by
C
C     IPBR(0,IB)  = JB for parent branch.
C         (1,IB)  = JB for 1-st daughter branch.
C         (2,IB)  = JB for 2-nd daughter branch.
C   where
C
C     JB = -1,...,-NP+4 : subsequent branches.
C        = 0            : the first branch from CM to 2 branches.
C        = 1, 2         : intial state particles.
C        = 3, ..., NP   : final state particles.
C
C   which point to
C
C     XMASS(JB)   = q**2 for branch JB if JB =0 or <0.
C          (JB)   = m**2 for particle JB if JB > 0.
C     XANGL(1,JB) = cos(theta) for branch JB( JB =0 or <0 ).
C          (2,JB) = phi        for branch JB( JB =0 or <0 ).
C
CC**********************************************************************
C*
C*========================================------===
C* Subroutine UHPHSN(NP,IPBR,XMASS2,XANGL,WAT,PV)
C*========================================------===
C*
C* (Purpose)
C*    Calculate (NP-2)-body phase space.
C* (Inputs)
C*       NP          : (I*4) : # external particles.
C*       IPBR  (0,i) : (I*4) : parent # for i-th branch.
C*             (1,i) : (I*4) : 1-st daughter # for i-th branch.
C*             (2,i) : (I*4) : 2-nd daughter # for i-th branch.
C*       XMASS2(  j) : (R*4) : invariant mass squared for branch #j.
C*       XANGLE(1,j) : (R*4) : cos(theta) for branch #j.
C*             (2,j) : (R*4) : phi for branch #j.
C* (Output)
C*       WAT      : (R*4) : (flux factor)*(phase space weight).
C*       PV(*,i)  : (R*4) : 4-momentum of i-th particle, where i=1,2
C*                        : should be initial state particles.
C* (Update Record)
C*    92/04/19  K.Fujii          Original version.
C*
CC**********************************************************************
 
      SUBROUTINE UHPHSN(NP,IPBR,XMASS2,XANGL,WAT,PV)
 
      IMPLICIT     REAL*4  ( A-H, O-Z    )
      INTEGER*4    NP, IPBR(0:2,NP-3)
      REAL   *4    XMASS2(-NP+4:NP), XANGL(2,-NP+4:0), WAT, PV(0:3,NP)
C--
C  MXxNP is the maximum number of external particles.
C--
      PARAMETER    ( MXxNP = 20 )
      REAL   *4    QV(0:3,-MXxNP:MXxNP), EE(3,3,-MXxNP:MXxNP)
      DATA NCALL /0/
C--
C  Statement function.
C--
      BETA(X1,X2) = SQRT(MAX(0.,1-2*(X1+X2)+(X1-X2)**2))
C
C========< Entry Point >================================================
C
C--
C  Initialize CM frame axes.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL  = 1
         CALL UVZERO(9,EE(1,1,0))
         EE(1,1,0) = 1
         EE(2,2,0) = 1
         EE(3,3,0) = 1
         xPI  = ACOS(-1.)
         x2PI = 2*xPI
         x4PI = 4*xPI
      ENDIF
C--
C  Clear 4-momenta.
C--
      CALL UVZERO(4*NP,PV)
      CALL UVZERO(4*(2*MXxNP+1),QV)
C--
C  Prepare initial state momenta.
C     p1 + p2 ---> CM
C--
      S       = XMASS2(0)
      AM12    = XMASS2(1)
      AM22    = XMASS2(2)
      BT0     = BETA(AM12/S,AM22/S)
      RS      = SQRT(S)
      PBM     = (RS/2)*BT0
      EBM1    = SQRT(AM12+PBM**2)
      EBM2    = SQRT(AM22+PBM**2)
      QV(0,1) = EBM1
      QV(3,1) = PBM
      QV(0,2) = EBM2
      QV(3,2) = -PBM
      QV(0,0) = RS
C--
C  Loop over final state branches.
C--
      NBR   = NP - 3
      WAT   = x2PI/(x4PI**(3*NBR))/(2*S*BT0)
      MODE  = 1
      DO 1000 IBR = 1, NBR
         IPR  = IPBR(0,IBR)
         IDT1 = IPBR(1,IBR)
         IDT2 = IPBR(2,IBR)
         AM12 = XMASS2(IDT1)
         AM22 = XMASS2(IDT2)
         CALL UHPHS2(MODE,QV(0,IPR),AM12,AM22,EE(1,1,IPR),
     .               XANGL(1,IPR),XANGL(2,IPR),QV(0,IDT1),QV(0,IDT2),
     .               BT,EE(1,1,IDT1))
         IF ( BT.EQ.0. )                         GO TO 9999
         CALL UVCOPY(9,EE(1,1,IDT1),EE(1,1,IDT2))
         WAT  = WAT*BT
         MODE = 2
1000  CONTINUE
C--
C  Store final state particle momenta in PV.
C--
      CALL UVCOPY(4*NP,QV(0,1),PV(0,1))
C--
C  That's it.
C--
      RETURN
C--
C  Kinematically forbidden.
C--
9999  WAT = 0
      RETURN
      END
