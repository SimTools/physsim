CC**********************************************************************
C*
C*================================================------------==
C* Subroutine UHPHS2(MODE,Q,AM12,AM22,EB,CSTH,PHI,P1,P2,BTB,EA)
C*================================================------------==
C*
C* (Purpose)
C*    Calculates daughter 4-momenta (P1,P2) in the overall frame,
C*    when cos(theta) and phi (CSTH,PHI) of the first daughter in
C*    the helicity frame (EA) of the parent is given.
C*    EB is the helicity frame of the grandparent necessary to
C*    to specify azimuthal angle in the new frame (EA).
C*    Notice that the parent 4-momentum (Q) should be defined in
C*    the overall frame, since it will be used to define new Z
C*    axis.
C* (Inputs)
C*       MODE    : (I*4) : (1,2) = (no transf.,transf.)
C*       Q(0:3)  : (R*8) : parent 4-momentum.
C*       AM12    : (R*8) : mass**2 of 1-st daughter.
C*       AM22    : (R*8) : mass**2 of 2-nd daughter.
C*       EB(*,i) : (R*8) : old reference frame.
C*       CSTH    : (R*8) : cos(theta) of the 1-st daughter.
C*       PHI     : (R*8) : phi of the 1-st daughter.
C* (Outputs)
C*       P1(0:3) : (R*8) : 1-st daughter 4-mementum.
C*       P2(0:3) : (R*8) : 2-nd daughter 4-mementum.
C*       BTB     : (R*8) : beta_bar of the daughters in the parent
C*                       : rest frame.
C*       EA(*,i) : (R*8) : new reference frame.
C* (Relation)
C*    Calls RDOTxx, BOOSTx, UHSETF, and UBTRND.
C* (Update Record)
C*   90/06/22  K.Fujii	Original version.
C*   92/04/19  K.Fujii	Modify to set EA even when MODE = 1.
C*   95/03/25  K.Fujii	Internal calcuation with R*8.		
C*
CC**********************************************************************
 
      SUBROUTINE UHPHS2(MODE,Q,AM12,AM22,EB,CSTH,PHI,P1,P2,BTB,EA)
 
      IMPLICIT   REAL*8  ( A-H, O-Z )
      INTEGER*4  MODE
      REAL   *8  Q(0:3), AM12, AM22, EB(3,3), CSTH, PHI,
     .           P1(0:3), P2(0:3), BTB, EA(3,3)
C--
      REAL   *8  BETA2, X1, X2, AQ, AMQ, AMQ2, AP1
C--
C  Statement function.
C--
      BETA2(X1,X2) =  1 - 2*(X1+X2) + (X1-X2)**2
C
C========< Entry Point >================================================
C
C--
C  Calculate daughter 4-momenta in the parent helicity frame.
C--
      AQ   = SQRT(Q(1)**2+Q(2)**2+Q(3)**2)
      AMQ2 = (Q(0)-AQ)*(Q(0)+AQ)
      IF ( AMQ2.LE.0.D0 ) THEN
         PRINT *, ' >>>> Error in UHPHS2 >>>> '
         PRINT *, '   Q    = ', Q
         PRINT *, '   AMQ2 = ', AMQ2, ' AM12 = ', AM12, ' AM22 = ', AM22
         BTB = 0
         RETURN
      ENDIF
C--
      AMQ   = SQRT(AMQ2)
      BTB   = BETA2(AM12/AMQ2,AM22/AMQ2)
      IF ( BTB.LE.0.D0 ) THEN
C        PRINT *, ' >>>> Error in UHPHS2 >>>> '
C        PRINT *, '   Q    = ', Q
C        PRINT *, '   AMQ2 = ', AMQ2, ' AM12 = ', AM12, ' AM22 = ', AM22
         BTB = 0
         RETURN
      ENDIF
C--
      BTB   = SQRT(BTB)
      AP1   = (AMQ/2)*BTB
      SNTH  = SQRT((1-CSTH)*(1+CSTH))
      P1(0) = SQRT(AP1*AP1+AM12)
      P1(1) = AP1*SNTH*COS(PHI)
      P1(2) = AP1*SNTH*SIN(PHI)
      P1(3) = AP1*CSTH
C--
      P2(0) = AMQ - P1(0)
      P2(1) =     - P1(1)
      P2(2) =     - P1(2)
      P2(3) =     - P1(3)
C--
C  Transform the daughter 4-momenta to the frame in which the parent
C  4-momentum (Q) is defined.
C--
      IF ( MODE.EQ.1 ) THEN
         CALL UVCOPY(18,EB,EA)
      ELSE
         CALL UHSETF(Q(1),EB,EA)
         CALL UBTRND(P1(1),EA,P1(1))
         CALL UBTRND(P2(1),EA,P2(1))
         CALL BOOSTx(P1,Q,P1)
         CALL BOOSTx(P2,Q,P2)
      ENDIF
C--
      IF ( P1(0).LE.0.D0 .OR. P2(0).LE.0.D0 ) THEN
         BTB = 0
      ELSE
         AP = SQRT(P1(1)**2+P1(2)**2+P1(3)**2)
         IF ( P1(0).LT.AP ) P1(0) = SQRT(AP*AP+AM12)
         AP = SQRT(P2(1)**2+P2(2)**2+P2(3)**2)
         IF ( P2(0).LT.AP ) P2(0) = SQRT(AP*AP+AM22)
      ENDIF
C--
C  That's it.
C--
      RETURN
      END
