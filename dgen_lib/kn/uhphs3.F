CC**********************************************************************
C*
C*=======================================================------===
C* Subroutine UHPHS3(QPR,AM1,AM2,AM3,CS1,FI1,FI2,Q12,Q13,WAT,PV)
C*=======================================================------===
C*
C* (Purpose)
C*    Set 4-momenta for QPR ---> AM1 + AM2 + AM3.
C* (Inputs)
C*       QPR(*)  : (R*8) : parent 4-momentum.
C*       AM1     : (R*8) : 1st daughter mass.
C*       AM2     : (R*8) : 2nd daughter mass.
C*       AM3     : (R*8) : 3rd daughter mass.
C*       CS1     : (R*8) : cos(theta_1).
C*       FI1     : (R*8) : phi_1.
C*       FI2     : (R*8) : phi_2 measured from PR-3 plane.
C*       Q12     : (R*8) : invariant mass for 1-2 pair.
C*       Q13     : (R*8) : invariant mass for 1-3 pair.
C* (Output)
C*       WAT     : (R*8) : phase space weight.
C*       PV(*,1) : (R*8) : 1st daughter 4-momentum.
C*         (*,2) : (R*8) : 2nd daughter 4-momentum.
C*         (*,3) : (R*8) : 3rd daughter 4-momentum.
C* (Relation)
C*    Invokes UVZERO, UHSETF, UBTRND, and routines in helas_v1.
C* (Update Record)
C*    95/04/05  K.Fujii		Original version.
C*
CC**********************************************************************
 
      SUBROUTINE UHPHS3(QPR,AM1,AM2,AM3,CS1,FI1,FI2,Q12,Q13,WAT,PV)
       
      IMPLICIT REAL*8  ( A-H, O-Z )
#include "../sm/include/SMCONS.inc"
      REAL*8    QPR(0:3), AM1, AM2, AM3, CS1, FI1, FI2, Q12, Q13,
     .          WAT, PV(0:3,3)
      REAL*8    QV(0:3,3), EC(3,3), EB(3,3), EA(3,3)
      DATA NCALL /0/
C
C========< Entry Point >================================================
C
C--
C  Constants.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL   = 1
         CALL UVZERO(18,EC)
         EC(1,1) = 1
         EC(2,2) = 1
         EC(3,3) = 1
         FACT    = 1/(32*x2PI**5)
      ENDIF
C--
C  Zap 4-vectors.
C--
      CALL UVZERO(24,PV)
      CALL UVZERO(24,QV)
C--
C  Check if the parent mass is large enough.
C--
      QQ123 = AM1*AM1 + AM2*AM2 + AM3*AM3
      Q122  = Q12*Q12
      Q132  = Q13*Q13
C--      
      AMPR2 = QPR(0)*QPR(0)
      AMPR2 = AMPR2 - QPR(1)*QPR(1) - QPR(2)*QPR(2) - QPR(3)*QPR(3)
      IF ( AMPR2+QQ123.LE.Q122+Q132 ) 			GO TO 9999
C--
C  Set 4-momenta in the parent rest frame.
C     PV(*,1) = 4-momentum of 1st daughter
C       (*,2) = 4-momentum of 2nd daughter
C       (*,3) = 4-momentum of 3rd daughter
C--
      AMPR = SQRT(AMPR2)
      E2 = ( AMPR2 + AM2*AM2 - Q13*Q13 )/2/AMPR
      E3 = ( AMPR2 + AM3*AM3 - Q12*Q12 )/2/AMPR
      E1 = AMPR - E2 - E3
      P2 = ( E2 - AM2 )*( E2 + AM2 )
      IF ( P2.LT.0.D0 ) 				GO TO 9999
      P2 = SQRT(P2)
      P3 = ( E3 - AM3 )*( E3 + AM3 )
      IF ( P3.LT.0.D0 ) 				GO TO 9999
      P3 = SQRT(P3)
      P1 = ( E1 - AM1 )*( E1 + AM1 )
      IF ( P1.LT.0.D0 ) 				GO TO 9999
      P1 = SQRT(P1)
C--
      IF ( P1.GT.0.D0 ) THEN
         IF ( P2.GT.0.D0 ) THEN
            CSTH2 = -(P1*P1+P2*P2-P3*P3)/(2*P1*P2)
            IF ( ABS(CSTH2).GT.1.D0 ) CSTH2 = SIGN(1.D0,CSTH2)
         ELSE
            CSTH2 = 0  
         ENDIF
         IF ( P3.GT.0.D0 ) THEN
            CSTH3 = -(P1*P1+P3*P3-P2*P2)/(2*P1*P3)
            IF ( ABS(CSTH3).GT.1.D0 ) CSTH3 = SIGN(1.D0,CSTH3)
         ELSE
            CSTH3 = 0  
         ENDIF
      ELSE
         CSTH2 =  1
         CSTH3 = -1
      ENDIF
      SNTH2 = SQRT((1-CSTH2)*(1+CSTH2))
      SNTH3 = SQRT((1-CSTH3)*(1+CSTH3))
C--      
      PV(0,2) = E2
      PV(1,2) = P2*SNTH2*COS(FI2)
      PV(2,2) = P2*SNTH2*SIN(FI2)
      PV(3,2) = P2*CSTH2
C--      
      PV(0,3) = E3
      PV(1,3) = -PV(1,2)
      PV(2,3) = -PV(2,2)
      PV(3,3) = -(PV(3,2) + P1)
C--
C  Set reference frame.
C--
      IF ( ABS(CS1).LT.1.D0 ) THEN
         SN1 = SQRT((1-CS1)*(1+CS1))
      ELSE
         SN1 = 0
      ENDIF
      PV(1,1) = SN1*COS(FI1)
      PV(2,1) = SN1*SIN(FI1)
      PV(3,1) = CS1
C--
      CALL UHSETF(QPR(1),EC,EB)
      CALL UHSETF(PV(1,1),EB,EA)
C--
C  Transform to the parent-moving frame. 
C--      
      CALL UBTRND(PV(1,2),EA,PV(1,2))
      CALL UBTRND(PV(1,3),EA,PV(1,3))
C--
      CALL BOOSTx(PV(0,2),QPR(0),PV(0,2))
      CALL BOOSTx(PV(0,3),QPR(0),PV(0,3))
      CALL PSUBxx(QPR(0),PV(0,2),PV(0,1))
      CALL PSUBxx(PV(0,1),PV(0,3),PV(0,1))
C--
C  Calculate phase space weight.
C--
      WAT = FACT/AMPR2
      RETURN
C--
9999  WAT = 0      
C--
C  That's it.
C--
      RETURN
      END
