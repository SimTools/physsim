CC********************************************************************CC
C*                                                                    *C
C*==========================--==                                      *C
C*  Subroutine USTV4D(AM,P1,P2)                                       *C
C*==========================--==                                      *C
C*                                                                    *C
C* (Purpose)                                                          *C
C*    Sets 4-momentum P2 using 3-momentum P1 and mass AM.             *C
C* (Inputs)                                                           *C
C*      AM    = mass to assign to P1.                                 *C
C*      P1(3) = a 3-momentum of a particle in question.               *C
C* (Output)                                                           *C
C*      P2(4) = a 4-momentum with assinged mass and momentum.         *C
C* (Relation)                                                         *C
C*   Calls                                                            *C
C*      UCOPV3                                                        *C
C* (Update Record)                                                    *C
C*   6/21/89  K. Fujii  Original version.                             *C
C*                                                                    *C
CC********************************************************************CC
 
      SUBROUTINE USTV4D(AM,P1,P2)
 
      IMPLICIT   REAL*8  ( A-H, O-Z )
      REAL*8     P1(3), P2(4)
 
      CALL UVCOPY(6,P1(1),P2(1))
      P2(4) = AM*AM + UDOT3D(P1,P1)
      P2(4) = SQRT(P2(4))
      RETURN
      END
