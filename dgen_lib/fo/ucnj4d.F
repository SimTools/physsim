CC********************************************************************CC
C*                                                                    *C
C*=======================--====                                       *C
C*  Subroutine UCNJ4D(P1,P2)                                          *C
C*=======================--====                                       *C
C*                                                                    *C
C* (Purpose)                                                          *C
C*   Take Lorentz conjugate.                                          *C
C* (Input)                                                            *C
C*       P1(4) = an input 4-vector.                                   *C
C* (Output)                                                           *C
C*       P2(4) = its Lorentz conjugated 4-vector                      *C
C* (Relation)                                                         *C
C*   Calls no subroutines.                                            *C
C* (Update Record)                                                    *C
C*    7/29/85  K. Fujii  Original version.                            *C
C*                                                                    *C
CC********************************************************************CC
 
      SUBROUTINE UCNJ4D(P1,P2)
 
      IMPLICIT   REAL*8 ( A-H, O-Z )
      REAL*8     P1(4), P2(4)
      P2(1) = -P1(1)
      P2(2) = -P1(2)
      P2(3) = -P1(3)
      P2(4) =  P1(4)
      RETURN
      END
