CC********************************************************************CC
C*                                                                    *C
C*=========================---=====                                   *C
C*  Subroutine UADD4D(P1,P2,P12)                                      *C
C*=========================---=====                                   *C
C*                                                                    *C
C* (Purpose)                                                          *C
C*    Carries out 4-vector addition.                                  *C
C*          P1 + P2 = P12                                             *C
C* (Inputs)                                                           *C
C*      P1    = a 4-vector.                                           *C
C*      P2    = a 4-vector.                                           *C
C* (Output)                                                           *C
C*      P12   = resultant 4-vector.                                   *C
C* (Update Record)                                                    *C
C*   7/29/85  K. Fujii  Original version.                             *C
C*                                                                    *C
CC********************************************************************CC
 
      SUBROUTINE UADD4D(P1,P2,P12)
 
      IMPLICIT   REAL*8 ( A-H, O-Z )
      REAL*8     P1(4), P2(4), P12(4)
      DO 10 I = 1, 4
         P12(I) = P1(I) + P2(I)
10    CONTINUE
      RETURN
      END
