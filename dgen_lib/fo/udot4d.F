CC********************************************************************CC
C*                                                                    *C
C*===========-----=======                                             *C
C*  Function UDOT4D(A,B)                                              *C
C*===========-----=======                                             *C
C*                                                                    *C
C* (Purpose)                                                          *C
C*    routine to calculate invariant 4-scalar product.                *C
C* (Inputs)                                                           *C
C*      A(4), B(4) :  2 4-vectors.                                    *C
C* (Output)                                                           *C
C*      UDOT4D     :  A(4)*B(4)-UDOT3D(A(1),B(1))                     *C
C* (Update Record)                                                    *C
C*   7/29/85  K. Fujii  Original version.                             *C
C*                                                                    *C
CC********************************************************************CC
 
      DOUBLE PRECISION FUNCTION UDOT4D(A,B)
 
      IMPLICIT   REAL*8 ( A-H, O-Z )
      REAL*8     A(4), B(4)
      UDOT4D = A(4)*B(4)
      DO 10 I = 1, 3
         UDOT4D = UDOT4D - A(I)*B(I)
10    CONTINUE
      RETURN
      END
