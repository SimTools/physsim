CC********************************************************************CC
C*                                                                    *C
C* ==========-----=======                                             *C
C*  Function UDOT2D(X,Y)                                              *C
C* ==========-----=======                                             *C
C*                                                                    *C
C*    1) Purpose                                                      *C
C*        Calculates scalar-product of X and Y.                       *C
C*    2) I/O specification                                            *C
C*        Input                                                       *C
C*          X(3)   ; 2-dim. vector                                    *C
C*          Y(3)   ; 2-dim. vector                                    *C
C*        Output                                                      *C
C*          UDOT2D ; scalar-product of X and Y.                       *C
C*    3) Update record                                                *C
C*         6/12/85 TKSF Original version.                             *C
C*                                                                    *C
CC********************************************************************CC
 
      DOUBLE PRECISION FUNCTION UDOT2D(X,Y)
 
      IMPLICIT   REAL*8 ( A-H, O-Z )
      REAL*8     X(2), Y(2)
      UDOT2D = 0
      DO 10 I = 1, 2
         UDOT2D = UDOT2D + X(I)*Y(I)
10    CONTINUE
      RETURN
      END
