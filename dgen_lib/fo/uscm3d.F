CC********************************************************************CC
C*                                                                    *C
C* =======================--==                                        *C
C*  Subroutine USCM3D(A,X,AX)                                         *C
C* =======================--==                                        *C
C*                                                                    *C
C*    1) Purpose                                                      *C
C*        Multiplies X by a scalar A.                                 *C
C*    2) I/O specification                                            *C
C*        Input                                                       *C
C*          A      ; a scalar number                                  *C
C*          X(3)   ; 3-dim. vector                                    *C
C*        Output                                                      *C
C*          AX(3)  ; A*X(*)                                           *C
C*    4) Update record                                                *C
C*         6/12/85 TKSF Original version.                             *C
C*                                                                    *C
CC********************************************************************CC
 
      SUBROUTINE USCM3D(A,X,AX)
 
      IMPLICIT   REAL*8 ( A-H, O-Z )
      REAL*8     X(3), AX(3)
      DO 10 I = 1, 3
         AX(I) = A*X(I)
10    CONTINUE
      RETURN
      END
