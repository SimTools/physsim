CC********************************************************************CC
C*                                                                    *C
C* =======================---=                                        *C
C*   Subroutine UNCRSD(A,B,C)                                         *C
C* =======================---=                                        *C
C*                                                                    *C
C* (Purpose)                                                          *C
C*   calculates 3-vector normed-cross product.                        *C
C*           C=AXB, Then C = C/!C!                                    *C
C* (Inputs)                                                           *C
C*   A,B  : 2 3-vectors to make cross product.                        *C
C* (Output)                                                           *C
C*     C  : results of the cross product with normalization.          *C
C* (Relation)                                                         *C
C*   calls                                                            *C
C*      UABV3D, UCROSD, USCM3D.                                       *C
C* (Update Record)                                                    *C
C*   7/29/85  K. Fujii  Original version.                             *C
C*                                                                    *C
CC********************************************************************CC
 
      SUBROUTINE UNCRSD(A,B,C)
 
      IMPLICIT   REAL*8 ( A-H, O-Z )
      REAL*8     A(3), B(3), C(3)
C--
C  Calculate vector product.
C--
      AMAG = UABV3D(A)
      BMAG = UABV3D(B)
      CALL UCROSD(A,B,C)
      CMAG = UABV3D(C)
      IF ( CMAG.LT.(1.0D-12*AMAG*BMAG) )          GO TO 5
      CALL USCM3D(1/CMAG,C,C)
                                                 GO TO 10
C--
C  Error return.
C--
5     CONTINUE
      WRITE(6,7)A,B
7     FORMAT(' ***UNCRSD Error...input vectors parallel...input',
     .       ' dumped ***',/'  ',3D14.7,' ',3D14.7)
10    RETURN
      END
