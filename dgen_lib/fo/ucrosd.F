C
C ======================================
C   Double Precesion Utility Routines.
C ======================================
C
CC********************************************************************CC
C*                                                                    *C
C*=======================---===                                       *C
C*  Subroutine UCROSD(X,Y,Z)                                          *C
C*=======================---===                                       *C
C*                                                                    *C
C*    1) Purpose                                                      *C
C*        Calculates vector-product of X x Y = Z.                     *C
C*    2) I/O specification                                            *C
C*        Input                                                       *C
C*          X(3)   ; 3-dim. vector                                    *C
C*          Y(3)   ; 3-dim. vector                                    *C
C*        Output                                                      *C
C*          Z(3)   ; X x Y                                            *C
C*    3) Relation                                                     *C
C*        Calls                                                       *C
C*          UVCOPY                                                    *C
C*    4) Update record                                                *C
C*         6/12/85 TKSF Original version.                             *C
C*                                                                    *C
CC********************************************************************CC
 
      SUBROUTINE UCROSD(X,Y,Z)
 
      IMPLICIT   REAL*8 ( A-H, O-Z )
      REAL*8     X(3), Y(3), Z(3), A(3)
      A(1) = X(2)*Y(3) - X(3)*Y(2)
      A(2) = X(3)*Y(1) - X(1)*Y(3)
      A(3) = X(1)*Y(2) - X(2)*Y(1)
      CALL UVCOPY(6,A,Z)
      RETURN
      END
