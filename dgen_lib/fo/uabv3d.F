CC********************************************************************CC
C*                                                                    *C
C*===========------====                                               *C
C*  Function UABV3D(X)                                                *C
C*===========------====                                               *C
C*                                                                    *C
C*    1) Purpose                                                      *C
C*        Calculates the absolute value of X.                         *C
C*    2) I/O specification                                            *C
C*        Input                                                       *C
C*          X(3)   ; 3-dim. vector                                    *C
C*        Output                                                      *C
C*          UABV3D ; absolute value of X.                             *C
C*    3) Relation                                                     *C
C*        Calls                                                       *C
C*          UDOT3D                                                    *C
C*    4) Update record                                                *C
C*         6/12/85 TKSF Original version.                             *C
C*                                                                    *C
CC********************************************************************CC
 
      DOUBLE PRECISION FUNCTION UABV3D(X)
 
      IMPLICIT   REAL*8 ( A-H, O-Z )
      REAL*8     X(3)
 
      UABV3D = UDOT3D(X,X)
      UABV3D = SQRT(UABV3D)
      RETURN
      END
