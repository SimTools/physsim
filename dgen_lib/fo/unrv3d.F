CC********************************************************************CC
C*                                                                    *C
C* ====================---===                                         *C
C*  Subroutine UNRV3D(X,UX)                                           *C
C* ====================---===                                         *C
C*                                                                    *C
C*    1) Purpose                                                      *C
C*        Normalizes the given vector,X.                              *C
C*    2) I/O specification                                            *C
C*        Input                                                       *C
C*          X(3)   ; 3-dim. vector                                    *C
C*        Output                                                      *C
C*          UX(3)  ; unit vector in the direction of X.               *C
C*    3) Relation                                                     *C
C*        Calls                                                       *C
C*            UABV3D, USCM3D                                          *C
C*    4) Update record                                                *C
C*         6/12/85 TKSF Original version.                             *C
C*                                                                    *C
CC********************************************************************CC
 
      SUBROUTINE UNRV3D(X,UX)
 
      IMPLICIT   REAL*8 ( A-H, O-Z )
      REAL*8     X(3),UX(3)
 
      AX   = UABV3D(X)
      AX   = 1/AX
      CALL USCM3D(AX,X,UX)
      RETURN
      END
