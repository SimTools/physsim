CC********************************************************************CC
C*                                                                    *C
C* ========================--==                                       *C
C*  Subroutine UBTRND(PB,E,PA)                                        *C
C* ========================--==                                       *C
C*                                                                    *C
C*    1) Purpose                                                      *C
C*        Rotates PB to PA so that PB in the frame E is PA            *C
C*        in the standard frame.                                      *C
C*    2) I/O specification                                            *C
C*        Input                                                       *C
C*          PB(3)  ; 3-dim. vector                                    *C
C*          E(3,3) ; orthonormal basis of the frame in which          *C
C*                   PB is defined.                                   *C
C*        Output                                                      *C
C*          PA(3)  ; PB in the standard basis.                        *C
C*    3) Relation                                                     *C
C*        Calls                                                       *C
C*          UVCOPY                                                    *C
C*    4) Update record                                                *C
C*         6/12/85 TKSF Original version.                             *C
C*                                                                    *C
CC********************************************************************CC
 
      SUBROUTINE UBTRND(PB,E,PA)
 
      IMPLICIT   REAL*8 ( A-H, O-Z )
      REAL*8     PB(3), PA(3), P(3), E(3,3)
      DO 10 I = 1, 3
         P(I) = 0
      DO 10 II = 1, 3
         P(I) = P(I) + PB(II)*E(I,II)
10    CONTINUE
      CALL UVCOPY(6,P,PA)
      RETURN
      END
