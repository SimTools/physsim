CC********************************************************************CC
C*                                                                    *C
C* ====================---===                                         *C
C*  Subroutine USETRD(P,E)                                            *C
C* ====================---===                                         *C
C*                                                                    *C
C*    1) Purpose                                                      *C
C*        Setups a reference frame with E(*,3) along P and E(*,1)     *C
C*        perpendicular to P-E3(*) plane.                             *C
C*        (i.e. E(*,1) in E1-E2 plane)                                *C
C*    2) I/O specification                                            *C
C*        Input                                                       *C
C*          P(3)   ; direction vector of polar axis.                  *C
C*        Output                                                      *C
C*          E(*,I) ; I-th axis of the reference frame                 *C
C*    3) Relation                                                     *C
C*        Calls                                                       *C
C*          UNRV3D, UNCRSD                                            *C
C*    4) Update record                                                *C
C*         6/12/85 TKSF Original version.                             *C
C*         2/27/90 TKSF Modified to avoid round off error.            *C
C*                                                                    *C
CC********************************************************************CC
 
      SUBROUTINE USETRD(P,E)
 
      IMPLICIT   REAL*8 ( A-H, O-Z )
      REAL*8     E(3,3), E3(3), E2(3), P(3)
C
      E3(1) = 0
      E3(2) = 0
      E3(3) = 1
      CALL UNRV3D(P,E(1,3))
      E2(1) = E(1,3)
      E2(2) = E(2,3)
      E2(3) = 0
      AE2   = SQRT( E2(1)*E2(1) + E2(2)*E2(2) )
      IF ( AE2.LT.1.D-8 ) THEN
         E(1,1) = 1
         E(2,1) = 0
         E(3,1) = 0
      ELSE
         E2(1) = E2(1)/AE2
         E2(2) = E2(2)/AE2
         CALL UNCRSD(E2(1),E3(1),E(1,1))
      ENDIF
      CALL UNCRSD(E(1,3),E(1,1),E(1,2))
      RETURN
      END
