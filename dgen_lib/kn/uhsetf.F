CC**********************************************************************
C*
C*========================--==
C* Subroutine UHSETF(P,EB,EA)
C*========================--==
C*
C* (Purpose)
C*    Setups a reference frame with EA(*,3) along P and EA(*,2)
C*    perpendicular to P-EB(*,3) plane.
C* (Input)
C*       P (3)   : (R*8) : direction vector of new polar axis.
C*       EB(*,i) : (R*8) : i-th axis of the old reference frame.
C* (Output)
C*       EA(*,i) : (R*8) : i-th axis of the new reference frame.
C* (Relation)
C*    Calls UNRV3D, UNCRSD, UVCOPY, UVZERO, and UABV3D.
C* (Update Record)
C*    6/22/90  K.Fujii      Derived from USETRF in 'T#PB.FORTLIB.FORT'.
C*
CC**********************************************************************
 
      SUBROUTINE UHSETF(P,EB,EA)
 
      REAL*8     P(3), EB(3,3), EA(3,3)
      REAL*8     E(3,3), E2(3)
      DATA XMIN  / 1.D-12 /
C
C========< Entry Point >================================================
C
C--
C  Set new Z-axis.
C--
      CALL UNRV3D(P,E(1,3))
C--
C  Set new Y-axis.
C--
      CALL UCROSD(EB(1,3),E(1,3),E2)
      AE2   = UABV3D(E2)
      CSTH  = UDOT3D(E(1,3),EB(1,3))
C--
      IF ( AE2.LT.XMIN ) THEN
         CALL UVCOPY(18,EB,EA)
         IF ( CSTH.LT.0.D0 ) THEN
            CALL USCM3D(-1.D0,EB(1,3),EA(1,3))         
         ENDIF
         RETURN
      ELSE
         CALL UNRV3D(E2,E(1,2))
      ENDIF
C--
C  Set new X-axis.
C--
      CALL UNCRSD(E(1,2),E(1,3),E(1,1))
C--
C  Store E in EA.
C--
      CALL UVCOPY(18,E,EA)
C--
C  That's it.
C--
      RETURN
      END
