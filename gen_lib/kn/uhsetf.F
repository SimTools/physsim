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
C*       P (3)   : (R*4) : direction vector of new polar axis.
C*       EB(*,i) : (R*4) : i-th axis of the old reference frame.
C* (Output)
C*       EA(*,i) : (R*4) : i-th axis of the new reference frame.
C* (Relation)
C*    Calls UNRMV3, UNCROS, UVCOPY, UVZERO, and UABSV3.
C* (Update Record)
C*    6/22/90  K.Fujii      Derived from USETRF in 'T#PB.FORTLIB.FORT'.
C*    8/08/96  K.Fujii      Added protection for |P| = 0.
C*
CC**********************************************************************
 
      SUBROUTINE UHSETF(P,EB,EA)
 
      REAL*4     P(3), EB(3,3), EA(3,3)
      REAL*4     E(3,3), E2(3)
      DATA XMIN  / 1.E-6  /
      DATA ZERO  / 1.E-20 /
C
C========< Entry Point >================================================
C
C--
C  Set new Z-axis.
C--
      AP    = UABSV3(P)
      IF ( AP.LT.ZERO ) THEN
         CALL UVCOPY(9,EB,EA)
         RETURN
      ELSE
         CALL UNRMV3(P,E(1,3))
      ENDIF
C--
C  Set new Y-axis.
C--
      CALL UCROSS(EB(1,3),E(1,3),E2)
      CSTH  = UDOT3(E(1,3),EB(1,3))
      AE2   = UABSV3(E2)
C--
      IF ( AE2.LT.XMIN ) THEN
         CALL UVCOPY(9,EB,EA)
         IF ( CSTH.LT.0. ) THEN
            CALL USCLM3(-1.,EB(1,3),EA(1,3))
            CALL USCLM3(-1.,EB(1,1),EA(1,1))
         ENDIF
         RETURN
      ELSE
         CALL UNRMV3(E2,E(1,2))
      ENDIF
C--
C  Set new X-axis.
C--
      CALL UNCROS(E(1,2),E(1,3),E(1,1))
C--
C  Store E in EA.
C--
      CALL UVCOPY(9,E,EA)
C--
C  That's it.
C--
      RETURN
      END
