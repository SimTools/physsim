C***********************************************************************
C*
C*============================---==
C* Subroutine USVCKM(MODE,I,J,VFF)
C*============================---==
C*
C* (Ppurpose)
C*    This subroutine allows users to put or get a CKM matrix element
C*    in or from /SMPTAB/ without directly referring to it.
C* (Inputs)
C*      MODE     : (R*4) : (1,2)=(get,put).
C*      I        : (I*4) : (1,2,3)=(u,c,t).
C*      J        : (I*4) : (1,2,3)=(d,s,b).
C* (Input/Output)
C*      VFF      : (R*4) : V_ij
C* (Relation)
C*    Calls no subroutines.
C* (Update Record)
C*   95/06/07  K.Fujii	Original version.
C*
C***********************************************************************
 
      SUBROUTINE USVCKM(MODE,I,J,VFF)
 
      INTEGER  *4   MODE, I, J
      REAL     *4   VFF
C--
#include "../sm/include/SMPTAB.inc"
C
C========< Entry Point >================================================
C
C--
C  Branch on MODE.
C--
      IF ( I.LT.1 .OR. I.GT.3 .OR.
     .     J.LT.1 .OR. J.GT.3 ) THEN
         PRINT *, ' >>>> ERROR in USVCKM >>>>'
         PRINT *, '   Invalid I or J: (I,J) = (', I, ',', J,')'
         PRINT *, '   STOP immediately. '
         STOP
      ELSE
         IF ( MODE.EQ.1 ) THEN
            VFF = VKM(I,J,2)
         ELSE
            VKM(I,J,2) = VFF
         ENDIF
      ENDIF
C--
C  That's it.
C--
      RETURN
      END
