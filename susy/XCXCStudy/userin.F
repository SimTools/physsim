C* (Update Record)
C*   95/03/21  K.Fujii		Modified for BASES V5.1 in FULL.
C*   95/04/04  K.Fujii		Added ISRB flag.
C*   99/05/22  K.Fujii		Modified for JSF.
C*

      SUBROUTINE USERIN
 
      IMPLICIT REAL*8 (A-H,O-Z)
C--
#include "XCXCCN.inc"
#include "USRPRM.inc"
#include "USSPRM.inc"
#include "USMPRM.inc"
#include "BEMCNS.inc"
#include "BSHUFL.inc"
#include "BSITPR.inc"
C--
C  Local R*4 variables used in subroutine calls.
C--
      REAL*4        ALF
C
C========< Entry Point >================================================
C
C--
C  Initialize numerical constants and particle table.
C--
      ALF   = 1/ALFI
      CALL INSMCN(ALF,ALFS,AMSW,AMSZ,AMSH,AMST)
      CALL INSSCN(ALF,ALFS,AMSW,AMSZ,AM0,AMU,AM2,TANB,AMA)
C--
      IF ( GAMSW1.GT.0.0 ) THEN
         CALL MODGAM('SW',1,GAMSW1)
      ENDIF
C--
C  Set user constants.
C--
      ROOTS  = SQRTS
      SGEBM  = SGMEBM
      POLE   = POLEBM
      ISRBM  = ISRB
      PRINT *, ' ******< Job constants >*********** '
      PRINT *, '  '
      PRINT *, '     ROOTS  = ', ROOTS
      PRINT *, '     SGEBM  = ', SGEBM
      PRINT *, '     POLE   = ', POLE
      PRINT *, '     GAMSW1 = ', GAMSW1
      PRINT *, '     ISRBM  = ', ISRBM
      PRINT *, '  '
      PRINT *, ' ********************************** '
C--
C  That's it.
C--
      RETURN
      END
