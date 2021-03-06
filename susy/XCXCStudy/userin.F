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
#include "USSPRMP.inc"
#include "USMPRM.inc"
#include "BEMCNS.inc"
#include "BSHUFL.inc"
#include "BSITPR.inc"
C--
C  Local R*4 variables used in subroutine calls.
C--
      REAL*4        ALF
      integer*4     is
C
C========< Entry Point >================================================
C
C--
C  Initialize numerical constants and particle table.
C--
      ALF   = 1/ALFI
      CALL INSMCN(ALF,ALFS,AMSW,AMSZ,AMSH,AMST)
c      CALL INSSCN(ALF,ALFS,AMSW,AMSZ,AM0,AMU,AM2,TANB,AMA,ASFT)
      if (isgut.eq.1) then
         is = 0
      else
         is = 3
      endif
      CALL INSSCNP(is,ALF,ALFS,AMSW,AMSZ,AM0,AMU,
     .                am1,AM2,am3,TANB,AMA,ASFT)
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
      IMD1MN = IMD1LO
      IMD1MX = IMD1HI
      IMD2MN = IMD2LO
      IMD2MX = IMD2HI
      PRINT *, ' ******< Job constants >*********** '
      PRINT *, '  '
      PRINT *, '     ROOTS  = ', ROOTS
      PRINT *, '     SGEBM  = ', SGEBM
      PRINT *, '     POLE   = ', POLE
      PRINT *, '     GAMSW1 = ', GAMSW1
      PRINT *, '     ISRBM  = ', ISRBM
      PRINT *, '     IMD1MN = ', IMD1MN
      PRINT *, '     IMD1MX = ', IMD1MX      
      PRINT *, '     IMD2MN = ', IMD2MN
      PRINT *, '     IMD2MX = ', IMD2MX
      PRINT *, '  '
      PRINT *, ' ********************************** '
C--
C  That's it.
C--
      RETURN
      END
