C* (Update Record)
C*   95/03/21  K.Fujii		Modified for BASES V5.1 in FULL.
C*   95/04/04  K.Fujii		Added ISRB flag.
C*   04/11/23  K.Fujii		Modified for JSF.
C*

      SUBROUTINE USERIN
 
      IMPLICIT REAL*8 (A-H,O-Z)
C--
#include "XNXNCN.inc"
#include "USRPRM.inc"
#include "USSPRMP.inc"
#include "USMPRM.inc"
#include "BEMCNS.inc"
#include "BSHUFL.inc"
#include "BSITPR.inc"
C--
C  Local R*4 variables used in subroutine calls.
C--
      REAL*8        ALF, ASFT8(3)
      integer*4     is
      real*8        am1d, am2d, am3d
C
C========< Entry Point >================================================
C
C--
C  Initialize numerical constants and particle table.
C--
      ALF   = 1/ALFI
      CALL INSMCN(DBLE(ALF),DBLE(ALFS),DBLE(AMSW),DBLE(AMSZ),
     .            DBLE(AMSH),DBLE(AMST))
c     CALL INSSCN(DBLE(ALF),DBLE(ALFS),DBLE(AMSW),DBLE(AMSZ),
c    .            DBLE(AM0),DBLE(AMU),DBLE(AM2),DBLE(TANB),
c    .            DBLE(AMA),ASFT8)
      am2d = am2
      if (isgut.eq.1) then
        is = 0
      else
        is = 3
        am1d = am1
        am3d = am3
      endif
      CALL INSSCNP(is,DBLE(ALF),DBLE(ALFS),DBLE(AMSW),DBLE(AMSZ),
     .            DBLE(AM0),DBLE(AMU),am1d,am2d,am3d,DBLE(TANB),
     .            DBLE(AMA),ASFT8)
C--
      CALL SETGSZ
C--
      IF ( GAMSZ2.GT.0.0D0 ) THEN
         CALL MODGAM('SZ',2,GAMSZ2)
      ENDIF
C--
C  Set user constants.
C--
      ROOTS  = SQRTS
      SGEBM  = SGMEBM
      POLE   = POLEBM
      ISRBM  = ISRB
      IMD2MN = IMD2LO
      IMD2MX = IMD2HI
      PRINT *, ' ******< Job constants >*********** '
      PRINT *, '  '
      PRINT *, '     ROOTS  = ', ROOTS
      PRINT *, '     SGEBM  = ', SGEBM
      PRINT *, '     POLE   = ', POLE
      PRINT *, '     GAMSZ2 = ', GAMSZ2
      PRINT *, '     ISRBM  = ', ISRBM
      PRINT *, '     IMD2MN = ', IMD2MN
      PRINT *, '     IMD2MX = ', IMD2MX
      PRINT *, '  '
      PRINT *, ' ********************************** '
C--
C  That's it.
C--
      RETURN
      END
