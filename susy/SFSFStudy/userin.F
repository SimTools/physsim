C* (Update Record)
C*   95/02/12  K.Fujii		Original version on UNIX.
C*   95/03/21  K.Fujii		Modified for BASES V5.1 in FULL.
C*   95/05/05  K.Fujii		Special treatment for nonzero IDOTAU.
C*   99/05/16  K.Fujii		Modified for JSF.
C*

      SUBROUTINE USERIN
 
      IMPLICIT REAL*8 (A-H,O-Z)
C--
C  User's commons.
C--
#include "SFSFCN.inc"
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
C  Initialize final state handedness.
C     IGENSF = generation #.
C     IHANDM = (1,2) for (L,R) for sf-.
C     IHANDP = (1,2) for (L,R) for sf+.
C--
      IGENSF = IGNSF
      IHANDM = IHNDM
      IHANDP = IHNDP
      IDOTAU = IDOTU
      HTAUM  = HTUM
      HTAUP  = HTUP
C--
C  Initialize numerical constants and particle table.
C--
      ALF   = 1/ALFI
      CALL INSMCN(ALF,ALFS,AMSW,AMSZ,AMSH,AMST)
      CALL INSSCT(ALF,ALFS,AMSW,AMSZ,AM0,AMU,AM2,TANB,AMA,ASFT,IDOTAU)
C--
C  Set user constants.
C--
      ROOTS  = SQRTS
      SGEBM  = SGMEBM
      POLE   = POLEBM
      ISRBM  = ISRB
C--      
      PRINT *, ' ******< Job constants >*********** '
      PRINT *, '  '
      PRINT *, '     ROOTS  = ', ROOTS
      PRINT *, '     SGEBM  = ', SGEBM
      PRINT *, '     POLE   = ', POLE
      PRINT *, '     ISRBM  = ', ISRB
      PRINT *, '  '
      PRINT *, '     IGENSF = ', IGENSF
      PRINT *, '     IHANDM = ', IHANDM
      PRINT *, '     IHANDP = ', IHANDP
      PRINT *, '     IDOTAU = ', IDOTAU
      PRINT *, '     HTAUM  = ', HTAUM
      PRINT *, '     HTAUP  = ', HTAUP
      PRINT *, '  '
      PRINT *, ' ********************************** '
C--
C  That's it.
C--
      RETURN
      END
