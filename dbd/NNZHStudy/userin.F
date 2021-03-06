C* (Update Record)
C*   95/04/22  K.Fujii		Original version for NNZH.
C*   99/05/20  K.Fujii		Modified for JFS.
C*

      SUBROUTINE USERIN
 
      IMPLICIT REAL*8 (A-H,O-Z)
C--
#include "NNZHCNF.inc"
#include "USRPRM.inc"
#include "USMPRM.inc"
#include "BEMCNS.inc"
#include "BSHUFL.inc"
#include "BSITPR.inc"
C--
C  Local R*4 variables used in subroutine calls.
C--
      REAL*8        ALF
C
C========< Entry Point >================================================
C
C--
C  Initialize numerical constants and particle table.
C--
      ALF   = 1/ALFI
      CALL INSMCN(DBLE(ALF),DBLE(ALFS),DBLE(AMSW),DBLE(AMSZ),
     .            DBLE(AMSH),DBLE(AMST))
C--
C  Set user constants.
C--
      ROOTS  = SQRTS
      SGEBM  = SGMEBM
      POLE   = POLEBM
      ISRBM  = ISRB
      IMDMN  = IMDMIN
      IMDMX  = IMDMAX
C--
      PRINT *, ' ******< Job constants >*********** '
      PRINT *, '  '
      PRINT *, '     ROOTS  = ', ROOTS
      PRINT *, '     SGEBM  = ', SGEBM
      PRINT *, '     POLE   = ', POLE
      PRINT *, '     ISRBM  = ', ISRBM
#ifdef WITH_DBD_STANDARD
      PRINT *, ' ( Use Tims lumi function, if ISRBM > 100. ',
     >         ' mode=Mod(ISRBM,100)'
#endif
      PRINT *, '     IMDMN  = ', IMDMN
      PRINT *, '     IMDMX  = ', IMDMX
      PRINT *, '  '
      PRINT *, ' ********************************** '
C--
C  That's it.
C--
      RETURN
      END
