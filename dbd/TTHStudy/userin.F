C* (Update Record)
C*  95/03/28  K.Fujii	Original version for TTH.
C*  95/06/07  K.Fujii	Added ISRB switch.
C*  99/05/21  K.Fujii	Modified for JSF.
C*

      SUBROUTINE USERIN
 
      IMPLICIT REAL*8 (A-H,O-Z)
C--
#include "TTZCNF.inc"
#include "USRPRM.inc"
#include "USMPRM.inc"
#include "BEMCNS.inc"
#include "BSHUFL.inc"
#include "BSITPR.inc"
C--
C  Local R*4 variables used in subroutine calls.
C--
      REAL*4   ALF
C
C========< Entry Point >================================================
C
C--
C  Initialize numerical constants and particle table.
C--
      ALF   = 1/ALFI
      CALL INSMCN(ALF,ALFS,AMSW,AMSZ,AMSH,AMST)
C--
C  Set user constants.
C--
      ROOTS  = SQRTS
      POLE   = POLEBM
      POLP   = POLPBM
      SGEBM  = SGMEBM
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
      PRINT *, '     POLP   = ', POLP
      PRINT *, '     ISRBM  = ', ISRBM
#ifdef WITH_DBD_STANDARD
      PRINT *, ' ( Use Tims lumi function, if ISRBM > 100. ',
     >         ' mode=Mod(ISRBM,100)'
#endif
      PRINT *, '     IMD1MN = ', IMD1MN
      PRINT *, '     IMD1MX = ', IMD1MX      
      PRINT *, '     IMD2MN = ', IMD2MN
      PRINT *, '     IMD2MX = ', IMD2MX      
#ifdef WITH_DBD_STANDARD
      PRINT *, '     NDCYSL = ', NDCYSL, '(W decay mode to be selected)'
#endif
      PRINT *, '  '
      PRINT *, ' ********************************** '
C--
C  That's it.
C--
      RETURN
      END
