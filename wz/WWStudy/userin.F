C* (Update Record)
C*  95/03/12  K.Fujii	Converted to unix. This version
C*			uses the new INSMCN in xcxc_lib.
C*  95/03/21  K.Fujii	Now calls SPPARM.
C*  99/05/17  K.Fujii	Modified for JSF.
C*

      SUBROUTINE USERIN
  
      IMPLICIT REAL*8 (A-H,O-Z)
C--
#include "WWCNSF.inc"
#include "USRPRM.inc"
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
C--
C  Set user constants.
C--
      ROOTS  = SQRTS
      POLE   = POLEBM
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
      PRINT *, '     ISRBM  = ', ISRBM
      PRINT *, '     IMD1MN = ', IMD1MN
      PRINT *, '     IMD1MX = ', IMD1MX      
      PRINT *, '     IMD2MN = ', IMD2MN
      PRINT *, '     IMD2MX = ', IMD2MX      
      PRINT *, '  '
      PRINT *, ' ********************************** '
C--
      IF ( ROOTS.LT.AMSW+AMSW ) THEN
         PRINT *, ' ROOTS is too small. '
         PRINT *, '    AMW = ', AMSW
         PRINT *, '    AMZ = ', AMSZ
      ENDIF
C--
C  That's it.
C--
      RETURN
      END
