C* (Update Record)
C*   08/11/30  K.Fujii		Modified for BASES V5.1 in FULL.
C*

      SUBROUTINE USERIN
 
      IMPLICIT REAL*8 (A-H,O-Z)
C--
#include "ZHHCNSF.inc"
#include "USRPRM.inc"
#include "USMPRM.inc"
#include "BEMCNS.inc"
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
      SGEBM  = SGMEBM
      POLE   = POLEBM
      ISRBM  = ISRB
      IMD1MN = IMD1LO
      IMD1MX = IMD1HI
      PRINT *, ' ******< Job constants >*********** '
      PRINT *, '  '
      PRINT *, '     ROOTS  = ', ROOTS
      PRINT *, '     POLE   = ', POLE
      PRINT *, '     SGEBM  = ', SGEBM
      PRINT *, '     ISRBM  = ', ISRBM
      PRINT *, '     IMD1MN = ', IMD1MN
      PRINT *, '     IMD1MX = ', IMD1MX
      PRINT *, '     IMD2MN = ', IMD2MN
      PRINT *, '     IMD2MX = ', IMD2MX
      PRINT *, '  '
      PRINT *, ' ********************************** '
C--
      IF ( ROOTS.LT.AMSH+AMSZ ) THEN
         PRINT *, ' ROOTS is too small. '
         PRINT *, '    AMH = ', AMSH
         PRINT *, '    AMZ = ', AMSZ
      ENDIF
C--
C  That's it.
C--
      RETURN
      END
