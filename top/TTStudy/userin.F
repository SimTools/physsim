C*
C* (Update Record)
C* 95/03/16  K.Fujii	Modified for BASES V5.1 in FULL.
C* 95/03/21  K.Fujii	Now calls SPPARM.
C* 95/06/07  K.Fujii	Added ISRB switch.
C* 98/12/20  K.Fujii	Modified for JSF.
C*

      SUBROUTINE USERIN
 
      IMPLICIT REAL*8 (A-H,O-Z)
C--
      INCLUDE 'TTZCNF.inc'
      INCLUDE 'THRCNS.inc'
      INCLUDE 'USRPRM.inc'
#include "BEMCNS.inc"
#include "USMPRM.inc"
#include "BSHUFL.inc"
C--
C  Local R*4 variables used in subroutine calls.
C--
      REAL*4        ALF, VFF
C
C========< Entry Point >================================================
C
C--
C  Modify VKM.
C--
      CALL USVCKM(1,3,3,VFF)
      VFF = VFF*VKMT
      CALL USVCKM(2,3,3,VFF)
C--
C  Initialize numerical constants and particle table.
C--
      ALF   = 1/ALFI
      CALL INSMCN(ALF,ALFS,AMSW,AMSZ,AMSH,AMST)
C--
C  Set user constants.
C--
      ROOTS  = 2*AMST + DELTRS
      SGEBM  = SGMEBM
      POLE   = POLEBM
      ISRBM  = ISRB
      IMD1MN = IMD1LO
      IMD1MX = IMD1HI
      IMD2MN = IMD2LO
      IMD2MX = IMD2HI
      VKMTT  = VKMT
      BETAH  = BETH
      NRQCD  = NQCD
C--      
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
      PRINT *, '     VKMTT  = ', VKMTT
      PRINT *, '     BETAH  = ', BETAH
      PRINT *, '     NRQCD  = ', NRQCD
      PRINT *, '  '
      PRINT *, ' ********************************** '
C--
C  That's it.
C--
      RETURN
      END
