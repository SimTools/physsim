C*
C* (Update Record)
C* 2010/11/29  K.Fujii	Derived from TTStudy/userin.F
C*

      SUBROUTINE USERIN
 
      IMPLICIT REAL*8 (A-H,O-Z)
C--
      INCLUDE 'TBWCNF.inc'
      INCLUDE 'USRPRM.inc'
#include "BEMCNS.inc"
#include "USMPRM.inc"
#include "BSHUFL.inc"
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
      POLE   = POLEBM
      ISRBM  = ISRB
      IMD1MN = IMD1LO
      IMD1MX = IMD1HI
      IMD2MN = IMD2LO
      IMD2MX = IMD2HI
C--      
      PRINT *, ' ******< Job constants >*********** '
      PRINT *, '  '
      PRINT *, '     ROOTS  = ', ROOTS
      PRINT *, '     POLE   = ', POLE
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
