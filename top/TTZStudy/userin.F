C*
C* (Update Record)
C* 95/03/15  K.Fujii	Modified for BASES V5.1 in FULL.
C* 95/03/16  K.Fujii	Takes into account beam and bremsstrahlng
C*			as well as ISR and QED vtx correction:
C*                      ( 1+3*BTE/4 ) and
C*                      ( 1 + 2*ALF0/PI * ( PI**2/6 - 1/4 ) )
C*                      included.
C* 95/03/21  K.Fujii	Now calls SPPARM.
C* 95/06/07  K.Fujii	Added ISRB switch.
C* 99/05/21  K.Fujii	Modified for JSF.
C*

      SUBROUTINE USERIN
 
      IMPLICIT REAL*8 (A-H,O-Z)
C--
#include "USRPRM.inc"
#include "TTZCNF.inc"
#include "BEMCNS.inc"
#include "USMPRM.inc"
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
      SGEBM  = SGMEBM
      ISRBM  = ISRB
      IMDZMN = IMDZLO
      IMDZMX = IMDZHI
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
      PRINT *, '     IMDZMN = ', IMDZMN
      PRINT *, '     IMDZMX = ', IMDZMX      
      PRINT *, '     IMD1MN = ', IMD1MN
      PRINT *, '     IMD1MX = ', IMD1MX      
      PRINT *, '     IMD2MN = ', IMD2MN
      PRINT *, '     IMD2MX = ', IMD2MX      
      PRINT *, '  '
      PRINT *, ' ********************************** '
C--
      IF ( ROOTS.LT.AMST+AMST+AMSZ ) THEN
         PRINT *, ' ROOTS is too small. '
         PRINT *, '    AMT = ', AMST
         PRINT *, '    AMZ = ', AMSZ
      ENDIF
C--
C  That's it.
C--
      RETURN
      END
