      SUBROUTINE INSSCT(ALF,ALFS,WM,ZM,AM0,AMU,AM2,TANB,AMA,ASFT,ITAU)
 
      IMPLICIT    REAL*4 ( A-H, O-Z )
#include "SSCONS.inc"
#include "SSPTAB.inc"
#include "SSCUPL.inc"
      INTEGER*4   ITAU
      REAL   *4   ALF, ALFS, WM, ZM, AM0, AMU, AM2, TANB, AMA
C
C========< Entry Point >================================================
C
C--
C  Initialize SUSY parameters.
C--
      CALL INSSCN(ALF,ALFS,WM,ZM,AM0,AMU,AM2,TANB,AMA,ASFT)
C--
C  When ITAU is nonzero, reset sfermion and LSP masses.
C--
      IF ( ITAU.NE.0 )THEN
         SFM(1) = AM0
         SFM(2) = AM0
         SFM(3) = AM0
         SZM(1) = AMU
      ENDIF
C--
C  Print out modified parameters.
C--
      PRINT *, '             '
      PRINT *, ' ****** INSSCT MODIFIED /SSCONS/ AND /SSPTAB/ ********'
      PRINT *, '             '
      PRINT *, '    SFM(1) = ', SFM(1)
      PRINT *, '       (2) = ', SFM(2)
      PRINT *, '       (3) = ', SFM(3)
      PRINT *, '    SZM(1) = ', SZM(1)
      PRINT *, '        '
      PRINT *, ' *****************************************************'
C--
C  That's it.
C--
      RETURN
      END







