 
      SUBROUTINE GRQADJ(MODE,STEPS1,STEPS2,STEPP1,STEPP2)
 
      IMPLICIT REAL*8  ( A-H, O-Z )
      INTEGER*4       MODE, STEPS1, STEPS2, STEPP1, STEPP2
C--
C  Parameters.
C--
      COMMON /RKPARM/ DRS, DRP, EPSS, EPSP
      REAL*8          DRS, DRP, EPSS, EPSP
C--
      REAL*8          DRSOLD, DRSNEW, DRPOLD, DRPNEW,
     .                EPSNEW, EPPNEW,
     .                ADJ1, ADJ2
      INTEGER*4       STPS2O, STPP2O
C
C===========< Entry Point >=============================================
C
C--
C  Adjust integration steps and criteria.
C--
      IF ( MODE.LE.2 ) THEN
         DRSOLD = DRS
         DRPOLD = DRP
         STPS2O = STEPS2
         STPP2O = STEPP2
C--
         EPSNEW = ADJ2(EPSS,STEPS1,STEPS2)
         EPPNEW = ADJ2(EPSP,STEPP1,STEPP2)
C--
         EPSS   = EPSNEW
         EPSP   = EPPNEW
      ELSE
         DRSNEW = ADJ1(DRSOLD,STPS2O,DRS,STEPS2)
         DRPNEW = ADJ1(DRPOLD,STPP2O,DRP,STEPP2)
         EPSNEW = ADJ2(EPSS,STEPS1,STEPS2)
         EPPNEW = ADJ2(EPSP,STEPP1,STEPP2)
C--
         DRSOLD = DRS
         DRPOLD = DRP
         STPS2O = STEPS2
         STPP2O = STEPP2
C--
         DRS    = DRSNEW
         DRP    = DRPNEW
         EPSS   = EPSNEW
         EPSP   = EPPNEW
      END IF
C--
C  That's it.
C--
      RETURN
      END
