 
      SUBROUTINE GRQDSG(MODE,E,Q,DSG)
 
      IMPLICIT REAL*8  ( A-H, O-Z )
      INTEGER*4       MODE
      REAL*8          E, Q, DSG
C--
C  Parameters.
C--
      COMMON /TTPARM/ MT, ALFS, ALF, S2W, YCUT, RWT, GAMMAT, ALAMB,
     .                CF, MB, MZ, MW, EFFMT, TLETA1, TLETA2, KAPPA,
     .                PI, IMAG
      REAL   *8       MT, ALFS, ALF, S2W, YCUT, RWT, GAMMAT, ALAMB,
     .                CF, MB, MZ, MW, TLETA1, TLETA2, KAPPA,
     .                PI
      COMPLEX*16      EFFMT, IMAG
C--
      COMMON /RKPARM/ DRS, DRP, EPSS, EPSP
      REAL   *8       DRS, DRP, EPSS, EPSP
C--
C  Other variables.
C--
      INTEGER*4       STEPS1, STEPS2, STEPP1, STEPP2
      SAVE
C
C===========< Entry Point >=============================================
C
C--
C  Check if Geen function calculation is necessary.
C--
      IF ( MODE.LE.3 ) THEN
C--
C  Calculate S-wave Green's function.
C--
         CALL  QCDS1(E,DRS,EPSS,STEPS1)
         CALL  QCDS2(STEPS2)
         CALL  QCDS3
C--
C  Calculate P-wave Green's function.
C--
         CALL  QCDP1(E,DRP,EPSP,STEPP1)
         CALL  QCDP2(STEPP2)
         CALL  QCDP3
      ENDIF
C--
C  Calculate momentum distribution from Green's function.
C--
      CALL  GRQPTP(MODE,E,Q,DSG)
C--
C  Adjust loop control vairables for the next step.
C--
      IF ( MODE.LE.3 ) THEN
         CALL  GRQADJ(MODE,STEPS1,STEPS2,STEPP1,STEPP2)
      ENDIF
C--
C  That's it.
C--
      RETURN
      END
