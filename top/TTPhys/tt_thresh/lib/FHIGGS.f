C*  94/03/30  K.Fujii        New Higgs treatment according to
C*                           PLB 316 (1993) 360.
C*
C*
      SUBROUTINE FHIGGS(ALF,S2W,AMW,AMT,AMH,BTH,FACH)
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      INTEGER*4 MODE
      REAL*8    ALF, S2W, AMW, AMT, AMH, BTH, FACH
      DATA NCALL /0/
      SAVE
C
C========< Entry Point >================================================
C
C--
C  Initialization.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL = 1
         xPI   = ACOS(-1.D0)
         x4PI  = 4*xPI
      ENDIF
C--
C  Calculate L4.
C--
      R = (AMH/AMT)**2
      IF ( R.LE.4.D0 ) THEN
         EL4 = SQRT( R*(4-R) ) * ACOS( SQRT(R)/2 )
      ELSE
         EL4 = -SQRT(R*(R-4))*LOG((1+SQRT(1-4/R))/(1-SQRT(1-4/R)))/2
      ENDIF
C--
C  Calculate FTH.
C--
      FTH = -1.D0/12 * ( -12 + 4*R + (-12+9*R-2*R*R)*LOG(R)
     .                           + (2/R)*(-6+5*R-2*R*R)*EL4 )
C--
C  Calculate FACH.
C--
      FACH = 1 + BTH**2*ALF/(x4PI*S2W) * (AMT/AMW)**2
     .           * ( FTH - xPI*(AMT/AMH) )
      FACH = FACH*FACH
C--
C  That's it.
C--
      RETURN
      END
