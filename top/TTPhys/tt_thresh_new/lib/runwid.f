CC**********************************************************************
C*
C*  Real*8 function  RUNWID(E,P)
C*
CC**********************************************************************
 
      DOUBLE PRECISION FUNCTION RUNWID(E,P)
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      REAL   *8       E, P
C--
      COMMON /TTPARM/ MT, ALFS, ALF, S2W, YCUT, RWT, GAMMAT, ALAMB,
     .                CF, MB, MZ, MW, EFFMT, TLETA1, TLETA2, KAPPA,
     .                PI, IMAG
      REAL   *8       MT, ALFS, ALF, S2W, YCUT, RWT, GAMMAT, ALAMB,
     .                CF, MB, MZ, MW, TLETA1, TLETA2, KAPPA,
     .                PI
      COMPLEX*16      EFFMT, IMAG
C
C===========< Entry Point >=============================================
C
C--
C  Initialization.
C--
      RUNWID = 2*GAMMAT*( 1 + TLETA1*(E/MT) - TLETA2*(P/MT)**2 )
C--
C  That's it.
C--
      RETURN
      END
