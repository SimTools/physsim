CC**********************************************************************
C*      REAL*8  FUNCTION GAMMA0(ENERGY)
CC**********************************************************************
 
      DOUBLE PRECISION FUNCTION GAMMA0( ENERGY )
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      REAL   *8      ENERGY
C--
      COMMON /TTPARM/ MT, ALFS, ALF, S2W, YCUT, RWT, GAMMAT, ALAMB,
     .                CF, MB, MZ, MW, EFFMT, TLETA1, TLETA2, KAPPA,
     .                PI, IMAG
      REAL   *8       MT, ALFS, ALF, S2W, YCUT, RWT, GAMMAT, ALAMB,
     .                CF, MB, MZ, MW, TLETA1, TLETA2, KAPPA,
     .                PI
      COMPLEX*16      EFFMT, IMAG
C
C========< Entry Point >================================================
C
C--
C  Calculate the lowest Gamma_t at p = 0.
C--
      GAMMA0 = 2*GAMMAT*( 1 + TLETA1*ENERGY/MT )
C--
C  That's it.
C--
      RETURN
      END
