CC**********************************************************************
C*
C*  SUBROUTINE GTQCDP( RQCD )
C*
C*  This subroutine prepares the parameters related to the QCD
C*  potential, which will be used in defining the boundary conditions
C*  of the S-wave Green's function, and used in calculating the
C*         ImG(0;E)!p<Lambda
C*
C*  Input:
C*       RQCD
C*       pi, CF, EFFMT, Lambda; stored in common/param/...
C*
C*  Output:    stored in common/optbcd/...
C*       b0, b1, beta,
C*       opt = 2Lambda/pi + beta*( ff1 + b1/b0^2*ff2 )
C*
CC**********************************************************************
 
      SUBROUTINE  GTQCDP( RQCD )
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      REAL   *8       RQCD
C--
      COMMON /TTPARM/ MT, ALFS, ALF, S2W, YCUT, RWT, GAMMAT, ALAMB,
     .                CF, MB, MZ, MW, EFFMT, TLETA1, TLETA2, KAPPA,
     .                PI, IMAG
      REAL   *8       MT, ALFS, ALF, S2W, YCUT, RWT, GAMMAT, ALAMB,
     .                CF, MB, MZ, MW, TLETA1, TLETA2, KAPPA,
     .                PI
      COMPLEX*16      EFFMT, IMAG
C--
      COMMON /OPTBCD/ B0, B1, BETA, OPT
      REAL   *8       B0, B1
      COMPLEX*16      BETA, OPT
C--
      REAL   *8       FX
C
C========< Entry Point >================================================
C
C--
C  Initialize parameters.
C--
      B0   = 23.D0/6.D0
      B1   = 29.D0/6.D0
      BETA = PI*CF/B0 * EFFMT
      OPT  = 2*ALAMB/PI + BETA*FX( ALAMB*RQCD )
C--
C  That's it.
C--
      RETURN
      END
