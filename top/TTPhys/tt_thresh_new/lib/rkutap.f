CC*********************************************************************
C*
C*  SUBROUTINE RKUTAP(X,Y,Z,DX,DY,DZ)
C*
C*  This subroutine solves the differential equation
C*
C*     dy/dx = f_P(x,y,z),  dz/dx = g_P(x,y,z)
C*
C*  using the Lunge-Kutta method.  x is the real variable,
C*  while y, z, f_P, g_P are complex.
C*
C*  Inputs:  x, y, z, dx
C*  Outputs: dy, dz
C*
CC*********************************************************************
 
      SUBROUTINE RKUTAP(X,Y,Z,DX,DY,DZ)
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      REAL   *8   X, DX
      COMPLEX*16  Y, Z, DY, DZ
C--
      COMPLEX*16  D1, D2, D3, D4, E1, E2, E3, E4, FP, GP
      REAL   *8   HALFDX
C
C===========< Entry Point >=============================================
C
C--
C  Initialization.
C--
      HALFDX = .5D0*DX
C--
C  Calculate this step.
C--
      D1 = DX*FP( X, Y, Z )
      E1 = DX*GP( X, Y, Z )
C--
      D2 = DX*FP( X+HALFDX, Y+.5D0*D1, Z+.5D0*E1 )
      E2 = DX*GP( X+HALFDX, Y+.5D0*D1, Z+.5D0*E1 )
C--
      D3 = DX*FP( X+HALFDX, Y+.5D0*D2, Z+.5D0*E2 )
      E3 = DX*GP( X+HALFDX, Y+.5D0*D2, Z+.5D0*E2 )
C--
      D4 = DX*FP( X+DX, Y+D3, Z+E3 )
      E4 = DX*GP( X+DX, Y+D3, Z+E3 )
C--
      DY = ( D1 + 2*D2 + 2*D3 + D4 )/6
      DZ = ( E1 + 2*E2 + 2*E3 + E4 )/6
C--
C  That's it.
C--
      RETURN
      END
 
CC*********************************************************************
C*  COMPLEX*16 FUNCTION  FP(X,Y,Z)
CC*********************************************************************
 
      DOUBLE COMPLEX FUNCTION FP(X,Y,Z)
 
      REAL   *8   X
      COMPLEX*16  Y, Z
C
C===========< Entry Point >=============================================
C
C--
C  Calculate FP.
C--
      FP = Z
C--
C  That's it.
C--
      RETURN
      END
 
CC*********************************************************************
C*
C*  COMPLEX*16 FUNCTION  GP(X,Y,Z)
C*
C*  Inputs:
C*       x, y, z,
C*       ener : E = sqrt(s) - 2m_t ; stored in  common/ener/...
C*       mt, Gammat, imag;  stored in  common/param/...
C*  Output:
C*       GP = ( - eff_mt*( E + i*Gamma0/2 - V(x) ) + 2/x**2 )*y
C*
CC*********************************************************************
 
      DOUBLE COMPLEX FUNCTION GP(X,Y,Z)
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      REAL   *8       X
      COMPLEX*16      Y, Z
C--
      COMMON /ENER/   ENER
      REAL   *8       ENER
C--
      COMMON /TTPARM/ MT, ALFS, ALF, S2W, YCUT, RWT, GAMMAT, ALAMB,
     .                CF, MB, MZ, MW, EFFMT, TLETA1, TLETA2, KAPPA,
     .                PI, IMAG
      REAL   *8       MT, ALFS, ALF, S2W, YCUT, RWT, GAMMAT, ALAMB,
     .                CF, MB, MZ, MW, TLETA1, TLETA2, KAPPA,
     .                PI
      COMPLEX*16      EFFMT, IMAG
C--
      REAL   *8       V, GAMMA0
C
C===========< Entry Point >=============================================
C
C--
C  Calculate GP.
C--
      GP = ( - EFFMT*( ENER + IMAG*GAMMA0(ENER)/2 - V(X) )
     .                                          + 2/X**2 )*Y
C--
C  That's it.
C--
      RETURN
      END
