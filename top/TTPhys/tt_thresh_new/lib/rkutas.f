CC*********************************************************************
C*
C*  Subroutine RKUTAS(X,Y,Z,DX,DY,DZ)
C*
C*  This subroutine solves the differential equation
C*
C*     dy/dx = f_S(x,y,z),  dz/dx = g_S(x,y,z)
C*
C*  using the Runge-Kutta method.  x is the real variable,
C*  while y, z, f_S, g_S are complex.
C*
C*  Inputs:  x, y, z, dx
C*  Outputs: dy, dz
C*
CC*********************************************************************
 
      SUBROUTINE RKUTAS(X,Y,Z,DX,DY,DZ)
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      REAL   *8   X, DX
      COMPLEX*16  Y, Z, DY, DZ
C--
      REAL   *8   HALFDX
      COMPLEX*16  D1, D2, D3, D4, E1, E2, E3, E4, FS, GS
C
C===========< Entry Point >=============================================
C
C--
C  Calculate this step.
C--
      HALFDX = DX/2
C--
      D1 = DX*FS( X, Y, Z )
      E1 = DX*GS( X, Y, Z )
      D2 = DX*FS( X+HALFDX, Y+.5D0*D1, Z+.5D0*E1 )
      E2 = DX*GS( X+HALFDX, Y+.5D0*D1, Z+.5D0*E1 )
      D3 = DX*FS( X+HALFDX, Y+.5D0*D2, Z+.5D0*E2 )
      E3 = DX*GS( X+HALFDX, Y+.5D0*D2, Z+.5D0*E2 )
      D4 = DX*FS( X+DX, Y+D3, Z+E3 )
      E4 = DX*GS( X+DX, Y+D3, Z+E3 )
C--
      DY = ( D1 + 2*D2 + 2*D3 + D4 )/6
      DZ = ( E1 + 2*E2 + 2*E3 + E4 )/6
C--
C  That's it.
C--
      RETURN
      END
 
CC*********************************************************************
C*  COMPLEX*16 FUNCTION  FS( X, Y, Z )
CC*********************************************************************
 
      DOUBLE COMPLEX FUNCTION FS(X,Y,Z)
 
      REAL   *8   X
      COMPLEX*16  Y, Z
C
C===========< Entry Point >=============================================
C
C--
C  Calculate FS.
C--
      FS = Z
C--
C  That's it.
C--
      RETURN
      END
 
CC*********************************************************************
C*
C*  COMPLEX*16 FUNCTION  GS(X,Y,Z)
C*
C*  Inputs:
C*       x, y, z,
C*       enerS : E = sqrt(s) - 2m_t ; stored in  common/enerS/...
C*       mt, eff_mt, imag;  stored in common/param/...
C*  Output:
C*       g_S = - eff_mt*( E + i*Gamma0/2 - V(x) ) *y
C*
CC*********************************************************************
 
      DOUBLE COMPLEX FUNCTION GS(X,Y,Z)
 
      IMPLICIT REAL*8 ( A-H, O-Z )
C--
      REAL   *8       X
      COMPLEX*16      Y, Z
C--
      COMMON /TTPARM/ MT, ALFS, ALF, S2W, YCUT, RWT, GAMMAT, ALAMB,
     .                CF, MB, MZ, MW, EFFMT, TLETA1, TLETA2, KAPPA,
     .                PI, IMAG
      REAL   *8       MT, ALFS, ALF, S2W, YCUT, RWT, GAMMAT, ALAMB,
     .                CF, MB, MZ, MW, TLETA1, TLETA2, KAPPA,
     .                PI
      COMPLEX*16      EFFMT, IMAG
C--
      COMMON /ENERS/  ENERS
      REAL   *8       ENERS
C--
      REAL   *8  GAMMA0, V
C
C===========< Entry Point >=============================================
C
C--
C  Calculate GS.
C--
      GS = - EFFMT * ( ENERS + IMAG*GAMMA0(ENERS)/2 - V(X) ) *Y
C--
C  That's it.
C--
      RETURN
      END
