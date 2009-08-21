C***********************************************************************
C* In this subroutine the inputs are x,y,z,dx and the outputs are
C* dy and dz.
C***********************************************************************
 
      SUBROUTINE RKUTTA(X,Y,Z,DX,DY,DZ)
 
      REAL   *8  X, DX
      COMPLEX*16 Y, Z, DY, DZ
C--
      COMMON /PARAM1/ E, GAMMA0, ALF, ALFS, S2W, AMT, AMW, AMZ, GMZ,
     .                AMB, VTB, AMH, BTH, ALFST
      REAL   *8       E, GAMMA0, ALF, ALFS, S2W, AMT, AMW, AMZ, GMZ,
     .                AMB, VTB, AMH, BTH, ALFST
      COMMON /EFFMAS/ ZMU, BETA
      COMPLEX*16      ZMU, BETA
C--
      REAL   *8       V
C--
      REAL   *8  HALFDX
      COMPLEX*16 D1, D2, D3, D4, E1, E2, E3, E4, G
C--
      COMPLEX*16      ZI
      PARAMETER ( ZI = (0.D0,1.D0) )
C--
      G(X,Y) = - ZMU*( E + .5D0*ZI*GAMMA0 - V(X) )*Y
C
C========< Entry Point >================================================
C
      HALFDX = .5D0*DX
 
      D1 = DX*Z
      E1 = DX*G( X, Y )
 
      D2 = DX*( Z+.5D0*E1 )
      E2 = DX*G( X+HALFDX, Y+.5D0*D1 )
 
      D3 = DX*( Z+.5D0*E2 )
      E3 = DX*G( X+HALFDX, Y+.5D0*D2 )
 
      D4 = DX*( Z+E3 )
      E4 = DX*G( X+DX, Y+D3 )
 
      DY = ( D1 + 2*D2 + 2*D3 + D4 )/6
      DZ = ( E1 + 2*E2 + 2*E3 + E4 )/6
 
      RETURN
      END
