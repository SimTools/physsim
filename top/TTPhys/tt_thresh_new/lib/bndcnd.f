CC*********************************************************************
C*
C*  SUBROUTINE BNDCND( X0, Y1, Z1, Y2, Z2 )
C*
C*  This subroutine set up the boundary condition for
C*  the solutions g0 and g1 to the S-wave Schroedinger equation
C*  at x=x0:
C*
C*       y1(x0) = x0,  y2(x0) = 1
C*
C*       z1(x0) = 1,
C*       z2(x0) = dg1/dx = beta*[loglog(r_QCD/x0)
C*                       + b1/b0^2*(loglog(r_QCD/x0)+1}/log(r_QCD/x0)]
C*  Inputs:
C*       mt, eff_mt, CF, pi; stored in common/param/...
C*       RQCD      ; stored in common/QCDpar/...
C*       b0, b1, beta ; stored in common/optbcd/...
C*  Outputs:
C*       x0 : x where the bondary condition is imposed
C*       y1, y2 : g0(x0) and g1(x0), respctvly
C*       z1, z2 : dg0/dx and dg1/dx at x=x0, respctvly
C*
CC*********************************************************************
 
      SUBROUTINE BNDCND( X0, Y1, Z1, Y2, Z2 )
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      REAL   *8       X0
      COMPLEX*16      Y1, Y2, Z1, Z2
C--
      COMMON /TTPARM/ MT, ALFS, ALF, S2W, YCUT, RWT, GAMMAT, ALAMB,
     .                CF, MB, MZ, MW, EFFMT, TLETA1, TLETA2, KAPPA,
     .                PI, IMAG
      REAL   *8       MT, ALFS, ALF, S2W, YCUT, RWT, GAMMAT, ALAMB,
     .                CF, MB, MZ, MW, TLETA1, TLETA2, KAPPA,
     .                PI
      COMPLEX*16      EFFMT, IMAG
C--
      COMMON /QCDPAR/ ALFST, ALFSP, RQCD
      REAL*8          ALFST, ALFSP, RQCD
C--
      COMMON /OPTBCD/ B0, B1, BETA, OPT
      REAL   *8       B0, B1
      COMPLEX*16      BETA, OPT
C
C===========< Entry Point >=============================================
C
C--
C  Set starting point.
C--
      X0 = .0001D0/MT
C--
C  Set boundary conditions at x0.
C--
      Y1 = X0
      Y2 =  1
      Z1 =  1
      Z2 = BETA*( LOG(LOG(RQCD/X0)) + (B1/B0**2)
     .      *LOG(LOG(RQCD/X0)+1)/LOG(RQCD/X0) )
C--
C  That's it.
C--
      RETURN
      END
