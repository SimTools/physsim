CC*********************************************************************
C*
C*  SUBROUTINE QCDS1( E, DR, EPS, STEP1 )
C*
C*  This subroutine solves the S-wave Schroedinger equation:
C*
C*        dy/dx = z
C*        dz/dx = -eff_mt*(E+iGamma0/2-V(x))*y
C*
C*  using the Runge-Kutta method.
C*
C*  It determines the two solutions g0(x) and g1(x), which satisfy
C*  the following boundary conditions at the origin
C*
C*       g0(x)  = x +... ,
C*       dg0/dx = 1 +... ,
C*
C*       g1(x)  = 1 +... ,
C*       dg1/dx = beta*[ loglog(r_QCD/x)
C*                  + b1/b0^2*{loglog(r_QCD/x)+1}/log(r_QCD/x) + ...]
C*
C*  g0(x) and g1(x) are given in common/SFIRST/...
C*
C*  Input:
C*       E     : energy = sqrt(s) - 2*m_t  ; it will be stored in
C*                       common/enerS/...  to be used in the
C*                               function subroutine g_S
C*       dr    : integration step
C*       eps   : criterion for the convergence of B;
C*               integration continued while  !B_old/b_new-1! > eps
C*
C*  Output: xx(nn), g0(nn), g1(nn), and nnmax are stored in /SFIRST/.
C*
C*       step1  : number of steps required to obtain convergent
C*                B = - lim_{x->infty} y1/y2
C*       xx     : x
C*       g0, g1 : solutions satisfying the above boundary conditions
C*       nnmax  : number of points for g0 and g1
C*
CC*********************************************************************
 
      SUBROUTINE QCDS1( E, DR, EPS, STEP1 )
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      REAL*8          E, DR, EPS
      INTEGER*4       STEP1
C--
      COMMON /SFIRST/ XX, G0, G1, NNMAX
      REAL   *8       XX(5000)
      COMPLEX*16      G0(5000), G1(5000)
      INTEGER*4       NNMAX
C--
      COMMON /ENERS/  ENERS
      REAL   *8       ENERS
C--
      REAL   *8       X0, X, DX, CRITER
      INTEGER*4       MXSTPS, N, CLOCK, K, NN
      COMPLEX*16      Y1, Y2, Z1, Z2,
     .                DY1, DZ1, DY2, DZ2,
     .                B, BNEW
      PARAMETER ( MXSTPS = 40000 )
C
C===========< Entry Point >=============================================
C
C--
C  Initialization.
C--
      WRITE(*,*)
      WRITE(*,*) 'S-WAVE: E=', E
      ENERS = E
C--
C  Set up boundary conditions at x = x0.
C--
      CALL BNDCND( X0, Y1, Z1, Y2, Z2 )
      X     = X0
      XX(1) = X
      G0(1) = Y1
      G1(1) = Y2
C--
C  Do loop while  ( CRITER > EPS ).
C--
      DX     = DR
      N      = 0
      B      = 0
      NN     = 1
      CRITER = 100
C--
1000  CONTINUE
         DO 100 K = 1, 2
            DO 10 CLOCK = 1, 30
               CALL  RKUTAS(X,Y1,Z1,DX,DY1,DZ1)
               CALL  RKUTAS(X,Y2,Z2,DX,DY2,DZ2)
               Y1 = Y1 + DY1
               Z1 = Z1 + DZ1
               Y2 = Y2 + DY2
               Z2 = Z2 + DZ2
               X  = X + DX
               N  = N + 1
10          CONTINUE
            NN     = NN + 1
            XX(NN) = X
            G0(NN) = Y1
            G1(NN) = Y2
100      CONTINUE
         BNEW = - Y2/Y1
         CRITER = ABS(B/BNEW-1)
         B = BNEW
         IF ( N .GT. MXSTPS ) THEN
            WRITE(*,*)  'integration steps exceeded',MXSTPS
                                                 GO TO 9999
         END IF
      IF ( CRITER .GT. EPS .OR. N .LT. 500 )     GO TO 1000
C--
C  Outputs.
C--
      NNMAX = NN
      STEP1 = NN
C>>>
C     WRITE(*,*) 'STEP1, B, BNEW, CRITER'
C     WRITE(*,*)  STEP1
C     WRITE(*,*)  B, BNEW
C     WRITE(*,*)  CRITER
C>>>
C--
C  That's it.
C--
9999  RETURN
      END
