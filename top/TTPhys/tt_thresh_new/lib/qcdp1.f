CC*********************************************************************
C*
C*  SUBROUTINE QCDP1( E, DR, EPS, STEP1 )
C*
C*  This subroutine solves the P-wave Schroedinger equation:
C*
C*        dy/dx = z
C*        dz/dx = ( -eff_mt*(E+iGamma0/2-V(x)) + 2/x^2 )*y
C*
C*  using the Lunge-Kutta method.
C*
C*  It determines the two independent solutions
C*  y1(x) and y2(x), which are determined by the conditions at
C*  x=0.01 GeV^-1 :
C*
C*        y1 = ( 4.d-3, -7.d-4 ), z1 = ( 4.d-2, -2.d-2 )
C*        y2 = ( 2. d8, -3. d7 ), z2 = ( 2. d9, -1. d9 ) .
C*
C*  y1(x) and y2(x) are given in common/PFRSTL/... for x>0.01,
C*  and in common/PFRSTS/... for x<0.01 .
C*
C*  Input:
C*       E     : energy = sqrt(s) - 2*m_t  ;  it will be stored in
C*                                common/ener/... to be used in the
C*                                function subroutine  g_P .
C*       dr    : integration step used in the region x>0.01
C*       eps   : criterion for the convergence of B;
C*               integration continued while  !B_old/B_new-1! > eps
C*
C*  Output: x_L(num), y1_L(num), y2_L(num), and num_L
C*                               are stored in common/PFRSTL/...
C*          x_S(num), y1_S(num), y2_S(num), and num_S
C*                               are stored in common/PFRSTS/...
C*
C*       step1 : number of steps required in x >= 0.01 to obtain
C*               convergent  B = - lim_{x->infty} y1/y2
C*       x_L, x_S   : x < 0.01 and x > 0.01, respectively
C*       y1_L, y1_S : solution y1 corresponding to x_L and x_S, rsptvly
C*       y2_L, y2_S : solution y2 corresponding to x_L and x_S, rsptvly
C*       num_L : number of points for x >= 0.01
C*       num_S : number of points for x < 0.01
C*
CC*********************************************************************
 
      SUBROUTINE QCDP1( E, DR, EPS, STEP1 )
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      REAL   *8       E, DR, EPS
      INTEGER*4       STEP1
C--
      COMMON /PFRSTL/ XL, Y1L, Y2L, NUML
      REAL   *8       XL(5000)
      COMPLEX*16      Y1L(5000), Y2L(5000)
      INTEGER*4        NUML
C--
      COMMON /PFRSTS/ XS, Y1S, Y2S, NUMS
      REAL   *8       XS(5000)
      COMPLEX*16      Y1S(5000), Y2S(5000)
      INTEGER*4       NUMS
C--
      COMMON /ENER/   ENER
      REAL   *8       ENER
C--
      REAL*8          X, DX, CRITER
      INTEGER*4       MXSTPS, N, CLOCK, K, NUM
      COMPLEX*16      Y1, Y2, Z1, Z2,
     .                DY1, DZ1, DY2, DZ2,
     .                B, BNEW, R1, R2, R1NEW, R2NEW
C--
C  Maximum number of steps allowed.
C--
      PARAMETER ( MXSTPS = 40000 )
C
C===========< Entry Point >=============================================
C
C--
C  Initialization.
C--
      WRITE(*,*) 'P-WAVE: E=', E
      ENER = E
C--
C  FOR r >= 0.01.
C    Conditions at x = 0.01.
C--
      X       = 0.01D0
      Y1      = (4.D-3,-7.D-4)
      Z1      = (4.D-2,-2.D-2)
      Y2      = (2.D8,-3.D7)
      Z2      = (2.D9,-1.D9)
      XL(1)   = X
      Y1L(1)  = Y1
      Y2L(1)  = Y2
C--
C  Do loop for x > 0.01  continued while  ( CRITER > EPS ).
C--
      N      = 0
      B      = 0
      CRITER = 100
      DX     = DR
      NUM    = 1
1000  CONTINUE
         DO 100 K = 1, 2
            DO 10 CLOCK = 1, 20
               CALL  RKUTAP(X,Y1,Z1,DX,DY1,DZ1)
               CALL  RKUTAP(X,Y2,Z2,DX,DY2,DZ2)
               Y1 = Y1 + DY1
               Z1 = Z1 + DZ1
               Y2 = Y2 + DY2
               Z2 = Z2 + DZ2
               X  = X + DX
               N  = N + 1
10          CONTINUE
            NUM       = NUM + 1
            XL(NUM)   = X
            Y1L(NUM) = Y1
            Y2L(NUM) = Y2
100      CONTINUE
C--
C  Calculate criterion.
C--
         BNEW = - Y1/Y2
         CRITER = ABS(B/BNEW-1)
         B = BNEW
         IF ( N .GT. MXSTPS ) THEN
            WRITE(*,*)  'Integration steps exceeded',MXSTPS
            WRITE(*,*)  CRITER
         END IF
      IF ( CRITER .GT. EPS .OR. N .LT. 500 )     GO TO 1000
C--
C  End of Do.
C--
      NUML = NUM
      STEP1 = NUM
C>>>
C     WRITE(*,*) 'r>=0.01 completed: step1, B, B_new, criterion'
C     WRITE(*,*)  STEP1
C     WRITE(*,*)  B, BNEW
C     WRITE(*,*)  CRITER
C>>>
C--
C  For r < 0.01.
C    Conditions at x = 0.01.
C--
      X  = 0.01D0
      Y1 = (4.D-3,-7.D-4)
      Z1 = (4.D-2,-2.D-2)
      Y2 = (2.D8,-3.D7)
      Z2 = (2.D9,-1.D9)
C--
C  Do loop for x < 0.01  continued while  ( CRITER > 1.d-3 ).
C--
      N      = 0
      R1     = 0
      R2     = 0
      CRITER = 100
      NUM    = 0
2000  CONTINUE
         DX = -X*.02D0
         DO 200 K = 1, 2
            DO 20 CLOCK = 1, 10
               CALL  RKUTAP(X,Y1,Z1,DX,DY1,DZ1)
               CALL  RKUTAP(X,Y2,Z2,DX,DY2,DZ2)
               Y1 = Y1 + DY1
               Z1 = Z1 + DZ1
               Y2 = Y2 + DY2
               Z2 = Z2 + DZ2
               X  = X + DX
               N  = N + 1
20          CONTINUE
            NUM      = NUM + 1
            XS(NUM)  = X
            Y1S(NUM) = Y1
            Y2S(NUM) = Y2
200      CONTINUE
         R1NEW  = X*Y1
         R2NEW  = X*Y2
         CRITER = ABS(R1/R1NEW-1.)+ABS(R2/R2NEW-1.)
         R1     = R1NEW
         R2     = R2NEW
         IF ( N .GT. MXSTPS ) THEN
            WRITE(*,*)  'integration steps exceeded',MXSTPS
            WRITE(*,*)  CRITER
                                                 GO TO 9999
         END IF
      IF ( CRITER .GT. 1.D-3 )                   GO TO 2000
C--
C  End of Do.
C--
      NUMS = NUM
C>>>
C     WRITE(*,*) 'r<0.01 completed: n, r1, r2, criterion'
C     WRITE(*,*) N, R1, R2
C     WRITE(*,*) CRITER
C     WRITE(*,*) 'X =', X
C     WRITE(*,*) 'STEPS =', NUM
C>>>
C--
C  That's it.
C--
9999  RETURN
      END
