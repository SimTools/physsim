CC*********************************************************************
C*
C*  SUBROUTINE QCDP2( STEP2 )
C*
C*  This subroutine calculates the factor B_P from the
C*  asymptotic behavior of y1 and y2.
C*  Then determine f(x) which is the solution to the
C*  homogeneous P-wave Schroedinger equation with the
C*  boundary condition
C*
C*       f(0) = 1,    f(x) -> 0  as x->infty .
C*
C*  This solution is given by  f(x) = A [ y1(x) + B_P*y2(x) ].
C*
C*  Inputs:
C*       x_L, y1_L, y2_L : two independent solutions at x>=0.01
C*                         ; stored in  common/PFRSTL/...
C*       x_S, y1_S, y2_S : two independent solutions at x<0.01
C*                         ; stored in  common/PFRSTS/...
C*       num_L, num_S    : number of points for the solutions,
C*                         for x>=0.01 and x<0.01, respctvly.
C*  Outputs:
C*       dst_P(num), f(num)  : solution satisfying the above boundary
C*                         conditions ; stored in  common/PSECND/...
C*       num_max         : number of points for f(x) ;
C*                               also stored in common/PSECND/...
C*       step2 : number of steps required for sufficiently convergent
C*               f(x) for x>=0.01
C*
CC*********************************************************************
 
      SUBROUTINE  QCDP2( STEP2 )
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      INTEGER*4       STEP2
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
      COMMON /PSECND/ DSTP, F, NUMMAX
      REAL   *8       DSTP(5000)
      COMPLEX*16      F(5000)
      INTEGER*4       NUMMAX
C--
      REAL   *8      NORM, NORMAX, CRITER
      COMPLEX*16     A, BP, BTMP, T
      INTEGER*4      NUM, K
C--
C  Set criterion for truncating the Green's function f(x).
C--
      PARAMETER ( CRITER = 1.D-5 )
C
C===========< Entry Point >=============================================
C
C--
C  Determination of BP = - lim_{x->infty} y1(x)/y2(x).
C--
      T = 0
      DO 100 NUM = NUML-4, NUML
         BTMP = - Y1L(NUM)/Y2L(NUM)
         T    = T + BTMP
100   CONTINUE
      BP = T/5
C--
C  Determination of A = lim_{x->0} 1/{ X*[ y1 + BP*y2 ] }.
C--
      A = 1/XS(NUMS)/( Y1S(NUMS) + BP*Y2S(NUMS) )
C--
C  Determination of f(x) for x < 0.01.
C--
      NUM = 0
      DO 200 K = NUMS, 1, -1
         NUM       = NUM + 1
         DSTP(NUM) = XS(K)
         F(NUM)    = A*( Y1S(K) + BP*Y2S(K) )
C>>>
C        WRITE(20,*) DSTP(NUM), DSTP(NUM)**2*ABS(F(NUM))**2
C>>>
200   CONTINUE
C--
C  Determination of f(x) for x >= 0.01.
C--
      K      = 0
      NORMAX = 0
3000  CONTINUE
         K = K + 1
         NUM = NUM + 1
         DSTP(NUM) = XL(K)
         F(NUM) = A*( Y1L(K) + BP*Y2L(K) )
         NORM = DSTP(NUM)**2*ABS(F(NUM))**2
         NORMAX = MAX( NORM, NORMAX )
C>>>
C        WRITE(20,*) DSTP(NUM), NORM
C>>>
         IF ( NORM .LT. CRITER*NORMAX )          GO TO 400
      IF ( K .LT. NUML )                         GO TO 3000
      WRITE(*,*) 'f(r) does not become sufficiently small.'
C--
C  Loop end.
C--
400   STEP2  = K
      NUMMAX = INT((NUM-1)/2)*2 + 1
C>>>
C     WRITE(*,*) 'total number of points, used number of points'
C     WRITE(*,*) NUMS+NUML, NUMMAX
C>>>
C--
C  That's it.
C--
      RETURN
      END
