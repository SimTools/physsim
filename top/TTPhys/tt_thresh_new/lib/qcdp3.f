CC*********************************************************************
C*
C*  SUBROUTINE QCDP3
C*
C*  This subroutine calculates the Fourier Transform
C*  of F^l(x;E) :
C*
C*                                  _infty
C*                         eff_mt  /
C*         tilde{F}(p;E) = ------ /  dr  f(r)*delfn(pr),
C*                          p^2 _/ 0
C*
C*  where the Fourier kernel delfn(x) is defined as
C*
C*         delfn(x) = (sin(x)-x*cos(x))/x .
C*
C*  tilde{F}(p;E) is given at 101 points for 0 =< p =< Lambda
C*
C*  Inputs:
C*       dst_P(num), f(num) : function f(r) ;
C*       num_max : number of points of f(r) ;
C*                               stored in common/PSECND/...
C*       eff_mt  : top quark mass ;
C*       Lambda  : maximum momentum ;
C*                               stored in common/TTPARM/..
C*  Outputs:
C*       p_P(k), tilf(k)  : Green's function in the momentum space
C*                       = tilde{F}(p;E)  ;
C*                               stored in common/PTHIRD/...
C*
C*  Remark:  We set p_P(0)=0, tilf(0)=0 .
C*
CC*********************************************************************
 
      SUBROUTINE QCDP3
 
      IMPLICIT REAL*8 ( A-H, O-Z )
C--
      COMMON /PSECND/ DSTP, F, NUMMAX
      REAL   *8       DSTP(5000)
      COMPLEX*16      F(5000)
      INTEGER*4       NUMMAX
C--
      COMMON /TTPARM/ MT, ALFS, ALF, S2W, YCUT, RWT, GAMMAT, ALAMB,
     .                CF, MB, MZ, MW, EFFMT, TLETA1, TLETA2, KAPPA,
     .                PI, IMAG
      REAL   *8       MT, ALFS, ALF, S2W, YCUT, RWT, GAMMAT, ALAMB,
     .                CF, MB, MZ, MW, TLETA1, TLETA2, KAPPA,
     .                PI
      COMPLEX*16      EFFMT, IMAG
C--
      COMMON /PTHIRD/ PP, TILF
      REAL   *8       PP(0:100)
      COMPLEX*16      TILF(0:100)
C--
      INTEGER*4       N, K
      REAL   *8       Q, R1, R2, R3, DELFN, CRITER
      COMPLEX*16      F1, F2, F3, NTGRND
C--
      PARAMETER ( CRITER = 1.D-3 )
C
C===========< Entry Point >=============================================
C
C--
C  For p = 0.
C--
C>>>
C      WRITE (40,*)  0., 0.
C>>>
      PP(0)   = 0
      TILF(0) = 0
C--
C  Do loop for  0 < p =< Lambda in 100 steps.
C--
      DO 1000 N=1, 100
          Q      = ALAMB*REAL(N)/100
          NTGRND = 0.
          R1     = DSTP(1)
          F1     = F(1)
          DO 100 K = 2, NUMMAX-1, 2
             R2 = DSTP(K)
             F2 = F(K)
             R3 = DSTP(K+1)
             F3 = F(K+1)
             IF ( ((R3-R2)/(R2-R1)-1) .GT. CRITER )
     .          WRITE(*,*)  'INTEGRATION STEPS NOT EQUAL'
             NTGRND = NTGRND
     .                + ( F1*DELFN(Q*R1) + 4*F2*DELFN(Q*R2)
     .                    + F3*DELFN(Q*R3) )* (R3-R1)
             R1 = R3
             F1 = F3
100       CONTINUE
C--
          PP(N)   = Q
          TILF(N) = EFFMT/Q**2 /6*NTGRND
C>>>
C          WRITE (40,*)  Q, 4*PI*Q**2*ABS(TILF(N))**2 *Q
C>>>
1000  CONTINUE
C--
C  That's it.
C--
      RETURN
      END
