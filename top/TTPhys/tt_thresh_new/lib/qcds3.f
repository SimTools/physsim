CC*********************************************************************
C*
C*  SUBROUTINE QCDS3
C*
C*  This subroutine calculates the Fourier Transform
C*  of G(x;E) :
C*
C*                                   _infty
C*                         eff_mt   /
C*         tilde{G}(p;E) = ------- /  dr  gg(r)*sin(pr).
C*                            p  _/ 0
C*
C*  tilde{G}(p;E) is given at 101 points for 0 =< p =< Lambda
C*
C*  Inputs:
C*       dst_S(num), gg(num) : function gg(r) ;
C*       point_max : number of points of gg(r) ;
C*                               stored in common/SSECND/...
C*       eff_mt    : complex effective mass ; stored in common/param/.
C*       Lambda : maximum momentum ;
C*                               stored in common/param/...
C*  Outputs:
C*       p_S(k), tilg(k)  : Green's function in the momentum space
C*                       = tilde{G}(p;E)  ;
C*                               stored in common/STHIRD/...
C*
C*  Remark:  We set p_S(0)=0, tilg(0)=0.
C*
CC*********************************************************************
 
      SUBROUTINE QCDS3
 
      IMPLICIT REAL*8 ( A-H, O-Z )
C--
      COMMON /SSECND/ DSTS, GG, PNTMAX
      REAL   *8       DSTS(5000)
      COMPLEX*16      GG(5000)
      INTEGER*4       PNTMAX
C--
      COMMON /TTPARM/ MT, ALFS, ALF, S2W, YCUT, RWT, GAMMAT, ALAMB,
     .                CF, MB, MZ, MW, EFFMT, TLETA1, TLETA2, KAPPA,
     .                PI, IMAG
      REAL   *8       MT, ALFS, ALF, S2W, YCUT, RWT, GAMMAT, ALAMB,
     .                CF, MB, MZ, MW, TLETA1, TLETA2, KAPPA,
     .                PI
      COMPLEX*16      EFFMT, IMAG
C--
      COMMON /STHIRD/ PS, TILG
      REAL   *8       PS(0:100)
      COMPLEX*16      TILG(0:100)
C--
      INTEGER*4       N, K
      REAL   *8       Q, R1, R2, R3, CRITER
      COMPLEX*16      GL1, GL2, GL3, NTGRND, I
 
 
      PARAMETER ( CRITER = 1.D-3 )
 
C
C===========< Entry Point >=============================================
C
C--
C  For p = 0.
C--
      PS(0)   = 0
      TILG(0) = 0
C--
C  Do loop for  0 < p =< ALAMB in 100 steps.
C--
      DO 1000 N = 1, 100
         Q      = ALAMB*REAL(N)/100
         NTGRND = 0
         R1     = DSTS(1)
         GL1    = GG(1)
         DO 100 K = 2, PNTMAX-1, 2
            R2 = DSTS(K)
            GL2 = GG(K)
            R3 = DSTS(K+1)
            GL3 = GG(K+1)
            IF ( ((R3-R2)/(R2-R1)-1) .GT. CRITER )
     .         WRITE(*,*)  'Integration steps not equal'
            NTGRND = NTGRND
     .          + ( GL1*SIN(Q*R1) + 4*GL2*SIN(Q*R2) + GL3*SIN(Q*R3) )
     .            * (R3-R1)
            R1     = R3
            GL1    = GL3
100       CONTINUE
          I = EFFMT/Q /6*NTGRND
C>>>
C         WRITE (40,*)  q, q**2*abs(I)**2
C         WRITE (50,*)  q, real(I), aimag(I), abs(q*I)**2*.002
C>>>
          PS(N)   = Q
          TILG(N) = I
1000  CONTINUE
C>>>
C      WRITE (40,*)  'JOIN'
C      WRITE (50,*)  'JOIN'
C>>>
C--
C  That's it.
C--
      RETURN
      END
