CC*********************************************************************
C*
C*  Subroutine QCDS2( STEP2 )
C*
C*  This subroutine calculates the factor B from the
C*  asymptotic behavior of g0 and g1.
C*  Then determine gg(x) which is the solution to the
C*  homogeneous S-wave Schroedinger equation with the
C*  boundary condition
C*
C*       gg(0) = 1,    gg(x) -> 0  as x->infty .
C*
C*  This solution is given by  gg(x) = g1(x) + B*g0(x) .
C*
C*  Inputs:      stored in common/SFIRST/...
C*       xx, g0, g1 : two solutions determined in subroutine QCDS1 .
C*       NNMAX      : number of points for the solutions g0 and g1
C*  Outputs:
C*       dst_S(nn), gg(nn) : solution satisfying the above boundary
C*                           conditions ; stored in  common/SSECND/.
C*       point_max : number of points for gg(x) ;
C*                               also stored in common/SSECND/...
C*       step2 : number of steps required for sufficiently convergent
C*               gg(x)
C*       B     : constant related to Im G(0;E) ; stored in common/ImG/
C*
CC*********************************************************************
 
      SUBROUTINE  QCDS2( STEP2 )
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      INTEGER*4       STEP2
C--
      COMMON /SFIRST/ XX, G0, G1, NNMAX
      REAL   *8       XX(5000)
      COMPLEX*16      G0(5000), G1(5000)
      INTEGER*4       NNMAX
C--
      COMMON /SSECND/ DSTS, GG, PNTMAX
      REAL   *8       DSTS(5000)
      COMPLEX*16      GG(5000)
      INTEGER*4       PNTMAX
C--
      COMMON /IMG/    B
      COMPLEX*16      B
C--
      REAL   *8       NORM, NORMAX, CRITER
      COMPLEX*16      BTMP, T
C--
C  Set criterion for truncating the Green's function GG(X).
C--
      PARAMETER ( CRITER = 1.D-7 )
C
C===========< Entry Point >=============================================
C
C--
C  Determination of B = - lim_{x->infty} G1/G0.
C--
      T = 0
      DO 100 NN = NNMAX-4, NNMAX
         BTMP = - G1(NN)/G0(NN)
         T = T + BTMP
100   CONTINUE
      B = T/5
C--
C  Loop until GG converges.
C--
      NN     = 0
      NORMAX = 0
C--
200   CONTINUE
         NN       = NN + 1
         GG(NN)   = G1(NN) + B*G0(NN)
         DSTS(NN) = XX(NN)
         NORM     = ABS(GG(NN))**2
         NORMAX   = MAX( NORM, NORMAX )
C>>>
C        WRITE (35,*) DSTS(NN), REAL(GG(NN)), AIMAG(GG(NN))
C        WRITE (36,*) DSTS(NN), NORM
C>>>
          IF ( NORM .LT. CRITER*NORMAX )         GO TO 300
      IF ( NN .LT. NNMAX )                       GO TO 200
C--
C  No convergence.
C--
      WRITE(*,*) 'g>(r) does not become sufficiently small.'
C--
C  Outputs.
C--
 300  STEP2  = NN
      PNTMAX = INT((NN-1)/2)*2 + 1
C>>>
C     WRITE(*,*) 'STEP1, STEP2'
C     WRITE(*,*)  NNMAX, STEP2
C>>>
C--
C  That's it.
C--
      RETURN
      END
