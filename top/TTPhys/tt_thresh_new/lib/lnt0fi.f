CC**********************************************************************
C*
C*  Real*8 function  LNT0FI(n)
C*
C*  This function subroutine calculates the loop integral
C*
C*                   _ Lambda
C*                  /
C*        LNT(n) = /  dq Im[ conj{tilG(q;E)}*tilG(p;E) ] * w0(q/p)
C*               _/ 0
C*
C*  with weight function
C*
C*                   x         x
C*        w0(x) = ------- - -------
C*                !x - 1!    x + 1
C*
CC**********************************************************************
 
      DOUBLE PRECISION FUNCTION  LNT0FI(N)
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      INTEGER*4       N
C--
      COMMON /STHIRD/ PS, TILG
      REAL   *8       PS(0:100)
      COMPLEX*16      TILG(0:100)
C--
      COMMON /TTPARM/ MT, ALFS, ALF, S2W, YCUT, RWT, GAMMAT, ALAMB,
     .                CF, MB, MZ, MW, EFFMT, TLETA1, TLETA2, KAPPA,
     .                PI, IMAG
      REAL   *8       MT, ALFS, ALF, S2W, YCUT, RWT, GAMMAT, ALAMB,
     .                CF, MB, MZ, MW, TLETA1, TLETA2, KAPPA,
     .                PI
      COMPLEX*16      EFFMT, IMAG
C--
      REAL   *8       SUM, DQ, SIMP1, SIMP2, SIMP3, W0
      REAL   *8       INT1, INT2
C--
      COMMON /KFLAGS/ KFNOFI 
      INTEGER*4       KFNOFI
C--
C  Statement function.
C--
      W0(X) = X/ABS( X - 1 ) - X/( X + 1 )
C
C===========< Entry Point >=============================================
C
C--
C  Check if FI corrections requested.
C--
      IF ( KFNOFI.EQ.1 ) THEN
         LNT0FI = 0
         RETURN
      ENDIF
C--
C  Initialization.
C--
      DQ    = ALAMB/50
      SUM   = 0
      SIMP1 = 0
C--
C  Numerical integration starts here.
C--
C>>>
C     DO 100  M = 1, 99, 2
      DO 100  M = 1, 98, 2
C>>>
          IF ( M .EQ. N ) THEN
               INT1  = W0( PS(M-1)/PS(N) )
     .                 *DIMAG( CONJG(TILG(M-1))*TILG(N) )
               INT2  = W0( PS(M+1)/PS(N) )
     .                 *DIMAG( CONJG(TILG(M+1))*TILG(N) )
               SIMP2 = ( INT1 + INT2 )/2
          ELSE
               SIMP2 = W0( PS(M)/PS(N) )
     .                 *DIMAG( CONJG(TILG(M))*TILG(N) )
          END IF
C--
          IF ( M+1 .EQ. N ) THEN
               INT1  = W0( PS(M)/PS(N) )
     .                 *DIMAG( CONJG(TILG(M))*TILG(N) )
               INT2  = W0( PS(M+2)/PS(N) )
     .                 *DIMAG( CONJG(TILG(M+2))*TILG(N) )
               SIMP3 = ( INT1 + INT2 )/2
          ELSE
               SIMP3 = W0( PS(M+1)/PS(N) )
     .                 *DIMAG( CONJG(TILG(M+1))*TILG(N) )
          END IF
          SUM   = SUM + ( SIMP1 + 4*SIMP2 + SIMP3 )
          SIMP1 = SIMP3
100   CONTINUE
      LNT0FI = SUM*DQ/6
C--
C  That's it.
C--
      RETURN
      END
