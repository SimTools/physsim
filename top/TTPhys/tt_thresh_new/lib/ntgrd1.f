CC**********************************************************************
C*
C*  Real*8 function NTGRD1(m,n)
C*
C*  This function subroutine calculates
C*
C*          Re[ conj{tilG(q;E)}*tilG(p;E) ] * w1(q/p)
C*
C*  with
C*                   2x + 1          ! x + 1 !
C*          w1(x) = -------- - x*log !-------!  .
C*                   x  + 1          ! x - 1 !
C*
C*      PS(M) : integration variable = q
C*      PS(N) : top momentum         = p
C*
CC**********************************************************************
 
      DOUBLE PRECISION FUNCTION  NTGRD1(M,N)
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      INTEGER*4       M, N
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
C  Statement function.
C--
      W1(X) = ( 2*X + 1 )/( X + 1 ) - X*LOG(ABS( (X+1)/(X-1) ))
C
C===========< Entry Point >=============================================
C
C--
C  Calculate integrand.
C--
      REMN = REAL( CONJG(TILG(M))*TILG(N) )
      IF ( M .EQ. N ) THEN
            EPS = ALAMB/100/PS(N)
            Y1  = REAL( CONJG(TILG(M-1))*TILG(N) )
            Y2  = REMN
            Y3  = REAL( CONJG(TILG(M+1))*TILG(N) )
            SUM = Y3 + Y1
            DIF = Y3 - Y1
            Z2  = 1/32.D0*(
     .             ( EPS*(SUM-6*Y2) - (2+8/EPS**2)*DIF
     .               + 16/EPS**3*(SUM-2*Y2) )
     .                         *LOG( (2-EPS)/(2+EPS) )
     .             + 16*Y2*LOG( (4-EPS**2)/EPS**2 )
     .             + 4*(SUM+10*Y2) - 8/EPS*DIF
     .                           + 16/EPS**2*(SUM-2*Y2) )
            NTGRD1 = REMN * 1.5D0 - Z2
      ELSE
            NTGRD1 = REMN * W1( PS(M)/PS(N) )
      END IF
C--
C  That's it.
C--
      RETURN
      END
