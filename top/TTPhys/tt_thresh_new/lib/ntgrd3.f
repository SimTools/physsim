CC**********************************************************************
C*
C*   Real*8 functin  NTGRD3(m,n)
C*
C*   This function subroutine calculates
C*                                                1
C*           Re[ conj{tilG(q;E)}*tilG(p;E) ] * -------
C*                                              x - 1
C*   PS(m) : integration variable = q
C*   PS(n) : top quark momentum   = p
C*
CC**********************************************************************
 
      DOUBLE PRECISION FUNCTION  NTGRD3(M,N)
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      INTEGER*4       M, N
C--
      COMMON /STHIRD/ PS, TILG
      REAL   *8       PS(0:100)
      COMPLEX*16      TILG(0:100)
C
C===========< Entry Point >=============================================
C
C--
C  Calculate integrand.
C--
      X      = PS(M)/PS(N)
      REMN   = REAL( CONJG(TILG(M))*TILG(N) )
      NTGRD3 = REMN/( X - 1 )
C--
C  That's it.
C--
      RETURN
      END
