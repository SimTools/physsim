CC**********************************************************************
C*
C*   Real*8 function NTGRD2(m,n)
C*
C*   This function subroutine calculates
C*                                                1
C*      P.V.{Re[ conj{tilG(q;E)}*tilG(p;E) ] * -------}
C*                                              x - 1
C*   PS(M) : integration variable = q
C*   PS(N) : top momentum         = p
C*
CC**********************************************************************
 
      DOUBLE PRECISION FUNCTION  NTGRD2(M,N)
 
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
      IF ( M .EQ. N )  THEN
            PRREMN = REAL( CONJG(TILG(M-1))*TILG(N) )
     .               - REAL( CONJG(TILG(N))*TILG(N) )
            X      = PS(M-1)/PS(N)
            BEFORE = PRREMN/( X - 1 )
            PRREMN = REAL( CONJG(TILG(M+1))*TILG(N) )
     .               - REAL( CONJG(TILG(N))*TILG(N) )
            X      = PS(M+1)/PS(N)
            AFTER  = PRREMN/( X - 1 )
            NTGRD2 = ( BEFORE + AFTER )/2
      ELSE
            PRREMN = REAL( CONJG(TILG(M))*TILG(N) )
     .               - REAL( CONJG(TILG(N))*TILG(N) )
            X      = PS(M)/PS(N)
            NTGRD2 = PRREMN/( X - 1 )
      END IF
C--
C  That's it.
C--
      RETURN
      END
