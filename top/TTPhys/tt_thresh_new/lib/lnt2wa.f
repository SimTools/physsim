CC**********************************************************************
C*
C*   Real*8 function  LNT2WA(n)
C*
C*   This function subroutine calculates the loop integration
C*                     _ Lambda
C*                    /
C*         LNT(n) =  /  dq  q * !tilG(q;E)!^2 * Phs2(q,p)
C*                 _/ 0
C*
C*
CC**********************************************************************
 
      DOUBLE PRECISION FUNCTION  LNT2WA(N)
 
      IMPLICIT REAL*8 ( A-H, O-Z )
      INTEGER*4     N
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
 
      REAL   *8       PHS2, INT1, INT2
C
C===========< Entry Point >=============================================
C
C--
C  Initialization.
C--
      DQ    = ALAMB/50
      SUM   = 0
      SIMP1 = 0
C--
C  Start numerical integration.
C--
      DO 100 M = 1, 99, 2
          IF ( M .EQ. N ) THEN
             INT1  = PS(M)*ABS(TILG(M-1))**2 * PHS2( PS(M-1), PS(N) )
             INT2  = PS(M)*ABS(TILG(M+1))**2 * PHS2( PS(M+1), PS(N) )
             SIMP2 = ( INT1 + INT2 )/2
          ELSE
             SIMP2 = PS(M)*ABS(TILG(M))**2   * PHS2( PS(M), PS(N) )
          END IF
          IF ( M+1 .EQ. N ) THEN
             INT1  = PS(M)*ABS(TILG(M))**2   * PHS2( PS(M), PS(N) )
             INT2  = PS(M)*ABS(TILG(M+2))**2 * PHS2( PS(M+2), PS(N) )
             SIMP3 = ( INT1 + INT2 )/2
          ELSE
             SIMP3 = PS(M)*ABS(TILG(M+1))**2 * PHS2( PS(M+1), PS(N) )
          END IF
          SUM   = SUM + ( SIMP1 + 4*SIMP2 + SIMP3 )
          SIMP1 = SIMP3
100   CONTINUE
      LNT2WA = SUM*DQ/6
C--
C  That's it.
C--
      RETURN
      END
