CC**********************************************************************
C*
C*  Real*8 function LNT1FI(n)
C*
C*  This function subroutine calculates the loop integral
C*
C*                   _ Lambda
C*                  /
C*        LNT(n) = /  dq Re[ conj{tilG(q;E)}*tilG(p;E) ] * w(q/p)
C*               _/ 0
C*
C*  with weight function
C*
C*                     1       2x + 1          ! x + 1 !
C*        w(x) = Pr.------- + -------- - x*log !-------!
C*                   x - 1     x  + 1          ! x - 1 !
C*
C*
C*  In order to take the principal value, integral region is
C*  classified into the singular region ( n-2<=m<=n+2 ) and
C*  other region.
C*  Also, the endpoints of the integral region is adjusted
C*  such that the point "q = p" correspond to simp2.
C*
CC**********************************************************************
 
      DOUBLE PRECISION FUNCTION  LNT1FI(N)
 
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
      REAL   *8       NTGRD1, NTGRD2, NTGRD3
C--
      COMMON /KFLAGS/ KFNOFI 
      INTEGER*4       KFNOFI
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
      DQ = ALAMB/50
      IF ( (N.LE.2) .OR. (N.GT.98) ) THEN
          WRITE(*,*)  'INVALID N'
          STOP
      END IF
      SUM = 0
C--
C  Adjust end point.
C--
      IF  ( MOD(N,2) .EQ. 1 )  THEN
             M = 1
      ELSE
             M = 2
      END IF
C--
C  Integration loop starts here.
C--
100   CONTINUE
         IF ( ABS(M-N) .LE. 2 ) THEN
               SIMP1 = NTGRD1(M-1,N) + NTGRD2(M-1,N)
               SIMP2 = NTGRD1(M,N)   + NTGRD2(M,N)
               SIMP3 = NTGRD1(M+1,N) + NTGRD2(M+1,N)
               SUM = SUM + ( SIMP1 + 4*SIMP2 + SIMP3 )
C>>>
C        WRITE(20,*) M-1, PS(M-1), SIMP1
C        WRITE(20,*) M,   PS(M),    SIMP2
C        WRITE(20,*) M+1, PS(M+1), SIMP3
C>>>
               M = M + 2
         ELSE
               SIMP1 = NTGRD1(M-1,N) + NTGRD3(M-1,N)
               SIMP2 = NTGRD1(M,N)   + NTGRD3(M,N)
               SIMP3 = NTGRD1(M+1,N) + NTGRD3(M+1,N)
               SUM   = SUM + ( SIMP1 + 4*SIMP2 + SIMP3 )
C>>>
C        WRITE(20,*) M-1, PS(M-1), SIMP1
C        WRITE(20,*) M,   PS(M),    SIMP2
C        WRITE(20,*) M+1, PS(M+1), SIMP3
C>>>
               M = M + 2
         END IF
      IF ( M .LE.98 )                            GO TO 100
C--
      LNT1FI = SUM*DQ/6
C--
C  That's it.
C--
      RETURN
      END
