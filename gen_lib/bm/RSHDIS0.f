C
C  *INCLUDE UDSRCH
C        CALL ISRINT(1,300.D0)
C        WRITE(20,'(''NEWFRAME'')')
C        WRITE(20,'(''SET ORDER X Y'')')
C        DO 10 I = 0, 100
C           X = 0.01*I
C           CALL RSHDIS(X,1,Y)
C           WRITE(20,'(2E15.8)') Y, X
C  10    CONTINUE
C        WRITE(20,'(''JOIN'')')
C        END
C*
C* This version for bremsstrahlung only.
C*                   _
C*                  / 1
C*    sig_RC(s) =  /   dz b [ z**(b-1)(1+3b/4) - (1-z/2) ] sig((1-z)s)
C*               _/ 0
C*
      SUBROUTINE RSHDIS(X,ITYP,Y)
 
      IMPLICIT   REAL*8 ( A-H, O-Z )
      INTEGER*4  ITYP
      REAL   *4  X, Y
C--
      PARAMETER       ( MXxDAT = 500 )
      COMMON   /RSHISR/ BMDATA(2,0:MXxDAT), BTE
      REAL*8            BMDATA,             BTE
C
C========< Entry Point >================================================
C
C--
C  Decide which bin X belongs to.
C--
      CALL UDSRCH(MXxDAT,2,1,BMDATA(1,1),DBLE(X),IP)
C--
C  Calculate Y.
C     FRAC = MAX(1.E-6,1-FRAC)
C     FRAC = -LOG(FRAC)
C--
      IF ( IP.NE.0 .AND. IP.NE.MXxDAT ) THEN
         XL  = BMDATA(1,IP)
         DX  = BMDATA(1,IP+1) - XL
         DF  = (BMDATA(2,IP+1)-BMDATA(2,IP))/DX
         VAL = BMDATA(2,IP) + (X-XL)*DF
      ELSE
         VAL = BMDATA(2,IP)
      ENDIF
      Y   = SQRT( 1 - VAL**(1/BTE) )
C--
C  That's it.
C--
      RETURN
      END
