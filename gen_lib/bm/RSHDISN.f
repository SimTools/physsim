C*
C* This version for the beam energy width effect studies.
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
      PARAMETER  ( MXxDAT = 100 )
      REAL   *8  BMDATA(2,0:MXxDAT)
      DATA ( BMDATA(I,0),I=1,2) /  2*0.D0 /
C      INCLUDE 'include/xF0070.inc'
C      INCLUDE 'include/xF0050.inc'
C      INCLUDE 'include/xF0035.inc'
C      INCLUDE 'include/xF0020.inc'
C      INCLUDE 'include/xF0010.inc'
C      INCLUDE 'include/xF0005.inc'
C      INCLUDE 'include/xY1070.inc'
C      INCLUDE 'include/xY1050.inc'
C      INCLUDE 'include/xY1035.inc'
C      INCLUDE 'include/xY1020.inc'
C      INCLUDE 'include/xY1010.inc'
      INCLUDE 'include/xY1005.inc'
C      INCLUDE 'include/xY1000.inc'
C
C========< Entry Point >================================================
C
C--
C  Decide which bin X belongs to.
C--
      CALL UDSRCH(MXxDAT,2,1,BMDATA(1,1),DBLE(X),IP)
C--
C  Calculate Y.
C     FRAC = MAX(1.E-6,1+5*SGEB-FRAC)
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
      Y   = 1 + SGEB - EXP(-VAL)
C--
C  That's it.
C--
      RETURN
      END

