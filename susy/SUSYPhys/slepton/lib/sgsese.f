C*
C*    AMZ       :(R*8 ): Z  mass
C*    AM (   1) :(R*8 ): se- mass
C*       (   2) :(R*8 ): se+ mass
C*    AMX(   i) :(R*8 ): i-th neutralino mass
C*    CX(1,1,i) :(C*16): ino-e(L)-se(L) coupling
C*      (2,1,i) :(C*16): ino-e(R)-se(L) coupling
C*      (1,2,i) :(C*16): ino-e(L)-se(R) coupling
C*      (2,2,i) :(C*16): ino-e(R)-se(R) coupling
C*    CGE(   1) :(R*8 ): e-gamma-e(L) coupling
C*       (   2) :(R*8 ): e-gamma-e(R) coupling
C*    CZE(   1) :(R*8 ): e-Z-e(L) coupling
C*       (   2) :(R*8 ): e-Z-e(R) coupling
C*    CGSE(  1) :(R*8 ): se-gamma-se(L) coupling
C*        (  2) :(R*8 ): se-gamma-se(R) coupling
C*    CZSE(  1) :(R*8 ): se-Z-se(L) coupling
C*        (  2) :(R*8 ): se-Z-se(R) coupling
C*
 
      SUBROUTINE SGSESE(MODE,AMZ,AM,AMX,CX,CGE,CZE,CGSE,CZSE,
     .                  RS,POLE,CSTH, SG)
 
      IMPLICIT   REAL*8  ( A-H, O-Z )
      INTEGER*4  MODE
      REAL   *8  RS, POLE, CSTH, AM(2), AMX(4),
     .           CGE(2), CZE(2), CGSE(2), CZSE(2), SG
      COMPLEX*16 CX(2,2,4)
C--
      COMPLEX*16 C
      INTEGER*4  IDSE(2,4), IHSE(2,4)
      DATA IDSE  /  1, 1,   2, 2,   2, 1,   1, 2 /
      DATA IHSE  / -1,-1,  +1,+1,  +1,-1,  -1,+1 /
C--
      DATA NCALL /0/
C--
C  Statement function.
C--
      BETA(X1,X2) = SQRT( 1 - 2*(X1+X2) + (X1-X2)**2 )
C
C========< Entry Point >================================================
C
C--
C  Constants.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL  = 1
         xPI    = ACOS(-1.D0)
         xGV2PB = 3.893796623D8
         FAC    = xGV2PB/(64*xPI*xPI)
         AMZ2   = AMZ**2
      ENDIF
C--
C  Set independent variables.
C--
      AM1  = AM(1)
      AM2  = AM(2)
      IF ( RS.LE.AM1+AM2 ) THEN
         SG = 0
         R  = 0
         RETURN
      ELSE
         BT   = BETA((AM1/RS)**2,(AM2/RS)**2)
      ENDIF
      S     = RS*RS
      P1    = BT*RS/2
      E1    = SQRT(AM1**2+P1**2)
      AI    = AM1**2 - RS*E1
      BI    = RS*P1
C--
C  Branch on the final state chirality combinations.
C     MODE = 1 : LL (SE+,SE-)
C          = 2 : RR
C          = 3 : LR
C          = 4 : RL
C--
      IE  = IDSE(1,MODE)
      IP  = IDSE(2,MODE)
      IHE = IHSE(1,MODE)
      IHP = IHSE(2,MODE)
C--
C  LL/RR cases.
C--
      IF ( MODE.LE.2 ) THEN
         CV    = ( CZE(2) + CZE(1) )/2
         CA    = ( CZE(2) - CZE(1) )/2
         CZ    = CZSE(IE)
C--
         A     = CGE(1)*CGSE(1) + CV*CZ* S/(S-AMZ2)
         B     =                  CA*CZ* S/(S-AMZ2)
C--
         T2SSL = 2*(A-B)**2 * 2*( 1 - CSTH**2/3 )
         T2SSR = 2*(A+B)**2 * 2*( 1 - CSTH**2/3 )
C--
         C     = 0
         DO 10 IX = 1, 4
            C = C + CX(IE,IE,IX)*CONJG(CX(IP,IP,IX))
     .               *S*ANG2(AI,BI,AMX(IX),CSTH)
10       CONTINUE
C--
         T2STL = (1-IHE)*DBLE(C)*(A-B)
         T2STR = (1+IHE)*DBLE(C)*(A+B)
C--
         C     = 0
         DO 20 IX = 1, 4
            DO 2 JX = 1, IX
               IF ( IX.EQ.JX ) THEN
                  SYM = 1
               ELSE
                  SYM = 2
               ENDIF
               C = C + CX(IE,IE,IX)*CONJG(CX(IP,IP,IX))
     .                *CX(IE,IE,JX)*CONJG(CX(IP,IP,JX))
     .                *S*S*ANG3(IX,JX,AI,BI,AMX(IX),AMX(JX),CSTH)
     .                *SYM
2           CONTINUE
20       CONTINUE
C--
         T2TTL = (1-IHE)*C/4
         T2TTR = (1+IHE)*C/4
C--
         T2L = ( T2SSL + T2STL + T2TTL )*BT**2/4
         T2R = ( T2SSR + T2STR + T2TTR )*BT**2/4
C--
C  LR/RL cases.
C--
      ELSE
         C     = 0
         DO 30 IX = 1, 4
            DO 3 JX = 1, IX
               IF ( IX.EQ.JX ) THEN
                  SYM = 1
               ELSE
                  SYM = 2
               ENDIF
               C = C + CX(IE,IE,IX)*CONJG(CX(IP,IP,IX))
     .                *CX(IE,IE,JX)*CONJG(CX(IP,IP,JX))
     .                *S*AMX(IX)*AMX(JX)
     .                *ANG4(IX,JX,AI,BI,AMX(IX),AMX(JX),CSTH)
     .                *SYM
3           CONTINUE
30       CONTINUE
C--
         T2L = (1-IHE)*C/2/2
         T2R = (1+IHE)*C/2/2
      ENDIF
C--
C  Calculate total coss section.
C--
      PL   = (1-POLE)/2
      PR   = (1+POLE)/2
      SG   = (FAC/S)*(T2L*PL+T2R*PR)*BT*2*xPI
C--
C  That's it.
C--
      RETURN
      END
 
      REAL*8 FUNCTION ANG2(A,B,AMX,Y)
 
      IMPLICIT   REAL*8 ( A-H, O-Z )
      REAL*8    A, B, AMX, Y
C
C========< Entry Point >================================================
C
      X    = (A-AMX**2)/B
      ANG2 = ( 2*X*Y + (1-X)*(1+X)*LOG((X+Y)/(X-Y)) )/B
      RETURN
      END
 
      REAL*8 FUNCTION ANG3(I,J,A,B,AMI,AMJ,Y)
 
      IMPLICIT   REAL*8 ( A-H, O-Z )
      INTEGER*4  I, J
      REAL   *8  A, B, AMI, AMJ, Y
C
C========< Entry Point >================================================
C
      IF ( I.EQ.J ) THEN
         X    = (A-AMI**2)/B
         ANG3 = ( -Y*( 1 + (1-X)*(1+X)/((Y-X)*(Y+X)) )
     .               + X*LOG((X+Y)/(X-Y)) )*2/B**2
      ELSE
         XI   = (A-AMI**2)/B
         XJ   = (A-AMJ**2)/B
         ANG3 = ( -2*Y + ((1-XI)*(1+XI)/(XJ-XI))*LOG((XI+Y)/(XI-Y))
     .                 + ((1-XJ)*(1+XJ)/(XI-XJ))*LOG((XJ+Y)/(XJ-Y)) )
     .          /B**2
      ENDIF
      RETURN
      END
 
      REAL*8 FUNCTION ANG4(I,J,A,B,AMI,AMJ,Y)
 
      IMPLICIT   REAL*8 ( A-H, O-Z )
      INTEGER*4  I, J
      REAL   *8  A, B, AMI, AMJ, Y
C
C========< Entry Point >================================================
C
      IF ( I.EQ.J ) THEN
         X    = (A-AMI**2)/B
         ANG4 = -2/((Y-X)*(Y+X))/B**2
      ELSE
         XI   = (A-AMI**2)/B
         XJ   = (A-AMJ**2)/B
         ANG4 = ( LOG((XI+Y)/(XI-Y))/(XJ-XI)
     .          + LOG((XJ+Y)/(XJ-Y))/(XI-XJ) )/B**2
      ENDIF
      RETURN
      END
