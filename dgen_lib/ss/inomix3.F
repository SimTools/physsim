C----------------------------------------------------------------------C
      SUBROUTINE INOMIX3( ZM, WM, SIN2W, ALPHE, 
     &                   HM1, M1,M2, TANB,
     &                   SWM, SZM,
     &                   OL, OR, ER, PHIL, PHIR,
     &                   ON, UN, ETA,
     &                   IFLG  )
C----------------------------------------------------------------------C
C=======================================================================
C   Diagonalize mass matrices of charged(neutral) gaugino & higgsino
C   and then calculate some constants such as couplings etc.
C  Ref.   Hikasa-san's note
C
C <<< Inputs & Output >>>
C    In:
C     IFLG = 0 : print out diagonalization results
C          .ne.0   no prints
C    Out:
C     IFLG < 0 : error
C     IFLG > 0 : OK
C
C <<< Inputs >>>  ---> all inputs are real*4
C ( standard model constants )
C     1) WM    : mass of a W-boson
C     2) SIN2W : weak mixing angle
C     3) ALPHE : fine structure constant
C     4) ALPHS : strong coupling constant
C ( free parameters )
C     1) HM1   : the higgsino mass term  > 0
C     2) G3MAS : gaugino mass for SU(3)  > 0
C     3) TANB=V2/V1 : ratio of the v.e.v. of two higgs doublets
C     ( gaugino mass relation at the weak scale )
C        G2MAS = (ALPHE/ALPHS/SIN2W) * G3MAS
C        G1MAS = (5/3)*TAN2W * G2MAS
C
C <<< outputs >>>
C REAL*8  ZM                  Z boson mass
C REAL*8  SWM(2), SZM(4)      masses of charginos & neutralinos ( > 0 )
C REAL*8  OL(2,2), OR(2,2), ER      matrices etc. for diagonalization
C REAL*8  PHIL, PHIR                mixing angles
C REL*4   ON(4,4)                   matrices etc. for diagonalization
C COMPLEX*16 UN(4,4), ETA(4)
C=============
C REAL*8  MIXCHG(2,2), MIXNEU(4,4)  mixing matrices
C=======================================================================

      IMPLICIT REAL*8 ( A-H, O-Z )
      PARAMETER ( PI = 3.1415927D0 )

      REAL*8  MIXCHG(2,2), MIXNEU(4,4)
      REAL*8  OL(2,2), OR(2,2), ER
      REAL*8  ON(4,4)
      COMPLEX*16 ETA(4), UN(4,4)
      REAL*8  SWM(2), SWM2(2)
      REAL*8  SZM(4), SZM1(4)
      REAL*8  M1,M2
C  for SSL-2
      REAL*8  WRKARY(10), EV(4,4), EI(4), VW(10), EV1(4,4), EI1(4)
Cxxx      REAL*8  WRKARY(10), EV(4,4), EI(4), VW(10), EV1(4,4), EI1(4)
      INTEGER*4 IORDER(4)

      REAL*8 TSTMAT(4,4)
C-----------------------------------------------------------------------
C check
      IF( G3MAS.LT.0 .OR. HM1.LT.0 ) THEN
       IFLG =-1
       RETURN
      ENDIF

      SQRT2 = SQRT(2.D0)

      COS2W = 1.D0-SIN2W
      TAN2W = SIN2W/COS2W
      SINW  = SQRT(SIN2W)
      COSW  = SQRT(COS2W)
      TANW  = SQRT(TAN2W)
      ZM    = WM/COSW
      WM2   = WM**2
      ZM2   = ZM**2
      HM12  = HM1**2
C  0 < beta < pi
C     BETA  = ATAN2( ABS(TANB), SIGN(1.,TANB) )
C     COSB  = COS(BETA)
C     SINB  = SIN(BETA)
C     COSDB = (COSB-SINB)*(COSB+SINB)
C     SINDB = 2.*SINB*COSB
C     SIN2B = SINB**2
C     COS2B = COSB**2
      TAN2B = TANB**2
      COS2B = 1.D0/(1.D0+TAN2B)
      COSB  = SIGN( SQRT(COS2B), TANB )
      SIN2B = 1.D0-COS2B
      SINB  = SQRT(SIN2B)
      COSDB = (COSB-SINB)*(COSB+SINB)
      SINDB = 2.D0*SINB*COSB

C  constants calculated by using of the given values
C    gaugino masses
      G2MAS = M2
      G1MAS = M1

      G2MAS2 = G2MAS**2
      G1MAS2 = G1MAS**2

C  chargino mass matrix ( 2*2 )
      MIXCHG(1,1) = G2MAS
      MIXCHG(2,1) = SQRT2*WM*SINB
      MIXCHG(1,2) = SQRT2*WM*COSB
      MIXCHG(2,2) = HM1

C  neutralino mass matrix ( 4*4 )
      MIXNEU(1,1) = G1MAS
      MIXNEU(2,1) = 0.D0
      MIXNEU(2,2) = G2MAS
      MIXNEU(3,1) = -ZM* SINW*COSB
      MIXNEU(3,2) =  ZM* COSW*COSB
      MIXNEU(3,3) = 0.D0
      MIXNEU(4,1) =  ZM* SINW*SINB
      MIXNEU(4,2) = -ZM* COSW*SINB
      MIXNEU(4,3) = -HM1
      MIXNEU(4,4) = 0.D0

      MIXNEU(1,2) = MIXNEU(2,1)
      MIXNEU(1,3) = MIXNEU(3,1)
      MIXNEU(1,4) = MIXNEU(4,1)
      MIXNEU(2,3) = MIXNEU(3,2)
      MIXNEU(2,4) = MIXNEU(4,2)
      MIXNEU(3,4) = MIXNEU(4,3)

C  calculation for chargino
C   mass eigenvalues
      DELTC2 = (G2MAS2+HM12 +2*WM2 )**2
     &        -4.D0*( G2MAS*HM1 -WM2*SINDB)**2
CCCC
      IF( DELTC2 .LT.0 ) THEN
        WRITE(6,*) 'G3MAS=',G3MAS,' HM1=',HM1, ' TANB=',TANB
        WRITE(6,*) 'G2MAS=',G2MAS,' WM=',WM, ' SINDB=',SINDB
        IFLG = -1
        RETURN
      ENDIF
      DELTC  = SQRT(DELTC2)
      SWM2(1) = 0.5D0 *(G2MAS2 +HM12 + 2.D0*WM2-DELTC)
C     SWM2(2) = 0.5 *(G2MAS2 +HM12 + 2.*WM2+DELTC)
      SWM2(2) = SWM2(1) + DELTC
      SWM(1) = SQRT(SWM2(1))
      SWM(2) = SQRT(SWM2(2))
      TANDPL = 2.D0*SQRT2*WM*(G2MAS*COSB+HM1*SINB)
     &         /(G2MAS2-HM12-2.D0*WM2*COSDB)
      TANDPR = 2.D0*SQRT2*WM*(G2MAS*SINB+HM1*COSB)
     &         /(G2MAS2-HM12+2.D0*WM2*COSDB)
      TANPL = ( SWM2(1) -G2MAS2 -2.D0*WM2*SIN2B)
     &         /SQRT2/WM/(G2MAS*COSB+HM1*SINB)
      TANPR = ( SWM2(1) -G2MAS2 -2.D0*WM2*COS2B)
     &         /SQRT2/WM/(G2MAS*SINB+HM1*COSB)
C  0 < phi_L < pi
      PHIL  = ATAN2( ABS(TANPL), SIGN(1.D0,TANPL) )
      COSPL = COS(PHIL)
      SINPL = SIN(PHIL)
      OL(1,1) =  COSPL
      OL(2,1) = -SINPL
      OL(1,2) =  SINPL
      OL(2,2) =  COSPL

C  take 0 < phi_R < pi and check the lighter mass eigen value
      PHIR  = ATAN2( ABS(TANPR), SIGN(1.D0,TANPR) )
      COSPR = COS(PHIR)
      SINPR = SIN(PHIR)
      OR(1,1) =  COSPR
      OR(2,1) = -SINPR
      OR(1,2) =  SINPR
      OR(2,2) =  COSPR
C diagonalize mass matrix with OL & OR and check the 1st eigen value
C and then re-define phi_R if the 1st eigen value is less than 0.
C ( the 1st eigen value corresponds to the lighter chargino )
      AM1=0.D0
      DO 5  I=1,2
      DO 5  J=1,2
5     AM1 = AM1 + OR(1,I)*MIXCHG(I,J)*OL(1,J)

      IF( AM1.LT.0 ) THEN
        SINPR = -SINPR
        COSPR = -COSPR
      ENDIF
      PHIR  = ATAN2( SINPR, COSPR )
      ER = SIGN(1.D0, G2MAS*HM1-WM2*SINDB)
      OR(1,1) =  COSPR
      OR(2,1) = -SINPR*ER
      OR(1,2) =  SINPR
      OR(2,2) =  COSPR*ER

      IF( IFLG .EQ.0 ) THEN
        WRITE(6,*) '%%%%% diagonalization results for chargino'
        WRITE(6,*) '   mixing matrix'
        WRITE(6,102) '  M_chargino=', (MIXCHG(1,I),I=1,2)
        WRITE(6,102) '             ', (MIXCHG(2,I),I=1,2)
        WRITE(6,101) ' SWM=',SWM
        WRITE(6,103) ' TANPL=',TANPL, ' TANPR=',TANPR,'  ER=',ER
        WRITE(6,103) ' TANDPL=',TANDPL, ' TANDPR=',TANDPR
        WRITE(6,103) ' PHIL =',PHIL/PI*180.D0,
     &               ' PHIR =',PHIR/PI*180.D0
        WRITE(6,102) '  OL=', ( OL(1,I),I=1,2)
     &              ,'  OR=', ( OR(1,I),I=1,2)
        WRITE(6,102) '     ', ( OL(2,I),I=1,2)
     &              ,'     ', ( OR(2,I),I=1,2)
        DO 1001 I=1,2
        DO 1002 J=1,2
        TSTMAT(I,J)=0
        DO 1003 K=1,2
        DO 1004 L=1,2
        TSTMAT(I,J)=TSTMAT(I,J)+ OR(I,K)*MIXCHG(K,L)*OL(J,L)
1004    CONTINUE
1003    CONTINUE
1002    CONTINUE
1001    CONTINUE
        WRITE(6,*) ' TSTMAT(I,J)= OR(I,K)*MIXCHG(K,L)*OL(J,L)'
        WRITE(6,201) '     ', (TSTMAT(1,J),J=1,2)
        WRITE(6,201) '     ', (TSTMAT(2,J),J=1,2)
      ENDIF
101   FORMAT( 1H ,    A, 2(F14.7,1X ) )
102   FORMAT( 1H ,    A, 2(F14.7,1X ),
     &           4X , A, 2(F14.7,1X )  )
103   FORMAT( 1H ,    4( A,F14.7 ) )

C  diagonalize neutralino mass matrix
      WRKARY(1) = MIXNEU(1,1)
      WRKARY(2) = MIXNEU(2,1)
      WRKARY(3) = MIXNEU(2,2)
      WRKARY(4) = MIXNEU(3,1)
      WRKARY(5) = MIXNEU(3,2)
      WRKARY(6) = MIXNEU(3,3)
      WRKARY(7) = MIXNEU(4,1)
      WRKARY(8) = MIXNEU(4,2)
      WRKARY(9) = MIXNEU(4,3)
      WRKARY(10)= MIXNEU(4,4)

C       ---- SSL-2 subroutine ----
      CALL DSEIG1( WRKARY, 4, EI1, EV1, 4, MEIG, VW, ICON )
Cxxx      CALL SEIG1( WRKARY, 4, EI1, EV1, 4, MEIG, VW, ICON )
C       EI1(i)   : eigen values
C       EV1(*,i) : eigen vectors corresponding to EI1(i)
      IF( ICON .GT. 0 ) THEN
C      WRITE(6,*) ' Error in MIXMAT : SSL-2 Error (SEIG1)',
C    &            ' Condition code=', ICON
C      STOP
       IFLG = -ICON
       RETURN
      ENDIF
C       re-arrange eigen-values in increasing (absolute value) order
      DO 10 I=1,4
       SZM(I) = ABS(EI1(I))
10    CONTINUE
      CALL SWPELM( SZM, IORDER, SZM1, 4 )
      IF( IFLG .EQ.-1) THEN
        WRITE(6,*) '%%%%% diagonalization results for neutralino'
        WRITE(6,*) '   mixing matrix'
        WRITE(6,201) '  M_neutralino=', (MIXNEU(1,I),I=1,4)
        WRITE(6,201) '               ', (MIXNEU(2,I),I=1,4)
        WRITE(6,201) '               ', (MIXNEU(3,I),I=1,4)
        WRITE(6,201) '               ', (MIXNEU(4,I),I=1,4)
        WRITE(6,201) '     SZM   =', SZM
        WRITE(6,201) '     SZM1  =', SZM1
        WRITE(6,201) '     EI1   =', EI1
        WRITE(6,202) '     IORDER=', IORDER
        WRITE(6,201) '     EV1(1,*)=',(EV1(1,I),I=1,4)
        WRITE(6,201) '     EV1(2,*)=',(EV1(2,I),I=1,4)
        WRITE(6,201) '     EV1(3,*)=',(EV1(3,I),I=1,4)
        WRITE(6,201) '     EV1(4,*)=',(EV1(4,I),I=1,4)
201     FORMAT(1H , A, 4(F14.7,2X) )
202     FORMAT(1H , A, 4(I14  ,2X) )
      ENDIF
C       treat the degenerate case
14    CONTINUE
      DO 11 I=1,4
       IOK = 0
       DO 12 J=1,4
         IF( IORDER(J) .EQ. I ) IOK=1
12     CONTINUE
       IF( IOK .NE. 1 ) THEN
         DO 13 K=1,4
           IF( IORDER(K)+1 .EQ. I ) THEN
             IORDER(K) = I
             GO TO 14
           ENDIF
13       CONTINUE
       ENDIF
11    CONTINUE

      DO 15 I=1,4
C        J : order index in increasing (abs.value) order.
        J = 5-IORDER(I)
        EI(J) = EI1(I)
        SZM(J)= EI1(I)
        IF( SZM(J).GT.0.D0 ) THEN
          ETA(J) = (1.D0, 0.D0)
        ELSE
          ETA(J) = (0.D0, 1.D0)
        ENDIF
        DO 16 K=1,4
          EV(K,J) = EV1(K,I)
          ON(J,K) = EV1(K,I)
16      CONTINUE
15    CONTINUE

      DO 17 I=1,4
      DO 17 J=1,4
17    UN(I,J)=ON(I,J)*ETA(I)

      DO 51 I=1,2
51    SWM(I) = ABS(SWM(I))
      DO 52 I=1,4
52    SZM(I) = ABS(SZM(I))

      IF( IFLG.EQ.0 ) THEN
        WRITE(6,*) '%%%%% diagonalization results for neutralino'
        WRITE(6,*) '   mixing matrix'
        WRITE(6,301) '  M_neutralino=', (MIXNEU(1,I),I=1,4)
        WRITE(6,301) '               ', (MIXNEU(2,I),I=1,4)
        WRITE(6,301) '               ', (MIXNEU(3,I),I=1,4)
        WRITE(6,301) '               ', (MIXNEU(4,I),I=1,4)
        WRITE(6,301) '     EI1   =', EI1
        WRITE(6,202) '     IORDER=', IORDER
        WRITE(6,301) '     SZM   =', SZM
        WRITE(6,*  ) '                    j--->'
        WRITE(6,301) '    EV1(i,j)=',(EV1(1,I),I=1,4)
        WRITE(6,301) '             ',(EV1(2,I),I=1,4)
        WRITE(6,301) '       i     ',(EV1(3,I),I=1,4)
        WRITE(6,301) '             ',(EV1(4,I),I=1,4)
        WRITE(6,*  ) ' '
        WRITE(6,302) '     ON(i,j)=',(ON(1,I),I=1,4),
     &  '    ETA=', ETA(1)
        WRITE(6,302) '             ',(ON(2,I),I=1,4),
     &  '    ETA=', ETA(2)
        WRITE(6,302) '       i     ',(ON(3,I),I=1,4),
     &  '    ETA=', ETA(3)
        WRITE(6,302) '             ',(ON(4,I),I=1,4),
     &  '    ETA=', ETA(4)
        DO 2001 I=1,4
        DO 2002 J=1,4
        TSTMAT(I,J)=0
        DO 2003 K=1,4
        DO 2003 L=1,4
        TSTMAT(I,J) = TSTMAT(I,J) +ON(I,K)*MIXNEU(K,L)*ON(J,L)
2003    CONTINUE
2002    CONTINUE
2001    CONTINUE
          WRITE(6,*) ' TSTMAT(I,J)= ON(I,K)*MIXNEU(K,L)*ON(J,L)'
          WRITE(6,201) '     ', (TSTMAT(1,J),J=1,4)
          WRITE(6,201) '     ', (TSTMAT(2,J),J=1,4)
          WRITE(6,201) '     ', (TSTMAT(3,J),J=1,4)
          WRITE(6,201) '     ', (TSTMAT(4,J),J=1,4)
        ENDIF
301       FORMAT(1H , A,  4(F14.7,2X) )
302       FORMAT(1H , A,  4(F14.7,2X), A, 2F3.0 )
CCCCCCCCCCCC
CCCCCCCCCCCC

      IFLG=1

      RETURN
      END
