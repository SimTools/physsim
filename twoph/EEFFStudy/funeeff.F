#define DRN ran
C....................................................................
C
C           E(+) E(-) INTO E(+) E(-) Qbar Q
C                                             JANUARY 29, 1987
C                            modified     Feb. 29, 1996
C           EXACT CALCULATION OF 6 DIAGRAMS
C      CROSS SECTION = FUNC IN NB=10**(-33) CM**2
C      THE FACTOR 16*E**8*CHARGE**8 IN TAU IS INCLUDED IN FACT.
C
C      THE VARIABLES X(4) AND X(5) ARE SCALED BY E.
C
C      Q12=-(P1-K1)**2/E**2 > 0,   Q22=-(P2-K2)**2/E**2 > 0
C
C      INDEPENDENT VARIABLES ARE Q12,Q22,Q10,Q20, Q32,COSTH AND PHI
C MOMENTA AND ANGLES CUT ACCORDINBG TO THE AMY DETECTOR CONDITTION.
C......................................................................

      FUNCTION FUNC(X)
 
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(10)
C-- 
#include "CUTS.inc"
C--
#include "KINEM.inc"
#include "KINEM2.inc"
#include "KINEMS.inc"
#include "BSHUFL.inc"
C--
      DIMENSION OMEGAX(2), AKX(2)
#ifdef MACHTEN_PPC
      REAL*8   QEXTD
      REAL*8   TMP
      QEXTD(TMP) = TMP
#endif
      EXTERNAL   XHFILL
C
C =====< Entry Point >==================================================
C
 
          ALQ32 = X(1)
         XCOSTH = X(2)
            PHI = X(3)
          ALQ10 = X(4)
          ALQ20 = X(5)
          ALQ12 = X(6)
          ALQ22 = X(7)
 
 
C-- CUT CONDITIONS ---------------------------------------------
C     PCUT=0.0
C     CCUT=1.0
C---------------------------------------------------------------
C     PCUT=5.0/E
C     CCUT=0.7071
C---------------------------------------------------------------
 
      FUNC=0.
      IBX =0
      IBCOS=0
      FUNCZ(1,1)=0.
      FUNCZ(2,1)=0.
      FUNCZ(1,2)=0.
      FUNCZ(2,2)=0.
      ICOS=0
C     ICOS=2 WHEN THE PHASE SPACE OF COSTH IS DEVIDED INTO TWO PART.
 
C     CHANGE OF THE VARIABLE FROM ALOG(Q10) TO Q10=1-OMEGA1 ETC
 
 
C Q10U WAS REPLACED TO APPLY THE MOMENTUM CUT ON TAGGED PARTICLE.
C                             Y.KURIHARA       13/JUN/89
      Q10U=MIN((1.-EME),(1.-PCUT),(1.-PTAGC))
C     Q10U=MIN((1.-EME),(1.-PCUT))
      Q10L=UME*(UME+EME)
      IF ((Q10U-Q10L).LT.0.0) RETURN
      Q10=Q10L*(Q10U/Q10L)**ALQ10
      OMEGA1=1.0D0-Q10
      AK1=DSQRT((OMEGA1-EME)*(OMEGA1+EME))
 
      A=4.*Q10+EME2
      B=(-2.*Q10+(2.-EME2+2.*UME2))*Q10+2.*UME2
      D=(Q10-UME2)**2-UME2*EME2
      Q20U=1.-EME
      IF (Q10.LT.(EME-EME2+2.*UME2)/(2.-EME)) Q20U=(B+2.*AK1*
     1    DSQRT(D))/A
      Q20U=MIN(Q20U,1.-PCUT)
      Q20L=(B-2.*AK1*DSQRT(D))/A
      IF ((Q20U-Q20L).LT.0.0) RETURN
      Q20=Q20L*(Q20U/Q20L)**ALQ20
      OMEGA2=1.0D0-Q20
      AK2=DSQRT((OMEGA2-EME)*(OMEGA2+EME))
 
      AJAC2=Q20*Q10*LOG(Q20U/Q20L)*LOG(Q10U/Q10L)
 
 
C     CHANGE OF THE VARIABLE FROM ALQ12 TO Q12 ETC
 
C      YY1=DLOG(OMEGA1+0.7071*BETHA*AK1-EME2)
C      YY2=DLOG(OMEGA2+0.7071*BETHA*AK2-EME2)
C      YZ1=2.*(YY1-DLOG(Q10*EME))
C      YZ2=2.*(YY2-DLOG(Q20*EME))
C      Q12=2.*EME2*Q10**2/(OMEGA1+BETHA*AK1-EME2)
C     1               *((OMEGA1+BETHA*AK1-EME2)/(EME*Q10))**(2.*ALQ12)
C      Q22=2.*EME2*Q20**2/(OMEGA2+BETHA*AK2-EME2)
C     1               *((OMEGA2+BETHA*AK2-EME2)/(EME*Q20))**(2.*ALQ22)
C      Q12U=2.*(-EME2+OMEGA1+CCUT*BETHA*AK1)
C      Q12L=((1.-CCUT**2)*(OMEGA1-EME2)**2+(CCUT*EME*Q10)**2)*4./Q12U
C      Q12L=2.*(-EME2+OMEGA1-CCUT*BETHA*AK1)
CCCCCCCC FOR TAGGING    CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C e+
C       Q12U=2.*(-EME2+OMEGA1+CTAGU_1*BETHA*AK1)
C       Q12L=2.*(-EME2+OMEGA1+CTAGL_1*BETHA*AK1)
       Q12U=2.*(-QEXTD(EME2)+QEXTD(OMEGA1)+QEXTD(CTAGU_1*BETHA*AK1))
       Q12L=2.*(-QEXTD(EME2)+QEXTD(OMEGA1)+QEXTD(CTAGL_1*BETHA*AK1))
       IF( ABS(CTAGL_1).EQ.ABS(CTAGU_1) ) THEN
         Q12L=( (1.-CTAGL_1**2)*(OMEGA1-EME2)**2 
     .          + (CTAGL_1*EME*Q10)**2   )*4./Q12U
       ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C e-
C       Q22U=2.*(-EME2+OMEGA2-CTAGL_2*BETHA*AK2)
C       Q22L=2.*(-EME2+OMEGA2-CTAGU_2*BETHA*AK2)
       Q22U=2.*(-QEXTD(EME2)+QEXTD(OMEGA2)-QEXTD(CTAGL_2*BETHA*AK2))
       Q22L=2.*(-QEXTD(EME2)+QEXTD(OMEGA2)-QEXTD(CTAGU_2*BETHA*AK2))
       IF( ABS(CTAGL_2).EQ.ABS(CTAGU_2) ) THEN
         Q22L=( (1.-CTAGL_2**2)*(OMEGA2-EME2)**2
     .          + (CTAGL_2*EME*Q20)**2   )*4./Q22U
       ENDIF


C       write(6,*) 'q12u=',q12u,' q12l=',q12l
C       write(6,*) 'q22u=',q22u,' q22l=',q22l
       if( q12l.le.0 ) return
       if( q22l.le.0 ) return

       Q12=Q12L*(Q12U/Q12L)**ALQ12
       Q22=Q22L*(Q22U/Q22L)**ALQ22
       AJAC1=LOG(Q12U/Q12L)*LOG(Q22U/Q22L)/(4.*BETHA**2)
 
 
C     IN THE FOLLOWING SINTH1=SIN(THETA1)*AK1, SINTH2=SIN(THETA2)*AK2
 
      IF(Q12*(1.-EME2-Q10-0.25D0*Q12)-EME2*Q10**2.LT.0 ) RETURN
      IF(Q22*(1.-EME2-Q20-0.25D0*Q22)-EME2*Q20**2.LT.0 ) RETURN
      SINTH1=DSQRT(Q12*(1.-EME2-Q10-0.25D0*Q12)-EME2*Q10**2)/BETHA
      SINTH2=DSQRT(Q22*(1.-EME2-Q20-0.25D0*Q22)-EME2*Q20**2)/BETHA
      SINT12=2.*SINTH1*SINTH2
      COSTH1=-(1.-EME2-Q10-0.5D0*Q12)/BETHA
      COSTH2=+(1.-EME2-Q20-0.5D0*Q22)/BETHA
      COSPHI=COS(PHI)
      SINPHI=SIN(PHI)
 
 
C     CHANGE OF VARIABLE FROM ALQ32 TO Q32
 
      PSI0=2.*Q10*Q20+2.*(Q10+0.5D0*Q12)*(Q20+0.5D0*Q22)/BETHA**2
     .    -Q12-Q22
      Q32U=4.*(1.-EME)**2
      Q32L=4.*UME2
C      write(6,*) 'q32u=',q32u,' psi0=',psi0,' sint12=',sint12
      Q32U=DMIN1(Q32U,PSI0+SINT12)
      W2MIN=WMIN*WMIN
      Q32L=DMAX1(Q32L,PSI0-SINT12,W2MIN)
      IF ((Q32U-Q32L).LE.0.) RETURN
      Q32=Q32L*(Q32U/Q32L)**ALQ32
      AJAC4=Q32*LOG(Q32U/Q32L)
 
 
      PSI=PSI0-Q32
      DDD=(SINT12+PSI)*(SINT12-PSI)
      IF (DDD.LE.0.) RETURN
      DDDD=DSQRT(DDD)
      COSPH2=PSI/SINT12
      AJAC5=2./DDDD
      IF(DRN(DUM).GT.0.5D0) DDDD=-DDDD
      SINPH2=DDDD/SINT12
 
 
C     CHANGE OF VARIABLE FROM XCOSTH TO COSTH
 
      A=0.5D0*Q32
      B=-Q10-Q20
       CS=(SINTH1+0.5D0*PSI/SINTH1)*COSPHI+0.5D0*DDDD*SINPHI/SINTH1
       CC=(Q10-Q20+0.5D0*(Q12-Q22))/BETHA
       CG2=B**2-A**2/UME2
 
           COSTHU=1.D0
           COSTHL=-1.D0
 
       IF (CG2.GT.0.) THEN
         CG=DSQRT(CG2)
         CF=CS**2+CC**2
         CD=CF-CG2
               IF (CD.LE.0.) THEN
                    RETURN
               ELSE
                    CD=DSQRT(CD)
                    XMMAX=(CC*CG+CS*CD)/CF
                    XMMIN=(CC*CG-CS*CD)/CF
                    XPMAX=(-CC*CG+CS*CD)/CF
                    XPMIN=(-CC*CG-CS*CD)/CF
 
                        IF (DABS(CC).GT.CG) THEN
                            ICOS=2
 
                            IF (CS.GT.0.) THEN
                               IF (CC.GT.0.) THEN
                                   COSXX1=XMMIN
                                   COSXX2=XPMIN
                                     ELSE
                                       COSXX1=XPMAX
                                       COSXX2=XMMAX
                                          ENDIF
                            ELSE
                                       IF (CC.GT.0.) THEN
                                     COSXX1=XMMIN
                                     COSXX2=XPMIN
                                   ELSE
                                 COSXX1=XPMAX
                                 COSXX2=XMMAX
                               ENDIF
                            ENDIF
 
                        ELSE
 
                             IF (CS.GT.0.) THEN
                                 COSTHU=XMMAX
                                 COSTHL=XMMIN
                                    ELSE
                                 COSTHU=XPMIN
                                 COSTHL=XPMAX
                             END IF
 
                        ENDIF
           ENDIF
 
       ENDIF
 
 
 
      DO 800 II=1,2
      IF (ICOS.EQ.0) THEN
           COSTHU=MIN(COSTHU,CMUCUT)
           COSTHL=MAX(COSTHL,-CMUCUT)
           COSTH=0.5D0*((COSTHU-COSTHL)*XCOSTH+COSTHU+COSTHL)
           AJAC3=0.5D0*(COSTHU-COSTHL)
           IF (AJAC3.LT.0.0) RETURN
      ELSE
           IF (II.EQ.1) THEN
           COSTH=0.5D0*((CMUCUT-COSXX1)*XCOSTH+CMUCUT+COSXX1)
           AJAC3=0.5D0*(CMUCUT-COSXX1)
           IF (AJAC3.LT.0.0) GO TO 800
                ELSE
           COSTH=0.5D0*((CMUCUT+COSXX2)*XCOSTH -CMUCUT+COSXX2)
           AJAC3=0.5D0*(CMUCUT+COSXX2)
           IF (AJAC3.LT.0.0) GO TO 800
           END IF
      END IF
 
           SINTH=DSQRT((1.-COSTH)*(1.+COSTH))
C---- SPRING
           COSTHZ(II)=COSTH
           SINTHZ(II)=SINTH
C-----------
 
      C=CS*SINTH+CC*COSTH
      DEN=(C+B)*(C-B)
      DEN1=1./DEN
      D=A**2+DEN*UME2
C     IF (D.LE.0.0D0) RETURN
      D=DSQRT(ABS(D))
      AKX(1) = (A*C+B*D)*DEN1
      AKX(2) = (A-B*UME)*(A+B*UME)/(A*C+B*D)
      OMEGAX(1)=(A*B+C*D)*DEN1
      OMEGAX(2)=-(A**2+C**2*UME2)/(A*B+C*D)
         IF (A*C.GT.0.) THEN
         AKX(1)=(A-B*UME)*(A+B*UME)/(A*C-B*D)
         AKX(2)=(A*C-B*D)*DEN1
         OMEGAX(1)=-(A**2+C**2*UME2)/(A*B-C*D)
         OMEGAX(2)=(A*B-C*D)*DEN1
         END IF
 
 
C
 
      DO 100 I=1,2
      OMEGA=OMEGAX(I)
      AK=AKX(I)
C----SPRING
      OMEGAZ(I,II)=OMEGA
      AKZ(I,II)=AK
C----------
      IF (AK.LT.0. ) GO TO 100
      IF (OMEGA.GT.(1.-EME)) GO TO 100
      IF (OMEGA.LT.0.) GO TO 100
      AJAC0=0.5D0*AK/D
 
 
 
 
      P3Q1=OMEGA*Q10+AK*(SINTH1*SINTH*COSPHI
     1                   +(Q10+0.5D0*Q12)/BETHA*COSTH)
      P3Q2=0.5D0*Q32-P3Q1
 
 
C.....................................................................
C
C     IN CONTRACTING THE ELECTRON AND POSITRON CURRENTS WITH
C     MUON CURRENT, USE THE TRICK SUCH THAT  EACH TERM CONTAINS
C     EXPLICITLY (1+COSTH1),(1-COSTH2) OR EME2.
C
C......................................................................
 
 
 
C     S1=SS1/BETHA/Q10,   S2=SS2/BETHA/Q20
 
      S1=(0.5D0*Q12/Q10 + EME2)/BETHA
      S2=(0.5D0*Q22/Q20 + EME2)/BETHA
      S1Q2=2.*Q20+0.5D0*Q22-0.5D0*(Q12+Q22+Q32)/Q10
      S2Q1=2.*Q10+0.5D0*Q12-0.5D0*(Q12+Q22+Q32)/Q20
      S1P3=-S1*COSTH*AK-AK*SINTH1*SINTH*COSPHI/Q10
      S2P3= S2*COSTH*AK
     .             -AK*SINTH*0.5*(PSI*COSPHI+DDDD*SINPHI)/SINTH1/Q20
      S1S2=-2.-EME2-0.5D0*(Q22/Q20+Q12/Q10-(Q12+Q22+Q32)/Q10/Q20)
      S1Q1=Q12*(1.-0.5D0*Q10)/Q10
      S2Q2=Q22*(1.-0.5D0*Q20)/Q20
      S12=EME2-Q12*OMEGA1/Q10**2
      S22=EME2-Q22*OMEGA2/Q20**2
 
 
      P1P2=2.*(1.-Q10-Q20)-EME2+0.5D0*Q32
      P1P3=OMEGA1*OMEGA-AK*SINTH1*SINTH*COSPHI
     .    +OMEGA1*AK*(1.+0.5D0*Q12/Q10)*COSTH/BETHA-AK*COSTH*S1
      P2P3=2.*OMEGA-0.5D0*Q32-P1P3
C     P1P3=OMEGA1*OMEGA-AK*SINTH1*SINTH*COSPHI
C    .     +AK*COSTH*(OMEGA1-0.5D0*Q12-EME2)/BETHA
      P1Q2=2.*Q20-0.5D0*(Q12+Q32)
      P2Q1=2.*Q10-0.5D0*(Q22+Q32)
      S1P2=-EME2-2.*Q20-0.5D0*(Q22-(Q22+Q32)/Q10)
      S2P1=-EME2-2.*Q10-0.5D0*(Q12-(Q12+Q32)/Q20)
 
 
      S1P1=EME2-0.5D0*OMEGA1*Q12/Q10
      S2P2=EME2-0.5D0*OMEGA2*Q22/Q20
 
 
          D1=Q22+2.*P3Q2
          D2=Q12+2.*P3Q1
          D3=Q22+2.*P1Q2
          D4=-Q12-Q32-2.*P1Q2
          D5=Q12+2.*P2Q1
          D6=-Q22-Q32-2.*P2Q1
 
 
      OMEGA4=Q10+Q20-OMEGA
      AK4=DSQRT((OMEGA4-UME)*(OMEGA4+UME))
 
      IF (ABS((CC+AK*COSTH)/AK4).GE.CMUCUT) GO TO 100
 
      IF(OMEGA .LE.PMUCUT) GO TO 100
      IF(OMEGA4.LE.PMUCUT) GO TO 100
 
      CALL FUNC0(TAU0)
      TAUA=TAU0/D1**2/D2**2/Q12/Q22
 
 
      CALL FUNC1(TAU1,TAU2)
      CALL FUNC2(TAU3)
      TAUB=(TAU1*Q12/Q22/(D3*D4)**2+TAU2*Q22/Q12/(D5*D6)**2
     .                            +TAU3/(D3*D4*D5*D6))/Q32**2
 
 
      CALL FUNC3(TAU4,TAU5)
      CALL FUNC4(TAU6,TAU7)
      TAUC=((TAU4/Q22/(D3*D4)+TAU5/Q12/(D5*D6))/D1
     .    + (TAU6/Q22/(D3*D4)+TAU7/Q12/(D5*D6))/D2)/Q32
 
 
 
C     TAU=TAUA+TAUB+TAUC
C         PERIPHERAL           BREMSS.        INTERFERENCE  COLLOR
C                 !                 !               !       FACTOR
      TAU= (TAUA*CHARGE**4 + TAUB*CHARGE**2 - TAUC*CHARGE**3)*3
      IF(CHARGE.LT.-1.D0/3.D0-.0001D0)TAU=TAU/3.D0
C
      FUNCX= FACT*TAU*AJAC1*AJAC2*AJAC3*AJAC4*AJAC5*AJAC0*AK
      FUNCZ(I,II)=FUNCX
      FUNC=FUNC+FUNCX
C--------------------------------------------------- HISTOGRAM PACKAGES
                 CALL XHFILL( 'h01',COSTH2/AK2, FUNCX )
                 CALL XHFILL( 'h02',COSTH1/AK1, FUNCX )
                 CALL XHFILL( 'h03',OMEGA2, FUNCX )
                 CALL XHFILL( 'h04',OMEGA1, FUNCX )
                 CALL XHFILL( 'h05',Q32, FUNCX )
                 CALL XHFILL( 'h06',COSTH, FUNCX )
                 CALL XHFILL( 'h07',OMEGA, FUNCX )
                 CALL XHFILL( 'h08',DSQRT(Q32), FUNCX )
                 CALL XHFILL( 'h09',DSQRT(Q32)*E, FUNCX )
C--------------------------------------------------------------------
  100 CONTINUE
 
      IF (ICOS.EQ.0) GO TO 900
  800 CONTINUE
  900 CONTINUE
      RSUM=DRN(DUM)*(FUNCZ(1,1)+FUNCZ(2,1)+FUNCZ(1,2)+FUNCZ(2,2))
      FSUM=0.
      DO 910 II=1,2
      DO 910 I =1,2
      IF(FUNCZ(I,II).LE.0) GOTO 910
      FSUM=FSUM+FUNCZ(I,II)
      IF(RSUM.LE.FSUM) THEN
         IBX=I
         IBCOS=II
C---- ONLY SPRING
         FUNC=FUNCZ(I,II)
C                CALL XHFILL( 'h01',COSTH2/AK2, FUNC )
C                CALL XHFILL( 'h02',COSTH1/AK1, FUNC )
C        CALL XHFILL('h01',OMEGA1,FUNC)
C        CALL XHFILL('h02',OMEGA2,FUNC)
C        CALL DHFILL('hd01', ALQ12,        OMEGAZ(I,II), FUNC)
C        CALL DHFILL('hd02', OMEGAZ(I,II), Q12,          FUNC)
C        CALL DHFILL('hd03', OMEGAZ(I,II), OMEGA1,       FUNC)
         RETURN
      ENDIF
 910  CONTINUE
      RETURN
      END

