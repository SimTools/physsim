C======================================================================C
C subroutines to give couplings of supersymmetric particles
C   Caution !!!!
C     outputs are complex*8
C   Outputs are always those of G_couple of the following Lagrangean
C      Lagrangean = G_couple * (field_i) * (field_j)....
C   Even if field_i and field_j are the same Majorana particles,
C   factor 1/2 is also included. In other word, if L is
C      Lagrangean = (1/2) * G_couple *(Maj)*(Maj)
C   output is given as (1/2)*G_couple.
C======================================================================C
C----------------------------------------------------------------------C
      SUBROUTINE  INOGBS( ALPHE, SIN2W,
     &                    PHIL, PHIR, ER, ON, ETA,
     &                    GNCW, GCCZ, GNNZ )
      IMPLICIT REAL*8 (A-H, O-Z)
C----------------------------------------------------------------------C
C Purpose: give couplings of Ino -Gauge_boson
C Input  : ALPHE  e**2/4/pi         !
C          SIN2W  weak mix. angle   !  S.M. parameters
C          PHIL, PHIR (real*4) mixing angle for chargino
C          ER         (real*4) phase for chargino
C          ON(i,j) (real*4)(i,j=1..4)diagonalizing matrix for chargino
C          ETA(i)  (complex*8) (i=1..4)    phase factor
C Output : GNCW(i,j,k) (i=1..2, j=1..4, k=1..2 )
C          GCCZ(i,j,k) (i=1..2, j=1..2, k=1..2 )
C          GNNZ(i,j,k) (i=1..2, j=1..4, k=1..4 )  ( j =/< k )
C              i=1/left-handed,  i=2/right-handed
C
C----------------------------------------------------------------------C
      REAL*8     ON(4,4)
      COMPLEX*16 ETA(4)
 
      COMPLEX*16 GNCW(2,4,2)
      COMPLEX*16 GCCZ(2,2,2), GV(2,2), GA(2,2)
      COMPLEX*16 GNNZ(2,4,4), CA(4,4), CV(4,4)
      COMPLEX*16 IMGI, CMPLX
      DATA IMGI  /(0.D0,1.D0)/
C----------------------------------------------------------------------C
      DO 1 I=1,4
        IF( ETA(I).NE.1. .AND. ETA(I).NE.IMGI ) THEN
          WRITE(6,*) ' ETA is wrong, ETA=',ETA
          STOP
        ENDIF
1     CONTINUE
 
      PI = ACOS(0.)*2.
      SQRT2 = SQRT(2.)
 
      COS2W = 1.-SIN2W
      COSW  = SQRT(COS2W)
      SINW  = SQRT(SIN2W)
      TANW  = SINW/COSW
      TAN2W = TANW**2
 
      GWEAK = SQRT(ALPHE*4.*PI/SIN2W)
      GZ    = GWEAK/COSW
 
      SINPL = SIN(PHIL)
      SINPR = SIN(PHIR)
      COSPL = COS(PHIL)
      COSPR = COS(PHIR)
      SIN2PL = SINPL**2
      COS2PL = COSPL**2
      SIN2PR = SINPR**2
      COS2PR = COSPR**2
 
C---
C Neutralino-Chargino-W_boson
C---
      DO 10 I=1,4
        GNCW(1,I,1) = ETA(I)*( ON(I,2)*COSPL +ON(I,3)*SINPL/SQRT2 )
     &                *(-GWEAK)
        GNCW(1,I,2) = ETA(I)*(-ON(I,2)*SINPL +ON(I,3)*COSPL/SQRT2 )
     &                *(-GWEAK)
        GNCW(2,I,1) =
     & -GWEAK*CONJG(ETA(I))*( ON(I,2)*COSPR -ON(I,4)*SINPR/SQRT2 )
        GNCW(2,I,2) =
     & -GWEAK*CONJG(ETA(I))*ER*(-ON(I,2)*SINPR -ON(I,4)*COSPR/SQRT2 )
10    CONTINUE
 
C---
C Chargino-Chargino-Z_boson
C---
      GV(1,1) = 0.25*( SIN2PL+ SIN2PR ) - COS2W
      GA(1,1) = 0.25*( SIN2PL- SIN2PR )
 
      GV(2,2) = 0.25*( COS2PL+ COS2PR ) - COS2W
      GA(2,2) = 0.25*( COS2PL- COS2PR )
 
      GV(1,2) = 0.25*( COSPL*SINPL + ER*COSPR*SINPR )
      GA(1,2) = 0.25*( COSPL*SINPL - ER*COSPR*SINPR )
      GV(2,1) = GV(1,2)
      GA(2,1) = GA(1,2)
 
      DO 20 I=1,2
      DO 20 J=1,2
        GCCZ(1,I,J) = GV(I,J)+ GA(I,J)
        GCCZ(2,I,J) = GV(I,J)- GA(I,J)
20    CONTINUE
      DO 21 I=1,2
      DO 21 J=1,2
      DO 21 K=1,2
21    GCCZ(I,J,K) = -GCCZ(I,J,K)*GZ
 
C---
C Neutralino-Neutralino-Z_boson
C---
      DO 30 I=1,4
      DO 30 J=1,4
       COM = ON(I,3)*ON(J,3) - ON(I,4)*ON(J,4)
       CMPLX = ETA(I)*CONJG(ETA(J))
       IF( REAL(CMPLX) .EQ. 0 ) THEN
         CV(I,J) = 0.5*IMGI *IMAG(CMPLX) * COM
         CA(I,J) = 0.
       ELSE
         CV(I,J) = 0.
         CA(I,J) = 0.5 *REAL(CMPLX) * COM
       ENDIF
30    CONTINUE
      DO 31 I=1,4
      DO 31 J=1,4
        GNNZ(1,I,J) = CV(I,J)+CA(I,J)
        GNNZ(2,I,J) = CV(I,J)-CA(I,J)
31    CONTINUE
      DO 32 I=1,2
      DO 32 J=1,4
      DO 32 K=1,4
      IF(J.NE.K ) THEN
      GNNZ(I,J,K) = -GNNZ(I,J,K)*GZ
      ELSE
      GNNZ(I,J,K) = -GNNZ(I,J,K)*GZ/2.
      ENDIF
32    CONTINUE
 
 
      RETURN
      END
C----------------------------------------------------------------------C
      SUBROUTINE  INOFSF( ALPHE, SIN2W, WM,
     &                    EMAS, UMAS, DMAS,
     &                    TANB, PHIL, PHIR, ER, ON, ETA,
     &                    GCESNL, GCNSEL, GCNSER,
     &                    GCUSDL, GCUSDR, GCDSUL, GCDSUR,
     &                    GNNSNL, GNESEL, GNESER,
     &                    GNUSUL, GNUSUR, GNDSDL, GNDSDR )
C----------------------------------------------------------------------C
C purpose : give couplings for Ino -fermion -Sfermion
C Input  : ALPHE  e**2/4/pi
C          SIN2W (real*4)  weak mixing angle
C          WM    (real*4)  W-boson mass
C          EMAS, UMAS, DMAS (real*4) fermion mass
C               ( charged-lepton, up-type, down-type )
C          TANB  (real*4)  SUSY parameter tan(beta)=v2/v1
C          PHIL, PHIR      chargino mix. angle
C          ER              chargino mix. phase
C          ON(i,j) (real*4)(i,j=1..4)diagonalizing matrix for chargino
C          ETA(i)  (complex*8) (i=1..4)    phase factor
C Output :
C    GCESNL(i,j) (i=1..2, j=1..2 ) / chargino(j)-electron-Sneutrino_L
C    GCNSEL(i,j) (i=1..2, j=1..2 ) / chargino(j)-neutrino-Selectron_L
C    GCNSER(i,j) (i=1..2, j=1..2 ) / chargino(j)-neutrino-Selectron_R
C    GCDSUL(i,j) (i=1..2, j=1..2 ) / chargino(j)-d_quark-Su_quark_L
C    GCDSUR(i,j) (i=1..2, j=1..2 ) / chargino(j)-d_quark-Su_quark_R
C    GCUSDL(i,j) (i=1..2, j=1..2 ) / chargino(j)-u_quark-Sd_quark_L
C    GCUSDR(i,j) (i=1..2, j=1..2 ) / chargino(j)-u_quark-Sd_quark_R
C    GNNSNL(i,j) (i=1..2, j=1..4 ) / neutralino(j)-neutrino-Sneutrino_L
C    GNESEL(i,j) (i=1..2, j=1..4 ) / neutralino(j)-electron-Selectron_L
C    GNESER(i,j) (i=1..2, j=1..4 ) / neutralino(j)-electron-Selectron_R
C    GNUSEL(i,j) (i=1..2, j=1..4 ) / neutralino(j)-u_quark-Su_quark_L
C    GNUSER(i,j) (i=1..2, j=1..4 ) / neutralino(j)-u_quark-Su_quark_R
C    GNDSEL(i,j) (i=1..2, j=1..4 ) / neutralino(j)-u_quark-Su_quark_L
C    GNDSER(i,j) (i=1..2, j=1..4 ) / neutralino(j)-u_quark-Su_quark_R
C              i=1/left-handed-coupling,  i=2/right-handed
C----------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-H, O-Z)
      REAL*8 ON(4,4)
      COMPLEX*16 ETA(4)
      COMPLEX*16 GCESNL(2,2), GCNSEL(2,2), GCNSER(2,2),
     &           GCUSDL(2,2), GCUSDR(2,2), GCDSUL(2,2), GCDSUR(2,2)
      COMPLEX*16 GNNSNL(2,4), GNESEL(2,4), GNESER(2,4),
     &           GNUSUL(2,4), GNUSUR(2,4), GNDSDL(2,4), GNDSDR(2,4)
C----------------------------------------------------------------------C
      PI = ACOS(0.)*2.
      SQRT2 = SQRT(2.)
      GWEAK = SQRT(ALPHE*4.*PI/SIN2W)
      COS2W = 1.-SIN2W
      TANW  = SQRT(SIN2W/COS2W)
C  0 < beta < pi
C     BETA  = ATAN2( ABS(TANB), SIGN(1.,TANB) )
C     BETA  = DATAN2( DABS(DBLE(TANB)), DSIGN(1.D0,DBLE(TANB)) )
C     COSB  = COS(BETA)
C     SINB  = SIN(BETA)
      TAN2B = TANB**2
      COS2B = 1./(1.+TAN2B)
      COSB  = SIGN( SQRT(COS2B), TANB )
      SIN2B = 1.-COS2B
      SINB  = SQRT(SIN2B)
      COSPL = COS(PHIL)
      SINPL = SIN(PHIL)
      COSPR = COS(PHIR)
      SINPR = SIN(PHIR)
 
      WMCOSB = WM*COSB
      WMSINB = WM*SINB
      WMSQR2 = WM*SQRT(2.)
      WSQ2CB = WMSQR2*COSB
      WSQ2SB = WMSQR2*SINB
 
C---
C Chargino-fermion-Sfermion
C---
C chargino-electron-Sneutrino
      GCESNL(1,1) = -GWEAK*COSPR
      GCESNL(2,1) =  GWEAK*EMAS*SINPL/WSQ2CB
      GCESNL(1,2) =  GWEAK*ER*SINPR
      GCESNL(2,2) =  GWEAK*EMAS*COSPL/WSQ2CB
C chargino-neutrino-Selectron_L
      GCNSEL(1,1) = -GWEAK*COSPL
      GCNSEL(2,1) =  0.
      GCNSEL(1,2) =  GWEAK*SINPL
      GCNSEL(2,2) =  0.
C chargino-neutrino-Selectron_R
      GCNSER(1,1) =  GWEAK*EMAS*SINPL/WSQ2CB
      GCNSER(2,1) =  0.
      GCNSER(1,2) =  GWEAK*EMAS*COSPL/WSQ2CB
      GCNSER(2,2) =  0.
C chargino-d_quark-Su_quark_L
      GCDSUL(1,1) = -GWEAK*COSPR
      GCDSUL(2,1) =  GWEAK*DMAS*SINPL/WSQ2CB
      GCDSUL(1,2) =  GWEAK*ER*SINPR
      GCDSUL(2,2) =  GWEAK*DMAS*COSPL/WSQ2CB
C chargino-d_quark-Su_quark_R
      GCDSUR(1,1) =  GWEAK*UMAS*SINPR/WSQ2SB
      GCDSUR(2,1) =  0.
      GCDSUR(1,2) =  GWEAK*UMAS*ER*COSPR/WSQ2SB
      GCDSUR(2,2) =  0.
C chargino-u_quark-Sd_quark_L
      GCUSDL(1,1) = -GWEAK*COSPL
      GCUSDL(2,1) =  GWEAK*UMAS*SINPR/WSQ2SB
      GCUSDL(1,2) =  GWEAK*SINPL
      GCUSDL(2,2) =  GWEAK*UMAS*ER*COSPR/WSQ2SB
C chargino-u_quark-Sd_quark_R
      GCUSDR(1,1) =  GWEAK*DMAS*SINPL/WSQ2CB
      GCUSDR(2,1) =  0.
      GCUSDR(1,2) =  GWEAK*DMAS*COSPL/WSQ2CB
      GCUSDR(2,2) =  0.
 
C---
C Neutralino-fermion-Sfermion
C---
      UCHRG= 2./3.
      DCHRG=-1./3.
      DO 10 I=1,4
        GNNSNL(1,I) = CONJG(ETA(I))*(ON(I,2)-ON(I,1)*TANW)
C       GNNSNL(1,I) = CONJG(ETA(I))*(ON(I,2))
        GNNSNL(2,I) = 0.
        GNESEL(1,I) =-CONJG(ETA(I))*(ON(I,2)+ON(I,1)*TANW)
        GNESEL(2,I) = ETA(I)*EMAS*ON(I,3)/WMCOSB
        GNESER(1,I) = CONJG(ETA(I))*EMAS*ON(I,3)/WMCOSB
        GNESER(2,I) = 2.*ETA(I)*ON(I,1)*TANW
C>>> 96/08/01  Bug fixup due to nojirim
C        GNUSUL(1,I) = CONJG(ETA(I))
C     &               *( ON(I,2)+(2.*UCHRG-1.)*ON(I,1)*TANW )
C        GNUSUL(2,I) = ETA(I)*UMAS*ON(I,3)/WMSINB
C        GNUSUR(1,I) = CONJG(ETA(I))*UMAS*ON(I,3)/WMSINB
C        GNUSUR(2,I) =-2.*ETA(I)*UCHRG*TANW*ON(I,1)
        GNUSUL(1,I) = CONJG(ETA(I))
     &               *( ON(I,2)+(2.*UCHRG-1.)*ON(I,1)*TANW )
        GNUSUL(2,I) = ETA(I)*UMAS*ON(I,4)/WMSINB
        GNUSUR(1,I) = CONJG(ETA(I))*UMAS*ON(I,4)/WMSINB
        GNUSUR(2,I) =-2.*ETA(I)*UCHRG*TANW*ON(I,1)
C>>>
        GNDSDL(1,I) = CONJG(ETA(I))
     &               *(-ON(I,2)+(2.*DCHRG+1.)*ON(I,1)*TANW )
        GNDSDL(2,I) = ETA(I)*DMAS*ON(I,3)/WMCOSB
        GNDSDR(1,I) = CONJG(ETA(I))*DMAS*ON(I,3)/WMCOSB
        GNDSDR(2,I) =-2.*ETA(I)*DCHRG*TANW*ON(I,1)
10    CONTINUE
 
      GWSQR2 = GWEAK/SQRT2
      DO 20 I=1,2
      DO 20 J=1,4
        GNNSNL(I,J) = -GWSQR2*GNNSNL(I,J)
        GNESEL(I,J) = -GWSQR2*GNESEL(I,J)
        GNESER(I,J) = -GWSQR2*GNESER(I,J)
        GNUSUL(I,J) = -GWSQR2*GNUSUL(I,J)
        GNUSUR(I,J) = -GWSQR2*GNUSUR(I,J)
        GNDSDL(I,J) = -GWSQR2*GNDSDL(I,J)
        GNDSDR(I,J) = -GWSQR2*GNDSDR(I,J)
20    CONTINUE
 
 
      RETURN
      END
C----------------------------------------------------------------------C
      SUBROUTINE  HGSGBS( ALPHE, SIN2W, WM, TANB, TANA,
     &                    GH0WW, GH0ZZ,
     &                    GH0HCW, GH0PHZ, GHHZ, GHHA,
     &                    GHHWW, GH02W2, GH0HWZ, GH0HWA,
     &                    GHHZZ, GH02Z2, GHHZA, GHHAA )
C----------------------------------------------------------------------C
C purpose : give  Higss- Gauge_boson couplings
C inputs  : ALPHE, SIN2W, WM / S.M. parameters
C           TANB / ratio of vev. ( v2/v1 )
C           TANA / tan of Higgs boson mix. angle
C outputs : GH0WW(i) (i=1..2 ) / Higgs_0(i)-W-W coupling
C           GH0ZZ(i) (i=1..2 ) / Higgs_0(i)-Z-Z coupling
C           GH0HCW(i)(i=1..3 ) / Higgs_0(i)-Higgs_(-)- W(dag.)
C                         (i=3 is for pseudoscaler Higgs )
C           GH0PHZ(i)(i=1..2 ) / Higgs_0(i)-P_Higgs_0 -Z
C           GHHZ               / Higgs_(+) -Higgs_(-) -Z
C           GHHA               / Higgs_(+) -Higgs_(-) -A
C           GHHWW              / Higgs_(+) -Higgs_(-) -W-W
C           GH02W2             / Higgs_0(i) -Higgs_0(i) -W-W
C                                ( same for all H_0(i) )
C           GH0HWZ(i)(i=1..3)  / Higgs_0(i) -Higgs_(-) -W(dag.)-Z
C           GH0HWA(i)(i=1..3)  / Higgs_0(i) -Higgs_(-) -W(dag.)-Z
C           GHHZZ              / Higgs_(+) -Higgs_(-) -Z-Z
C           GH02Z2             / Higgs_0(i) -Higgs_0(i) -Z-Z
C                                ( same for all H_0(i) )
C           GHHZA              / Higgs_(+) -Higgs_(-) -Z-A
C           GHHAA              / Higgs_(+) -Higgs_(-) -A-A
C----------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-H, O-Z)
      COMPLEX*16 GH0WW(2), GH0ZZ(2), GH0HCW(3), GH0PHZ(2),
     &           GH0HWZ(3), GH0HWA(3)
      COMPLEX*16 IMGI
      DATA IMGI  /(0.D0,1.D0)/
C----------------------------------------------------------------------C
      PI = ACOS(0.)*2.
      GWEAK = SQRT(ALPHE*4.*PI/SIN2W)
      EC    = SQRT(ALPHE*4.*PI)
 
      COS2W = 1.-SIN2W
      COSW  = SQRT(COS2W)
      ZM    = WM/COSW
      GZ    = GWEAK/COSW
 
C     ALPHA = ATAN2( ABS(TANA), SIGN(1.,TANA) )
C     ALPHA = DATAN2( DABS(DBLE(TANA)), DSIGN(1.D0,DBLE(TANA)) )
      TAN2A = TANA**2
      COS2A = 1./(1.+TAN2A)
      COSA  = SIGN( SQRT(COS2A), TANA )
      SIN2A = 1.-COS2A
      SINA  = SQRT(SIN2A)
C     BETA  = ATAN2( ABS(TANB), SIGN(1.,TANB) )
C     BETA  = DATAN2( DABS(DBLE(TANB)), DSIGN(1.D0,DBLE(TANB)) )
      TAN2B = TANB**2
      COS2B = 1./(1.+TAN2B)
      COSB  = SIGN( SQRT(COS2B), TANB )
      SIN2B = 1.-COS2B
      SINB  = SQRT(SIN2B)
 
C     DELT = ALPHA-BETA
C     COSD = COS(DELT)
C     SIND = SIN(DELT)
      COSD = COSA*COSB + SINA*SINB
C>>> 96/08/01  Bug fixup due to nojirim
C      SIND = SINA*COSB + COSA*SINB
      SIND = SINA*COSB - COSA*SINB
C>>>
 
      GH0WW(1) = -GWEAK*SIND *WM
      GH0WW(2) =  GWEAK*COSD *WM
 
      GH0ZZ(1) = -0.5*GZ*SIND *ZM
      GH0ZZ(2) =  0.5*GZ*COSD *ZM
 
      GH0HCW(1) = 0.5*GWEAK*COSD *IMGI
      GH0HCW(2) = 0.5*GWEAK*SIND *IMGI
      GH0HCW(3) =-0.5*GWEAK
 
      GH0PHZ(1) =-0.5*GZ*COSD
      GH0PHZ(2) =-0.5*GZ*SIND
 
      GHHZ      = -(-0.5+SIN2W)*GZ *IMGI
      GHHA      = EC *IMGI
 
      GHHWW     = 0.5*GWEAK**2
      GH02W2    = 0.25*GWEAK**2
      GH0HWZ(1) =-0.5*GWEAK*GZ*SIN2W*COSD
      GH0HWZ(2) =-0.5*GWEAK*GZ*SIN2W*SIND
      GH0HWZ(3) =-0.5*GWEAK*GZ*SIN2W *IMGI
      GH0HWA(1) = 0.5*GWEAK*EC*COSD
      GH0HWA(2) = 0.5*GWEAK*EC*SIND
      GH0HWA(3) = 0.5*GWEAK*EC *IMGI
      GHHZZ     = GZ**2*(-0.5+SIN2W)**2
      GH02Z2    = GZ**2/8.
      GHHZA     =-2.*GZ*EC*(-0.5+SIN2W)
      GHHAA     = EC**2
 
 
      RETURN
      END
C----------------------------------------------------------------------C
      SUBROUTINE  INOHGS( ALPHE, SIN2W, TANB, TANA,
     &                    PHIL, PHIR, ER, ON, ETA,
     &                    GCCH, GCNH, GNNH )
C----------------------------------------------------------------------C
C purpose : give couplings for Ino-Higgs couplings
C inputs  : ALPHE, SIN2W     / S.M. parameters
C           TANB             / vev. ratio
C           TANA             / tan of Higgs mix. angle
C           PHIL, PHIR, ER   / chargino mix. angles & phase
C           ON(4,4), ETA(4)  / neutralino mix. angles & phase
C outputs : GCCH(i,j,k,l) (i,j,k=1..2 ) C.ino(j)-C.ino(k)-H_0(l)
C                         (   l=1..3  )
C           GCNH(i,j,k) (i,j=1..2,  k=1..4 ) C.ino(j)-N.ino(k)-H_(-)
C           GNNH(i,j,k,l)(i=1..2, j,k=1..4 ) N.ino(j)-N.ino(k)-H_0(l)
C                         (   l=1..3  )
C              i=1/left-handed,  i=2/right-handed
C----------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-H, O-Z)
      REAL*8     ON(4,4)
      COMPLEX*16 ETA(4)
      COMPLEX*16 GSCCH(2,2), GCCH(2,2,2,3), GPCCP(2)
      COMPLEX*16 GCNH(2,2,4), GNNH(2,4,4,3), SS(3),PP(3), SS0, PP0
      COMPLEX*16 IMGI, CMPLX
      DATA IMGI  /(0.D0,1.D0)/
C----------------------------------------------------------------------C
      PI = ACOS(0.)*2.
      SQRT2 = SQRT(2.)
      GWEAK = SQRT(ALPHE*4.*PI/SIN2W)
      COS2W = 1.-SIN2W
      SINW  = SQRT(SIN2W)
      COSW  = SQRT(COS2W)
      TANW  = SINW/COSW
      GZ    = GWEAK/COSW
 
C  0 < beta < pi
C     BETA  = ATAN2( ABS(TANB), SIGN(1.,TANB) )
C     BETA  = DATAN2( DABS(DBLE(TANB)), DSIGN(1.D0,DBLE(TANB)) )
C     COSB  = COS(BETA)
C     SINB  = SIN(BETA)
      TAN2B = TANB**2
      COS2B = 1./(1.+TAN2B)
      COSB  = SIGN( SQRT(COS2B), TANB )
      SIN2B = 1.-COS2B
      SINB  = SQRT(SIN2B)
C     ALPHA = ATAN2( ABS(TANA), SIGN(1.,TANA) )
C     ALPHA = DATAN2( DABS(DBLE(TANA)), DSIGN(1.D0,DBLE(TANA)) )
C     SINA  = SIN(ALPHA)
C     COSA  = COS(ALPHA)
      TAN2A = TANA**2
      COS2A = 1./(1.+TAN2A)
      COSA  = SIGN( SQRT(COS2A), TANA )
      SIN2A = 1.-COS2A
      SINA  = SQRT(SIN2A)
 
      COSPL = COS(PHIL)
      SINPL = SIN(PHIL)
      COSPR = COS(PHIR)
      SINPR = SIN(PHIR)
 
C---
C Chargino-Chargino-Higgs
C---
C  GSCCH(i,j) : i=chargino(i), j=Higgs(j)
      GSCCH(1,1) = -SINA*SINPL*COSPR +COSA*COSPL*SINPR
      GSCCH(2,1) =  ER*(SINA*COSPL*SINPR -COSA*SINPL*COSPR)
      GSCCH(1,2) =  COSA*SINPL*COSPR +SINA*COSPL*SINPR
      GSCCH(2,2) =  ER*(-COSA*COSPL*SINPR -SINA*SINPL*COSPR)
      GCCH(1,1,1,1) = GSCCH(1,1)
      GCCH(2,1,1,1) = GSCCH(1,1)
      GCCH(1,2,2,1) = GSCCH(2,1)
      GCCH(2,2,2,1) = GSCCH(2,1)
      GCCH(1,1,1,2) = GSCCH(1,2)
      GCCH(2,1,1,2) = GSCCH(1,2)
      GCCH(1,2,2,2) = GSCCH(2,2)
      GCCH(2,2,2,2) = GSCCH(2,2)
 
      GPCCP(1) = SINB*SINPL*COSPR +COSB*COSPL*SINPR
      GPCCP(2) = ER*(-SINB*COSPL*SINPR -COSB*SINPL*COSPR)
      GCCH(1,1,1,3) =-GPCCP(1)
      GCCH(2,1,1,3) = GPCCP(1)
      GCCH(1,2,2,3) =-GPCCP(2)
      GCCH(2,2,2,3) = GPCCP(2)
 
      GCCH(1,1,2,1) = -SINA*COSPL*COSPR -COSA*SINPL*SINPR
      GCCH(2,1,2,1) = ER*(SINA*SINPL*SINPR +COSA*COSPL*COSPR)
      GCCH(1,2,1,1) = GCCH(2,1,2,1)
      GCCH(2,2,1,1) = GCCH(1,1,2,1)
      GCCH(1,1,2,2) = COSA*COSPL*COSPR -SINA*SINPL*SINPR
      GCCH(2,1,2,2) = ER*(-COSA*SINPL*SINPR +SINA*COSPL*COSPR)
      GCCH(1,2,1,2) = GCCH(2,1,2,2)
      GCCH(2,2,1,2) = GCCH(1,1,2,2)
 
      GCCH(1,1,2,3) = -SINB*COSPL*COSPR +COSB*SINPL*SINPR
      GCCH(2,1,2,3) = ER*(-SINB*SINPL*SINPR +COSB*COSPL*COSPR)
      GCCH(1,2,1,3) = -GCCH(2,1,2,3)
      GCCH(2,2,1,3) = -GCCH(1,1,2,3)
 
      DO 5 I=1,2
      DO 5 J=1,2
      DO 5 K=1,2
      DO 5 L=1,3
5     GCCH(I,J,K,L) = -GCCH(I,J,K,L) *GWEAK/SQRT2
 
      DO 6 I=1,2
      DO 6 J=1,2
      DO 6 K=1,2
6     GCCH(I,J,K,3) =  GCCH(I,J,K,3) *IMGI
 
C---
C Chargino-Neutralino-Higgs
C---
      DO 10 I=1,4
       ONI21T = (ON(I,2)+ON(I,1)*TANW)/SQRT2
       GCNH(1,1,I) = COSB*( COSPR*ON(I,4) +SINPR*ONI21T )*CONJG(ETA(I))
       GCNH(2,1,I) = SINB*( COSPL*ON(I,3) -SINPL*ONI21T ) * ETA(I)
       GCNH(1,2,I) = COSB*ER*(-SINPR*ON(I,4)
     &               +COSPR*ONI21T )*CONJG(ETA(I))
       GCNH(2,2,I) = SINB*(-SINPL*ON(I,3) -COSPL*ONI21T )*ETA(I)
10    CONTINUE
 
      DO 15 I=1,2
      DO 15 J=1,2
      DO 15 K=1,4
15    GCNH(I,J,K) = GCNH(I,J,K) * (-GWEAK)
 
C---
C Neutralino-Neutralino-Higgs
C---
      DO 30 I=1,4
       ONI21T = (ON(I,2)-ON(I,1)*TANW)
       SS0= ETA(I)**2* ONI21T*( -SINA*ON(I,3)-COSA*ON(I,4) )
       GNNH(1,I,I,1) = SS0
       GNNH(2,I,I,1) = SS0
       SS0= ETA(I)**2* ONI21T*( COSA*ON(I,3)-SINA*ON(I,4) )
       GNNH(1,I,I,2) = SS0
       GNNH(2,I,I,2) = SS0
       PP0=-ETA(I)**2* ONI21T*(-SINB*ON(I,3)+COSB*ON(I,4) )
       GNNH(1,I,I,3) =-PP0*IMGI
       GNNH(2,I,I,3) = PP0*IMGI
30    CONTINUE
 
      DO 40 I=1,4
      DO 40 J=1,4
       IF(I.EQ.J) GO TO 40
       CMPLX  = ETA(I)*ETA(J)
       ONI21T = (ON(I,2)-ON(I,1)*TANW)
       ONJ21T = (ON(J,2)-ON(J,1)*TANW)
       FF1=   ONI21T*(-SINA*ON(J,3)-COSA*ON(J,4))
     &       +ONJ21T*(-SINA*ON(I,3)-COSA*ON(I,4))
       FF2=   ONI21T*( COSA*ON(J,3)-SINA*ON(J,4))
     &       +ONJ21T*( COSA*ON(I,3)-SINA*ON(I,4))
       FF3=   ONI21T*(-SINB*ON(J,3)+COSB*ON(J,4))
     &       +ONJ21T*(-SINB*ON(I,3)+COSB*ON(I,4))
       IF( REAL(CMPLX) .NE.0 ) THEN
        SS(1) = REAL(CMPLX)*FF1
        PP(1) = 0.
        SS(2) = REAL(CMPLX)*FF2
        PP(2) = 0.
        SS(3) = 0.
        PP(3) =-REAL(CMPLX)*FF3
       ELSE
        SS(1) = 0.
        PP(1) = IMAG(CMPLX)*FF1
        SS(2) = 0.
        PP(2) = IMAG(CMPLX)*FF2
        SS(3) = IMAG(CMPLX)*FF3
        PP(3) = 0.
       ENDIF
       GNNH(1,I,J,1) = SS(1) -IMGI*PP(1)
       GNNH(2,I,J,1) = SS(1) +IMGI*PP(1)
       GNNH(1,I,J,2) = SS(2) -IMGI*PP(2)
       GNNH(2,I,J,2) = SS(2) +IMGI*PP(2)
       GNNH(1,I,J,3) = SS(3) -IMGI*PP(3)
       GNNH(2,I,J,3) = SS(3) +IMGI*PP(3)
40    CONTINUE
 
      DO 50 I=1,2
      DO 50 J=1,4
      DO 50 K=1,4
      DO 50 L=1,3
50    GNNH(I,J,K,L) = GNNH(I,J,K,L)*(-GWEAK)/2.
 
 
      RETURN
      END
