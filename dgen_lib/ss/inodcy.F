C calculate chargino/neutralino decay widthes/BR
C based on SUSY.FORT(TESTBR23,TESTBR32)
C----------------------------------------------------------------------C
      SUBROUTINE INODCY( WM, SIN2W, ALPHE, ALPHS,
     &                   HM1, G3MAS, TANB,
     &                   EMAS, UMAS, DMAS,
     &                   SFMSL, SFMSR, SNEUM, PSHMAS,
     &                   SWM, SZM, HGCMAS, HG0MAS,
     &                   DCCINO, BRCINO, DCNINO, BRNINO )
C----------------------------------------------------------------------C
C Purpose: calculate decay widthes & brancing ratios for
C          chargino/neutralino decays
C Inputs:  WM, SIN2W, ALPHE, ALPHS / parms. for SM.
C          HM1, G3MAS, TANB        / parms. for MSSM.
C          EMAS, UMAS, DMAS        / masses for charged lepton,
C                                    u-type & d-type quarks
C          SFMSL, SFMSR, SNEUM     / masses for sleptons
C          PSHMAS                  / pseudo-scaler higgs mass
C Output: SWM, SZM                 / chargino,neutralino masses
C         DCCINO, DCNINO
C         BRCINO, BRNINO
C----------------------------------------------------------------------C
      IMPLICIT REAL*8 ( A-H, O-Z )
      REAL*8  OL(2,2), OR(2,2), ER
      REAL*8  ON(4,4)
      COMPLEX*16 ETA(4), UN(4,4), GL, GR
      COMPLEX*16 G1L,G1R, G2L,G2R
      REAL*8  SWM(2), SZM(4)
      REAL*8  DCCINO(40,4), DCNINO(40,2)
      REAL*8  BRCINO(40,4), BRNINO(40,2)
 
      COMPLEX*16 GCESNL(2,4), GCNSEL(2,4), GCNSER(2,4),
     &          GCUSDL(2,4), GCUSDR(2,4), GCDSUL(2,4), GCDSUR(2,4)
      COMPLEX*16 GNNSNL(2,4), GNESEL(2,4), GNESER(2,4),
     &          GNUSUL(2,4), GNUSUR(2,4), GNDSDL(2,4), GNDSDR(2,4)
 
      COMPLEX*16 GNCW(2,4,2)
      COMPLEX*16 GCCZ(2,2,2)
      COMPLEX*16 GNNZ(2,4,4)
 
      REAL*8 HG0MAS(2)
      COMPLEX*16 GCCH(2,2,2,3)
      COMPLEX*16 GCNH(2,2,4), GNNH(2,4,4,3)
C------------------------------------------------
C    No printout diagonalization results
      IFLG = 1
      CALL       INOMIX( ZM, WM, SIN2W, ALPHE, ALPHS,
     &                   HM1, G3MAS, TANB,
     &                   SWM, SZM,
     &                   OL, OR, ER, PHIL, PHIR,
     &                   ON, UN, ETA,
     &                   IFLG  )
 
C     require that chargino(1) is heavier than neutralino(1)
      IF( SWM(1).LT.SZM(1) ) RETURN
 
      DO 1 I=1,40
       DO 2 J=1,4
        DCNINO(I,J)=0.D0
        BRNINO(I,J)=0.D0
2      CONTINUE
       DO 3 J=1,2
        DCCINO(I,J)=0.D0
        BRCINO(I,J)=0.D0
3      CONTINUE
1     CONTINUE
 
C   calculate higgs mixing
      CALL  HGSMIX( SIN2W, WM, TANB, PSHMAS, HG0MAS, HGCMAS, TANA)
 
C  determine Sfermion-fermion-chargino/neutralino couplings
      CALL        INOFSF( ALPHE, SIN2W, WM,
     &                    EMAS, UMAS, DMAS,
     &                    TANB, PHIL, PHIR, ER, ON, ETA,
     &                    GCESNL, GCNSEL, GCNSER,
     &                    GCUSDL, GCUSDR, GCDSUL, GCDSUR,
     &                    GNNSNL, GNESEL, GNESER,
     &                    GNUSUL, GNUSUR, GNDSDL, GNDSDR )
 
C       determine ino-ino-gauge_boson  couplings
      CALL        INOGBS( ALPHE, SIN2W,
     &                    PHIL, PHIR, ER, ON, ETA,
     &                    GNCW, GCCZ, GNNZ )
 
C       determine ino-ino-Higgs couplings
      CALL       INOHGS( ALPHE, SIN2W, TANB, TANA,
     &                    PHIL, PHIR, ER, ON, ETA,
     &                    GCCH, GCNH, GNNH )
 
C       set SM couplings etc.
      GAMW = 2.25D0
      GAMZ = 2.53D0
      PI  = ACOS(0.D0)*2.D0
      EC2 = 4.D0*PI*ALPHE
      GWEAK = SQRT(EC2/SIN2W)
      COS2W = 1.D0-SIN2W
      GZ    = SQRT(EC2/COS2W/SIN2W)
C   f-f_bar-Z
      GZNL = -GZ*( 0.5D0 )
      GZNR = 0.D0
      GZEL = -GZ*(-0.5D0+SIN2W)
      GZER = -GZ*(     SIN2W)
      GZUL = -GZ*( 0.5D0-2.D0/3.D0*SIN2W)
      GZUR = -GZ*(    -2.D0/3.D0*SIN2W)
      GZDL = -GZ*(-0.5D0+1.D0/3.D0*SIN2W)
      GZDR = -GZ*(     1.D0/3.D0*SIN2W)
C   f1-f2_bar-W
      GWFL = -GWEAK/SQRT(2.D0)
      GWFR = 0.D0
 
C=============================
C neutralino decay
C=============================
C loop over neutralino_J
      DO 1000 J=1,4
 
C----
C 2-bodies decays
C----
C      calculate decay widthes of neutralino_4 2-body decays
C  neutralino_J --> chargino_i(-) + W(+)  + h.c.
      N = 0
      DO 10 I=1,2
      F1MAS = SZM(J)
      F2MAS = SWM(I)
      VMAS  = WM
      GL = GNCW(1,J,I)
      GR = GNCW(2,J,I)
      CALL F1F2VD( F1MAS, F2MAS, VMAS, GL, GR, DCYWDT )
      N = N+1
      DCNINO(N,J) = DCYWDT *2.D0
10    CONTINUE
C N is 2
 
C  neutralino_J --> chargino_i(-) + Higgs(+)  + h.c.
      DO 11 I=1,2
      F1MAS = SZM(J)
      F2MAS = SWM(I)
      SCMAS = HGCMAS
      GL = GCNH(1,I,J)
      GR = GCNH(2,I,J)
      CALL F1F2SC( F1MAS, F2MAS, SCMAS, GL, GR, DCYWDT )
      N = N+1
      DCNINO(N,J) = DCYWDT *2.D0
11    CONTINUE
C N is 4
 
C  neutralino_J --> neutralino_i + Higgs0_k
      DO 12 I=1,4
      DO 12 K=1,3
      F1MAS = SZM(J)
      F2MAS = SZM(I)
      IF( K.LE.2 ) THEN
        SCMAS = HG0MAS(K)
      ELSE
        SCMAS = PSHMAS
      ENDIF
      GL = GNNH(1,I,J,K)
      GR = GNNH(2,I,J,K)
      CALL F1F2SC( F1MAS, F2MAS, SCMAS, GL, GR, DCYWDT )
      N = N+1
      DCNINO(N,J) = DCYWDT
12    CONTINUE
C N is 16
 
C  neutralino_J --> neutralino_i + Z0
      DO 13 I=1,4
      F1MAS = SZM(J)
      F2MAS = SZM(I)
      VMAS = ZM
      GL = GNNZ(1,I,J)
      GR = GNNZ(2,I,J)
      CALL F1F2VD( F1MAS, F2MAS, VMAS, GL, GR, DCYWDT )
      N = N+1
      DCNINO(N,J) = DCYWDT
13    CONTINUE
C N is 20
 
C----
C 3-bodies decays
C----
C  neutralino_J --> neutralino_I + f f_bar
C    by virtual Z exchange, mass(f) is neglected
      AM1 = SZM(J)
      AMV = ZM
      GAMV= GAMZ
 
      DO 14 I=1,4
       AM2 = SZM(I)
       G1L = GNNZ(1,I,J)
       G1R = GNNZ(2,I,J)
C   1. (f f_bar) = (nu nu_bar) type
       G2L = GZNL
       G2R = GZNR
       NF = 3
       NCOL=1
       CALL FR3DCY( G1L,G1R, G2L,G2R, AM1,AM2, AMV,GAMV,  DCYWDT )
       N = N+1
       DCNINO(N,J) = DCYWDT*NF*NCOL
C   2. (f f_bar) = (e e_bar)   type
       G2L = GZEL
       G2R = GZER
       NF = 3
       NCOL=1
       CALL FR3DCY( G1L,G1R, G2L,G2R, AM1,AM2, AMV,GAMV,  DCYWDT )
       N = N+1
       DCNINO(N,J) = DCYWDT*NF*NCOL
C   3. (f f_bar) = (u u_bar)   type
       G2L = GZUL
       G2R = GZUR
       NF = 2
       NCOL=3
       CALL FR3DCY( G1L,G1R, G2L,G2R, AM1,AM2, AMV,GAMV,  DCYWDT )
       N = N+1
       DCNINO(N,J) = DCYWDT*NF*NCOL
C   4. (f f_bar) = (d d_bar)   type
       G2L = GZDL
       G2R = GZDR
       NF = 3
       NCOL=3
       CALL FR3DCY( G1L,G1R, G2L,G2R, AM1,AM2, AMV,GAMV,  DCYWDT )
       N = N+1
       DCNINO(N,J) = DCYWDT*NF*NCOL
14    CONTINUE
C N is 36
 
C  neutralino_J --> chargino_I + f1 f2_bar
C     by virtual W exchange, mass(f) is neglected
      AM1 = SZM(J)
      AMV = WM
      GAMV= GAMW
      G2L = GWFL
      G2R = GWFR
      DO 15 I=1,2
       AM2 = SWM(I)
C       be careful to give couplings
       G1L = AIMAG(GNCW(1,J,I))
       G1R = AIMAG(GNCW(2,J,I))
C   1.(f1 f2_bar) = (e nu_bar)  type
       NF = 3
       NCOL=1
       CALL FR3DCY( G1L,G1R, G2L,G2R, AM1,AM2, AMV,GAMV,  DCYWDT )
       N = N+1
       DCNINO(N,J) = DCYWDT*NF*NCOL
C   2.(f1 f2_bar) = (u d_bar)   type
       NF = 2
       NCOL=3
       CALL FR3DCY( G1L,G1R, G2L,G2R, AM1,AM2, AMV,GAMV,  DCYWDT )
       N = N+1
       DCNINO(N,J) = DCYWDT*NF*NCOL
15    CONTINUE
C N is 40
 
1000  CONTINUE
C end of neutralino loop
 
C=============================
C chargino decay
C=============================
C loop over chargino_K
      DO 1100 K=1,2
C----
C 2-bodies decays
C----
C      calculate decay widthes of chargino 2-body decays
C  chargino_K --> chargino_1 + Z0
      N=0
      IF( K.EQ.2 ) THEN
        F1MAS = SWM(2)
        F2MAS = SWM(1)
        VMAS = ZM
        GL = GCCZ(1,1,2)
        GR = GCCZ(2,1,2)
        CALL F1F2VD( F1MAS, F2MAS, VMAS, GL, GR, DCYWDT )
        N = N+1
        DCCINO(K,N) = DCYWDT
      ELSE
        N = N+1
        DCCINO(K,N) = 0
      ENDIF
C N is 1
 
C  chargino_K --> chargino_1 + H0_j
      IF( K.EQ.2 ) THEN
      F1MAS = SWM(2)
      F2MAS = SWM(1)
      DO 20 J = 1,3
        IF( J.LE.2 ) THEN
          SCMAS = HG0MAS(J)
        ELSE
          SCMAS = PSHMAS
        ENDIF
        GL = GCCH(1,2,1,J)
        GR = GCCH(2,2,1,J)
        CALL F1F2SC( F1MAS, F2MAS, SCMAS, GL, GR, DCYWDT )
        N = N+1
        DCCINO(K,N) = DCYWDT
20    CONTINUE
      ELSE
      DO 21 J=1,3
        N = N+1
        DCCINO(K,N) = 0
21    CONTINUE
      ENDIF
C N is 4
 
C  chargino_K --> neutralino_j + W(-)
      F1MAS = SWM(K)
      DO 22 J = 1,4
        F2MAS = SZM(J)
        VMAS = WM
        GL = GNCW(1,J,K)
        GR = GNCW(2,J,K)
        CALL F1F2VD( F1MAS, F2MAS, VMAS, GL, GR, DCYWDT )
        N = N+1
        DCCINO(K,N) = DCYWDT
22    CONTINUE
C N is 8
 
C  chargino_K --> neutralino_j + H(-)
      F1MAS = SWM(K)
      DO 23 J = 1,4
        F2MAS = SZM(J)
        SCMAS = HGCMAS
        GL = CONJG(GCNH(2,K,J))
        GR = CONJG(GCNH(1,K,J))
        CALL F1F2SC( F1MAS, F2MAS, SCMAS, GL, GR, DCYWDT )
        N = N+1
        DCCINO(K,N) = DCYWDT
23    CONTINUE
C N is 12
 
C----
C 3-bodies decays
C----
C  chargino_2 --> chargino_1 + f f_bar
C    by virtual Z exchange, mass(f) is neglected
      IF( K.EQ.2 ) THEN
       AM1 = SWM(2)
       AM2 = SWM(1)
       AMV = ZM
       GAMV= GAMZ
       G1L = GCCZ(1,1,2)
       G1R = GCCZ(2,1,2)
C   1. (f f_bar) = (nu nu_bar) type
       G2L = GZNL
       G2R = GZNR
       NF = 3
       NCOL=1
       CALL FR3DCY( G1L,G1R, G2L,G2R, AM1,AM2, AMV,GAMV,  DCYWDT )
       N = N+1
       DCCINO(K,N) = DCYWDT*NF*NCOL
C   2. (f f_bar) = (e e_bar)   type
       G2L = GZEL
       G2R = GZER
       NF = 3
       NCOL=1
       CALL FR3DCY( G1L,G1R, G2L,G2R, AM1,AM2, AMV,GAMV,  DCYWDT )
       N = N+1
       DCCINO(K,N) = DCYWDT*NF*NCOL
C   3. (f f_bar) = (u u_bar)   type
       G2L = GZUL
       G2R = GZUR
       NF = 2
       NCOL=3
       CALL FR3DCY( G1L,G1R, G2L,G2R, AM1,AM2, AMV,GAMV,  DCYWDT )
       N = N+1
       DCCINO(K,N) = DCYWDT*NF*NCOL
C   4. (f f_bar) = (d d_bar)   type
       G2L = GZDL
       G2R = GZDR
       NF = 3
       NCOL=3
       CALL FR3DCY( G1L,G1R, G2L,G2R, AM1,AM2, AMV,GAMV,  DCYWDT )
       N = N+1
       DCCINO(K,N) = DCYWDT*NF*NCOL
      ELSE
       DO 24 J=1,4
         N = N+1
         DCCINO(K,N) = 0
24     CONTINUE
      ENDIF
C N is 16
 
C  chargino_K --> neutralino_j + f1 f2_bar
C     by virtual W exchange, mass(f) is neglected
      AM1 = SWM(K)
      AMV = WM
      GAMV= GAMW
      G2L = GWFL
      G2R = GWFR
      DO 25 J=1,4
       AM2 = SZM(J)
       G1L = GNCW(1,J,K)
       G1R = GNCW(2,J,K)
C    (f1 f2_bar) = (e nu_bar)  type
       NF = 3
       NCOL=1
       CALL FR3DCY( G1L,G1R, G2L,G2R, AM1,AM2, AMV,GAMV,  DCYWDT )
       N = N+1
       DCCINO(K,N) = DCYWDT*NF*NCOL
C    (f1 f2_bar) = (u d_bar)   type
       NF = 2
       NCOL=3
       CALL FR3DCY( G1L,G1R, G2L,G2R, AM1,AM2, AMV,GAMV,  DCYWDT )
       L = L+1
       N = N+1
       DCCINO(K,N) = DCYWDT*NF*NCOL
25    CONTINUE
C N is 24
 
1100  CONTINUE
C end of chargino loop
 
 
 
 
      RETURN
      END
