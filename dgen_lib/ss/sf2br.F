C----------------------------------------------------------------------C
      SUBROUTINE SF2BR(  WM, SIN2W, ALPHE, ALPHS,
     &                   HM1, G3MAS, TANB,
     &                   EMAS, UMAS, DMAS,
     &                   SFMSL, SFMSR, SNEUM, PSMAS,
     &                   SWM, SZM, HGCMAS, HG0MAS,
     &                   DCWSEL, DCWSER, DCWSNL,
     &                   BRSEL, BRSER, BRSNL,
     &                   DCWIN1, BRINO1 )
C----------------------------------------------------------------------C
C Purpose: calculate decay widthes & brancing ratios for
C          2-body-sferiom decays such as
C            1.) charged slepton --> neutrino + chargino
C            2.) charged slepton --> charged lepton + neutralino
C            3.) sneutrino       --> charged lepton + chargino
C            4.) sneutrino       --> neutrino + neutralino
C Inputs:  WM, SIN2W, ALPHE, ALPHS / parms. for SM.
C          HM1, G3MAS, TANB        / parms. for MSSM.
C          EMAS, UMAS, DMAS        / masses for charged lepton,
C                                    u-type & d-type quarks
C          SFMSL, SFMSR, SNEUM     / masses for sleptons
C          PSMAS                   / pseudo-scaler higgs mass
C Output: SWM, SZM                 / chargino,neutralino masses
C         DCWSEL, DCWSER, DCWSNL   / decay widthes
C         BRSEL, BRSER, BRSNL      / branching ratios
C         DCWIN1, BRINO1           / 2-body decay widthes & branching
C                                    ratios of chargino(1)-->neutralino
C----------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-H, O-Z)
      REAL*8     OL(2,2), OR(2,2), ER
      REAL*8     ON(4,4)
      COMPLEX*16 ETA(4), UN(4,4), GL, GR
      REAL*8     SWM(2), SZM(4), INOMS
 
      COMPLEX*16 GCESNL(2,4), GCNSEL(2,4), GCNSER(2,4),
     &           GCUSDL(2,4), GCUSDR(2,4), GCDSUL(2,4), GCDSUR(2,4)
      COMPLEX*16 GNNSNL(2,4), GNESEL(2,4), GNESER(2,4),
     &           GNUSUL(2,4), GNUSUR(2,4), GNDSDL(2,4), GNDSDR(2,4)
 
      COMPLEX*16 GNCW(2,4,2)
      COMPLEX*16 GCCZ(2,2,2)
      COMPLEX*16 GNNZ(2,4,4)
 
      COMPLEX*16 GCCH(2,2,2,3)
      COMPLEX*16 GCNH(2,2,4), GNNH(2,4,4,3)
 
      REAL*8  DCWSEL(6), DCWSER(6), DCWSNL(6)
      REAL*8  BRSEL(6), BRSER(6), BRSNL(6)
 
      REAL*8 DCWIN1(8), BRINO1(8)
 
      REAL*8 HG0MAS(2)
C------------------------------------------------
 
C    No printout diagonalization results
      IFLG = 1
      CALL       INOMIX( ZM, WM, SIN2W, ALPHE, ALPHS,
     &                   HM1, G3MAS, TANB,
     &                   SWM, SZM,
     &                   OL, OR, ER, PHIL, PHIR,
     &                   ON, UN, ETA,
     &                   IFLG  )
 
C   calculate higgs mixing
      CALL  HGSMIX( SIN2W, WM, TANB, PSMAS, HG0MAS, HGCMAS, TANA)
 
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
 
      DO 5 I=1,8
        DCWIN1(I) = 0
        BRINO1(I) = 0
5     CONTINUE
 
      DO 11 I=1,6
        DCWSEL(I) = 0
        DCWSER(I) = 0
        DCWSNL(I) = 0
        BRSEL(I) = 0
        BRSER(I) = 0
        BRSNL(I) = 0
11    CONTINUE
 
C
C---
C  calculate decay widthes for 2-body-sfermion decay
C---
C left-handed Selectron decay
C   Se_L --> nu_e + chargino(i)
      FMAS  = 0.
      DO 10 I=1,2
        INOMS = SWM(I)
C      note the convention
        GR = GCNSEL(1,I)
        GL = GCNSEL(2,I)
        CALL SCF1F2( SFMSL, FMAS, INOMS, GL, GR, DCYWDT )
        DCWSEL(I) = DCYWDT
10    CONTINUE
C   Se_L --> e + neutralino(i)
      FMAS  = EMAS
      DO 20 I=1,4
        INOMS = SZM(I)
C      note the convention
        GR = GNESEL(1,I)
        GL = GNESEL(2,I)
        CALL SCF1F2( SFMSL, FMAS, INOMS, GL, GR, DCYWDT )
        DCWSEL(I+2) = DCYWDT
20    CONTINUE
 
C Right-handed Selectron decay
C   Se_R --> nu_e + chargino(i)
      FMAS  = 0.
      DO 30 I=1,2
        INOMS = SWM(I)
C      note the convention
        GR = GCNSER(1,I)
        GL = GCNSER(2,I)
        CALL SCF1F2( SFMSR, FMAS, INOMS, GL, GR, DCYWDT )
        DCWSER(I) = DCYWDT
30    CONTINUE
C   Se_R --> e + neutralino(i)
      FMAS  = EMAS
      DO 40 I=1,4
        INOMS = SZM(I)
C      note the convention
        GR = GNESER(1,I)
        GL = GNESER(2,I)
        CALL SCF1F2( SFMSR, FMAS, INOMS, GL, GR, DCYWDT )
        DCWSER(I+2) = DCYWDT
40    CONTINUE
 
C Sneutrino decay
C  snu_e --> electron + chargino(i)
      FMAS  = EMAS
      DO 60 I=1,2
        INOMS = SWM(I)
C      note the convention
        GR = GCESNL(1,I)
        GL = GCESNL(2,I)
        CALL SCF1F2( SNEUM, FMAS, INOMS, GL, GR, DCYWDT )
        DCWSNL(I) = DCYWDT
60    CONTINUE
C  snu_e --> neutrino + neutralino(i)
      FMAS  = 0.
      DO 65 I=1,4
        INOMS = SZM(I)
C      note the convention
        GR = GNNSNL(1,I)
        GL = GNNSNL(2,I)
        CALL SCF1F2( SNEUM, FMAS, INOMS, GL, GR, DCYWDT )
        DCWSNL(I+2) = DCYWDT
65    CONTINUE
 
C---
C  calculate branching ratios
C---
 
      TOTDL = 0
      TOTDR = 0
      TOTDN = 0
      DO 50 I=1,6
       TOTDL = TOTDL+DCWSEL(I)
       TOTDR = TOTDR+DCWSER(I)
       TOTDN = TOTDN+DCWSNL(I)
50    CONTINUE
 
      IF( TOTDL.GT.0 ) THEN
        DO 55 I=1,6
         BRSEL(I) = DCWSEL(I)/TOTDL
55      CONTINUE
      ENDIF
 
      IF( TOTDR.GT.0 ) THEN
        DO 56 I=1,6
         BRSER(I) = DCWSER(I)/TOTDR
56      CONTINUE
      ENDIF
 
      IF( TOTDN.GT.0 ) THEN
        DO 57 I=1,6
         BRSNL(I) = DCWSNL(I)/TOTDN
57      CONTINUE
      ENDIF
 
C--- calculate chargino_1 2-body decay widthes ( to neutralino )
C  chargino_1 --> neutralino_j + W(-)
      F1MAS = SWM(1)
      DO 70 J = 1,4
        F2MAS = SZM(J)
        VMAS = WM
        GL = GNCW(1,J,1)
        GR = GNCW(2,J,1)
        CALL F1F2VD( F1MAS, F2MAS, VMAS, GL, GR, DCYWDT )
        DCWIN1(J) = DCYWDT
70    CONTINUE
 
C  chargino_1 --> neutralino_j + H(-)
      F1MAS = SWM(1)
      DO 72 J = 1,4
        F2MAS = SZM(J)
        SCMAS = HGCMAS
        GL = CONJG(GCNH(2,2,J))
        GR = CONJG(GCNH(1,2,J))
        CALL F1F2SC( F1MAS, F2MAS, SCMAS, GL, GR, DCYWDT )
        DCWIN1(J+4) = DCYWDT
72    CONTINUE
 
 
      TOTD = 0
      DO 75 I=1,8
      IF( DCWIN1(I).LT.0 ) THEN
       WRITE(6,*) ' I=',I,' G3M=',G3MAS
      ENDIF
75    TOTD = TOTD+DCWIN1(I)
      IF( TOTD.GT.0 ) THEN
        DO 80 I=1,8
80      BRINO1(I) = DCWIN1(I)/TOTD
      ENDIF
 
 
 
      RETURN
      END
