C----------------------------------------------------------------------C
      SUBROUTINE SF2BR(  ALP, S2W, WM, ZM, SWM, SZM, SFM, SFM3M, SPHI, 
     &   SHM, FM3P,
     &   GNCW, GCCZ, GNNZ, GCNH,GNNH,
     &   GCESNL, GCNSEL, GCNSER,
     &   GNNSNL, GNESEL, GNESER,
     &   GCSNL3,GCNSE1,GCNSE2,
     &   GNSNL3,GNESE1,GNESE2,
     &   GMSF,GMSW,GMSZ,GMSF3)

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
C         DCWIN2, BRINO2           / 2-body decay widthes & branching
C                                    ratios of neutralino(1)-->neutralino
C----------------------------------------------------------------------C
      IMPLICIT REAL*8 ( A-H, O-Z )
      REAL*8  SWM(2), SZM(4), SFM(7), SHM(4), SFM3M(7), SPHI(3)
      REAL*8  FM3P(3),PI
      REAL*8   GMSF(7),GMSW(2),GMSZ(4),GMSF3(7)
      COMPLEX*16 GCESNL(2,2), GCNSEL(2,2), GCNSER(2,2) 
      COMPLEX*16 GNNSNL(2,4), GNESEL(2,4), GNESER(2,4)  
      COMPLEX*16 GCSNL3(2,2), GCNSE1(2,2), GCNSE2(2,2)
      COMPLEX*16 GNSNL3(2,4), GNESE1(2,4), GNESE2(2,4)
      COMPLEX*16 GNCW(2,4,2), GCCZ(2,2,2), GNNZ(2,4,4)
      COMPLEX*16 GCNH(2,2,4), GNNH(2,4,4,3) 
      COMPLEX*16 GL, GR
      REAL*8  INOMS
      REAL*8  DCWSEL(6), DCWSER(6), DCWSNL(6)
      REAL*8  BRSEL(6), BRSER(6), BRSNL(6)
      REAL*8  DCWSE1(6), DCWSE2(8), DCWSN3(6)
      REAL*8  BRSE1(6), BRSE2(8), BRSN3(6)
      REAL*8 DCWIN1(8), BRINO1(8)
      REAL*8 DCWIN2(4), BRINO2(4)
C------------------------------------------------
C     WEAK COUPLING
      PI = ACOS(0.D0)*2.D0
      GWEAK2 = ALP/S2W*4.D0*PI
      GWEAK =SQRT(GWEAK2)
      GZ2=GWEAK2/(1.D0-S2W)
      GZ=SQRT(GZ2)
C
      DO 5 I=1,8
        DCWIN1(I) = 0
        BRINO1(I) = 0
5     CONTINUE
      DO 6 I=1,4
        DCWIN2(I) = 0
        BRINO2(I) = 0
6     CONTINUE
      DO 11 I=1,6
        DCWSEL(I) = 0
        DCWSER(I) = 0
        DCWSNL(I) = 0
        BRSEL(I) = 0
        BRSER(I) = 0
        BRSNL(I) = 0
        DCWSE1(I) = 0
        DCWSN3(I) = 0
        BRSE1(I) = 0
        BRSN3(I) = 0
11    CONTINUE
      DO 12 I=1,8
        DCWSE2(I) = 0 
        BRSE2(I) = 0
12    CONTINUE

C---
C  calculate decay widthes for 2-body-sfermion decay
C---
C left-handed Selectron decay
C   Se_L --> nu_e + chargino(i)
      FMAS  = 0.D0
      DO 10 I=1,2
        INOMS = SWM(I)
C      note the convention
        GR = GCNSEL(1,I)
        GL = GCNSEL(2,I)
        CALL SCF1F2( SFM(2), FMAS, INOMS, GL, GR, DCYWDT )
C       write(*,*) gl,gr
        DCWSEL(I) = DCYWDT
10    CONTINUE
C   Se_L --> e + neutralino(i)
      FMAS  = EMAS
      DO 20 I=1,4
        INOMS = SZM(I)
C      note the convention
        GR = GNESEL(1,I)
        GL = GNESEL(2,I)
        CALL SCF1F2( SFM(2), FMAS, INOMS, GL, GR, DCYWDT )
C        write(*,*) gl,gr
        DCWSEL(I+2) = DCYWDT
20    CONTINUE
 
C Right-handed Selectron decay
C   Se_R --> nu_e + chargino(i)
      FMAS  = 0.D0
      DO 30 I=1,2
        INOMS = SWM(I)
C      note the convention
        GR = GCNSER(1,I)
        GL = GCNSER(2,I)
        write(*,*) gr, gl
        CALL SCF1F2( SFM(3), FMAS, INOMS, GL, GR, DCYWDT )
        DCWSER(I) = DCYWDT
30    CONTINUE
C   Se_R --> e + neutralino(i)
      FMAS  = 0
      DO 40 I=1,4
        INOMS = SZM(I)
C      note the convention
        GR = GNESER(1,I)
        GL = GNESER(2,I)
        write(*,*) gr,gl
        CALL SCF1F2( SFM(3), FMAS, INOMS, GL, GR, DCYWDT )
        DCWSER(I+2) = DCYWDT
40    CONTINUE
 
C Sneutrino decay
C  snu_e --> electron + chargino(i)
      FMAS  = 0
      DO 60 I=1,2
        INOMS = SWM(I)
C      note the convention
        GR = GCESNL(1,I)
        GL = GCESNL(2,I)
        CALL SCF1F2( SFM(1), FMAS, INOMS, GL, GR, DCYWDT )
        DCWSNL(I) = DCYWDT
60    CONTINUE
C  snu_e --> neutrino + neutralino(i)
      FMAS  = 0.D0
      DO 65 I=1,4
        INOMS = SZM(I)
C      note the convention
        GR = GNNSNL(1,I)
        GL = GNNSNL(2,I)
        CALL SCF1F2( SFM(1), FMAS, INOMS, GL, GR, DCYWDT )
        DCWSNL(I+2) = DCYWDT
65    CONTINUE
 
C
C---
C  calculate decay widthes for 2-body- 3rd generation sfermion decay
C---
C left-handed Selectron decay
C   STAU_1 --> neu_e + chargino(i)
      FMAS  = 0.D0
      DO 21 I=1,2
        INOMS = SWM(I)
C      note the convention
        GR = GCNSE1(1,I)
        GL = GCNSE1(2,I)
        write(*,*) gr,gl
        CALL SCF1F2( SFM3M(2), FMAS, INOMS, GL, GR, DCYWDT )
        DCWSE1(I) = DCYWDT
 21   CONTINUE
C   STAU_1 --> TAU + neutralino(i)
      FMAS  = EMAS
      DO 22 I=1,4
        INOMS = SZM(I)
C      note the convention
        GR = GNESE1(1,I)
        GL = GNESE1(2,I)
        write(*,*) gr,gl
        CALL SCF1F2( SFM3M(2), FMAS, INOMS, GL, GR, DCYWDT )
        DCWSE1(I+2) = DCYWDT
 22   CONTINUE
 
C Right-handed Selectron decay
C   STAU_2 --> NU_TAU + chargino(i)
      FMAS  = 0.D0
      DO 31 I=1,2
        INOMS = SWM(I)
C      note the convention
        GR = GCNSE2(1,I)
        GL = GCNSE2(2,I)
C        write(*,*) gr,gl
        CALL SCF1F2( SFM3M(3), FMAS, INOMS, GL, GR, DCYWDT )
        DCWSE2(I) = DCYWDT
 31   CONTINUE
C   Se_2 --> e + neutralino(i)
      FMAS  = FM3P(1)
      DO 41 I=1,4
        INOMS = SZM(I)
C      note the convention
        GR = GNESE2(1,I)
        GL = GNESE2(2,I)
C        write(*,*) gr,gl
        CALL SCF1F2( SFM3M(3), FMAS, INOMS, GL, GR, DCYWDT )
        DCWSE2(I+2) = DCYWDT
 41     CONTINUE
C   STAU_2--> STAU_1 Z
       CNPHI=COS(SPHI(1))
       SNPHI=SIN(SPHI(1))
       write(*,*) cnphi, snphi, 'cnphi, snphi'   
       GR = GZ*CNPHI*SNPHI*0.5D0
C       write(*,*) gr
       CALL SC1SC2V(SFM3M(3),SFM3M(2),ZM,GR,DCYDWT)
       DCWSE2(7) =DCYDWT
C  STAU_2 ---> SNUTAU_1,W        
       GR=GWEAK*SNPHI/SQRT(2.D0)
C       write(*,*) gr
       CALL SC1SC2V(SFM3M(3),SFM3M(1),WM,GR,DCYDWT)
       DCWSE2(7) =DCYDWT
       
C Sneutrino decay
C  snu_tau --> TAU + chargino(i)
      FMAS  = FM3P(1)
      DO 61 I=1,2
        INOMS = SWM(I)
C      note the convention
        GR = GCSNL3(1,I)
        GL = GCSNL3(2,I)
        CALL SCF1F2( SFM3M(1), FMAS, INOMS, GL, GR, DCYWDT )
        DCWSN3(I) = DCYWDT
 61     CONTINUE
C  snu_tau --> neutrino + neutralino(i)
      FMAS  = 0.D0
      DO 66 I=1,4
        INOMS = SZM(I)
C      note the convention
        GR = GNSNL3(1,I)
        GL = GNSNL3(2,I)
        CALL SCF1F2( SFM3M(1), FMAS, INOMS, GL, GR, DCYWDT )
        DCWSN3(I+2) = DCYWDT
 66   ENDDO
 
C---
C  calculate branching ratios
C---
 
      TOTDL = 0
      TOTDR = 0
      TOTDN = 0
      TOTD1 = 0
      TOTD2 = 0
      TOTD3 = 0
      DO 50 I=1,6
       TOTDL = TOTDL+DCWSEL(I)
       TOTDR = TOTDR+DCWSER(I)
       TOTDN = TOTDN+DCWSNL(I)
       TOTD1 = TOTD1+DCWSE1(I)
       TOTD3 = TOTD3+DCWSN3(I)
50    CONTINUE
      write(*,*) totdl, 'totdl'
      DO 51 I=1,8 
       TOTD2 = TOTD2+DCWSE2(I)
51    CONTINUE
             IF( TOTDL.GT.0 ) THEN
        DO 54 I=1,6
         BRSEL(I) = DCWSEL(I)/TOTDL
 54   CONTINUE
      ENDIF
 
      IF( TOTDR.GT.0 ) THEN
        DO 55 I=1,6
         BRSER(I) = DCWSER(I)/TOTDR
 55   CONTINUE
      ENDIF
 
      IF( TOTDN.GT.0 ) THEN
        DO 56 I=1,6
         BRSNL(I) = DCWSNL(I)/TOTDN
 56   CONTINUE
      ENDIF

      IF( TOTD1.GT.0 ) THEN
        DO 57 I=1,6
         BRSE1(I) = DCWSE1(I)/TOTD1
 57   CONTINUE
      ENDIF
 
      IF( TOTD2.GT.0 ) THEN
        DO 58 I=1,6
         BRSE2(I) = DCWSE2(I)/TOTD2
 58   CONTINUE
      ENDIF
 
      IF( TOTD3.GT.0 ) THEN
        DO 59 I=1,6
         BRSN3(I) = DCWSN3(I)/TOTD3
 59   CONTINUE
      ENDIF
C---  OUTPUT TO GMSF AND GMSF3
      GMSF(1)=TOTDN
      GMSF(2)=TOTDL
      GMSF(3)=TOTDR
      GMSF3(1)=TOTD3
      GMSF3(2)=TOTD1
      GMSF3(3)=TOTD2
      do i=1,3
      if(gmsf(i).eq.0D0)  gmsf(i)=5.D-3
      if(gmsf3(i).eq.0D0) gmsf3(i)=5.D-3
      enddo
C--  dummy values for squark WIDTH
      GMSF(4)=GMSF(2)
      GMSF(5)=GMSF(3)
      GMSF(6)=GMSF(2)
      GMSF(7)=GMSF(3)
      GMSF3(4)=GMSF3(2)
      GMSF3(5)=GMSF3(3)
      GMSF3(6)=GMSF3(2)
      GMSF3(7)=GMSF3(3)

 
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
      WRITE(*,*) 'FROM CHARGINO DECAY INOT CHARGED HIGGS'
C  chargino_1 --> neutralino_j + H(-)
      F1MAS = SWM(1)
      DO 72 J = 1,4
        F2MAS = SZM(J)
        SCMAS = SHM(4)
        GL = CONJG(GCNH(2,1,J))
        GR = CONJG(GCNH(1,1,J))
        WRITE(*,*) J,GL,GR
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
      GMSW(1)=TOTD
      IF ( GMSW(1).LE.1.D-10 ) GMSW(1) = 5.D-3
C--     SET DAMMY PARAMETER
      GMSW(2)=GMSW(1)      
 
C--- calculate neutralino_2 2-body decay widthes ( to neutralino )
C  neutralino_2 --> neutralino_1 + Z(-)
        F1MAS = SZM(2)
        F2MAS = SZM(1)
        VMAS = ZM
        GL = GNNZ(1,1,2)
        GR = GNNZ(2,1,2)
        CALL F1F2VD( F1MAS, F2MAS, VMAS, GL, GR, DCYWDT )
        DCWIN2(1) = DCYWDT
 
C  neutralino_2 --> neutralino_1 + h, or H, OR A
        DO 90 J=1,3
        F1MAS = SZM(2)
        F2MAS = SZM(1)
        SCMAS = SHM(J)
        GL = CONJG(GNNH(2,1,2,J))
        GR = CONJG(GNNH(1,1,2,J))
        CALL F1F2SC( F1MAS, F2MAS, SCMAS, GL, GR, DCYWDT )
        DCWIN2(1+J) = DCYWDT
 90     CONTINUE
      TOTD = 0
      DO 91 I=1,4
      IF( DCWIN2(I).LT.0 ) THEN
       WRITE(6,*) ' I=',I,' G3M=',G3MAS
      ENDIF
 91   TOTD = TOTD+DCWIN2(I)
      IF( TOTD.GT.0 ) THEN
        DO 92 I=1,4
 92        BRINO1(I) = DCWIN2(I)/TOTD
      ENDIF
      GMSZ(1)=0
      GMSZ(2)=TOTD
      IF ( GMSZ(2).LE.1.D-10 ) GMSZ(2) = 5.D-3
C   SET DAMMY PARAMETER
      GMSZ(3)=GMSZ(2)
      GMSZ(4)=GMSZ(2)
 
 
      RETURN
      END
