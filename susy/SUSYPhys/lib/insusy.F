C*     FM (1)     :(R*4): e mass.
C*        (2)     :(R*4): u mass.
C*        (3)     :(R*4): d mass.
C*     SFM(1)     :(R*4): snu_L mass.
C*        (2)     :(R*4): se_L mass.
C*        (3)     :(R*4): se_R mass.
C*        (4)     :(R*4): su_L mass.
C*        (5)     :(R*4): su_R mass.
C*        (6)     :(R*4): sd_L mass.
C*        (7)     :(R*4): sd_R mass.
C*     SWM(1-2)   :(R*4): chargino masses.
C*     SZM(1-4)   :(R*4): neutralino masses.
 
      SUBROUTINE INSUSY(AM0,AMU,AM2,TANB,AMA,ALF,ALFS,S2W,ZM,FM,
     .                  SFM,SWM,SZM,GMSF,GMSW,GMSZ,
     .                  GNCW,GCCZ,GNNZ,
     .                  GCESNL,GCNSEL,GCNSER,
     .                  GCUSDL,GCUSDR,GCDSUL,GCDSUR,
     .                  GNNSNL,GNESEL,GNESER,
     .                  GNUSUL,GNUSUR,GNDSDL,GNDSDR )
 
      IMPLICIT   REAL*4  ( A-H, O-Z )
      REAL   *4  AM0, AMU, AM2, TANB, AMA, ALF, ALFS, S2W, ZM, FM(3),
     .           SFM(7), SWM(2), SZM(4), GMSF(7), GMSW(2), GMSZ(4)
      COMPLEX*8  GNCW(2,4,2), GCCZ(2,2,2), GNNZ(2,4,4)
      COMPLEX*8  GCESNL(2,2), GCNSEL(2,2), GCNSER(2,2),
     .           GCUSDL(2,2), GCUSDR(2,2), GCDSUL(2,2), GCDSUR(2,2)
      COMPLEX*8  GNNSNL(2,4), GNESEL(2,4), GNESER(2,4),
     .           GNUSUL(2,4), GNUSUR(2,4), GNDSDL(2,4), GNDSDR(2,4)
C--
      REAL   *4  OL(2,2), OR(2,2), ER
      REAL   *4  ON(4,4)
      COMPLEX*8  ETA(4), UN(4,4)
C--
      REAL   *4  DCWSEL(6), DCWSER(6), DCWSNL(6)
      REAL   *4  BRSEL(6), BRSER(6), BRSNL(6)
      REAL   *4  DCWIN1(8), BRINO1(8)
      REAL   *4  HG0MAS(2)
C
C==================< Entry Point >======================================
C
C--
C  Set MSSM parameters.
C--
      WM    = ZM*SQRT(1-S2W)
C--
      SM    = AM0
      HM1   = AMU
      G2MAS = AM2
      G3MAS = G2MAS*(ALFS/ALF)*S2W
C--
C  Calculate sfermion masses.
C     SFM(I) = (SNL,SEL,SER,SUL,SUR,SDL,SDR)
C--
      CALL GTSFMS(AM0,AM2,AMU,TANB,ALF,S2W,ZM,SFM)
C--
C  Diagonalize mass matrices for inos.
C--
      CALL  INOMIX( ZM, WM, S2W, ALF, ALFS,
     .              HM1, G3MAS, TANB,
     .              SWM, SZM,
     .              OL, OR, ER, PHIL, PHIR,
     .              ON, UN, ETA,
     .              IFLG  )
C--
C  Calculate ino-ino-gauge boson couplings.
C--
      CALL  INOGBS( ALF, S2W,
     .              PHIL, PHIR, ER, ON, ETA,
     .              GNCW, GCCZ, GNNZ )
C--
C  Calculate ino-sfermion-fermion couplings.
C--
      CALL  INOFSF( ALF, S2W, WM,
     .              FM(1), FM(2), FM(3),
     .              TANB, PHIL, PHIR, ER, ON, ETA,
     .              GCESNL, GCNSEL, GCNSER,
     .              GCUSDL, GCUSDR, GCDSUL, GCDSUR,
     .              GNNSNL, GNESEL, GNESER,
     .              GNUSUL, GNUSUR, GNDSDL, GNDSDR )
C--
C  Calculate widths of SUSY particles.
C  Very temporary version which is by no means correct.
C--
      CALL SF2BR(  WM, S2W, ALF, ALFS,
     .             HM1, G3MAS, TANB,
     .             FM(1), FM(2), FM(3),
     .             SFM(2), SFM(3), SFM(1), AMA,
     .             SWM, SZM, HGCMAS, HG0MAS,
     .             DCWSEL, DCWSER, DCWSNL,
     .             BRSEL, BRSER, BRSNL,
     .             DCWIN1, BRINO1 )
C--
      GMSW(1) = DCWIN1(1) + DCWIN1(2) + DCWIN1(3) + DCWIN1(4)
     .        + DCWIN1(5) + DCWIN1(6) + DCWIN1(7) + DCWIN1(8)
      IF ( GMSW(1).LE.1.E-10 ) GMSW(1) = 5.E-3
      GMSW(1) = GMSW(1)
      GMSW(2) = GMSW(1)
      GMSZ(1) = 0
      GMSZ(2) = GMSW(1)
      GMSZ(3) = GMSW(1)
      GMSZ(4) = GMSW(1)
C--
      GMSF(1) = DCWSNL(1) + DCWSNL(2) + DCWSNL(3)
     .        + DCWSNL(4) + DCWSNL(5) + DCWSNL(6)
      GMSF(2) = DCWSEL(1) + DCWSEL(2) + DCWSEL(3)
     .        + DCWSEL(4) + DCWSEL(5) + DCWSEL(6)
      GMSF(3) = DCWSER(1) + DCWSER(2) + DCWSER(3)
     .        + DCWSER(4) + DCWSER(5) + DCWSER(6)
      GMSF(4) = GMSF(2)
      GMSF(5) = GMSF(3)
      GMSF(6) = GMSF(2)
      GMSF(7) = GMSF(3)
C--
C  That's it.
C--
      RETURN
      END
