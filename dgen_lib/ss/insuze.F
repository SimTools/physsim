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
 
      SUBROUTINE INSUZE(AMU,AM1,AM2,TANB,ALF,ALFS,S2W,ZM,FM,
     .                  SWM,SZM,
     .                  GNCW,GCCZ,GNNZ,
     .                  GCESNL,GCNSEL,GCNSER,
     .                  GCUSDL,GCUSDR,GCDSUL,GCDSUR,
     .                  GNNSNL,GNESEL,GNESER,
     .                  GNUSUL,GNUSUR,GNDSDL,GNDSDR )
 
      IMPLICIT   REAL*8  ( A-H, O-Z )
      REAL   *8  AMU, AM1, AM2, TANB, ALF, ALFS, S2W, ZM, FM(3),
     .           SWM(2), SZM(4)
      COMPLEX*16 GNCW(2,4,2), GCCZ(2,2,2), GNNZ(2,4,4)
      COMPLEX*16 GCESNL(2,2), GCNSEL(2,2), GCNSER(2,2),
     .           GCUSDL(2,2), GCUSDR(2,2), GCDSUL(2,2), GCDSUR(2,2)
      COMPLEX*16 GNNSNL(2,4), GNESEL(2,4), GNESER(2,4),
     .           GNUSUL(2,4), GNUSUR(2,4), GNDSDL(2,4), GNDSDR(2,4)
C--
      REAL   *8  OL(2,2), OR(2,2), ER
      REAL   *8  ON(4,4)
      COMPLEX*16 ETA(4), UN(4,4)
C--
      REAL   *8  DCWSEL(6), DCWSER(6), DCWSNL(6)
      REAL   *8  BRSEL(6), BRSER(6), BRSNL(6)
      REAL   *8  DCWIN1(8), BRINO1(8)
      REAL   *8  HG0MAS(2)
C
C==================< Entry Point >======================================
C
C--
C  Set MSSM parameters.
C--
      WM    = ZM*SQRT(1-S2W)
C--
      HM1   = AMU
      G1MAS = AM1
      G2MAS = AM2
C--
C  Diagonalize mass matrices for inos.
C--
      IFLG = 1
      CALL  INOMXX( ZM, WM, S2W, ALF, ALFS,
     .              HM1, G1MAS, G2MAS, TANB,
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
C  That's it.
C--
      RETURN
      END
