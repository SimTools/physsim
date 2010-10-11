      SUBROUTINE INIPRM(MODE,IFLG,AM0,AMU,AM2,TNB,
     .                  AMZ,AM,AMX,CX,CGE,CZE,CGSE,CZSE)
 
      IMPLICIT   REAL*4  ( A-H, O-Z )
      REAL   *8  AM0, AMU, AM2, TNB, AMZ, AM(2), AMX(4),
     .           CGE(2), CZE(2), CGSE(2), CZSE(2)
      COMPLEX*16 CX(2,2,4)
C--
      REAL*8     BETA, ALF, S2W, AMSF2(7)
C--
      REAL*4     OL(2,2), OR(2,2), ER
      REAL*4     ON(4,4)
      COMPLEX*8  ETA(4), UN(4,4)
      REAL*4     SWM(2), SZM(4)
C--
      COMPLEX*8  GCESNL(2,2), GCNSEL(2,2), GCNSER(2,2),
     .           GCUSDL(2,2), GCUSDR(2,2), GCDSUL(2,2), GCDSUR(2,2)
      COMPLEX*8  GNNSNL(2,4), GNESEL(2,4), GNESER(2,4),
     .           GNUSUL(2,4), GNUSUR(2,4), GNDSDL(2,4), GNDSDR(2,4)
C
C==================< Entry Point >======================================
C
C--
C  Set MSSM parameters.
C--
      ZM    = 91.17
      WM    = 80.0
      SIN2W = 0.23
      ALPHE = 1/128.
      ALPHS = 0.12
      ALF   = ALPHE
      S2W   = SIN2W
      AMZ   = ZM
C--
      SM    = AM0
      HM1   = AMU
      G2MAS = AM2
      G3MAS = G2MAS*(ALPHS/ALPHE)*SIN2W
      TANB  = TNB
      BETA  = ATAN(TNB)
C--
C  Calculate sfermion masses.
C     MODE = 1 : (LL) = (se+,se-)
C          = 2 : (RR)
C          = 3 : (LR)
C          = 4 : (RL)
C--
      CALL SFMASS(AM0,AM2,AMU,BETA,ALF,S2W,AMZ,AMSF2)
C--
      IF ( MODE.EQ.1 ) THEN
         AM(1) = SQRT(AMSF2(5))
         AM(2) = AM(1)
      ELSE IF ( MODE.EQ.2 ) THEN
         AM(1) = SQRT(AMSF2(7))
         AM(2) = AM(1)
      ELSE IF ( MODE.EQ.3 ) THEN
         AM(1) = SQRT(AMSF2(7))
         AM(2) = SQRT(AMSF2(5))
      ELSE IF ( MODE.EQ.4 ) THEN
         AM(1) = SQRT(AMSF2(5))
         AM(2) = SQRT(AMSF2(7))
      ELSE
         PRINT *, ' Invalid MODE = ', MODE
         STOP
      ENDIF
C--
C  Diagonalize mass matrices for inos.
C--
      CALL  INOMIX( ZM, WM, SIN2W, ALPHE, ALPHS,
     .              HM1, G3MAS, TANB,
     .              SWM, SZM,
     .              OL, OR, ER, PHIL, PHIR,
     .              ON, UN, ETA,
     .              IFLG  )
C--
      AMX(1) = SZM(1)
      AMX(2) = SZM(2)
      AMX(3) = SZM(3)
      AMX(4) = SZM(4)
C--
C  Calculate ino-sfermion-fermion couplings.
C--
      CALL  INOFSF( ALPHE, SIN2W, WM,
     .              EMAS, UMAS, DMAS,
     .              TANB, PHIL, PHIR, ER, ON, ETA,
     .              GCESNL, GCNSEL, GCNSER,
     .              GCUSDL, GCUSDR, GCDSUL, GCDSUR,
     .              GNNSNL, GNESEL, GNESER,
     .              GNUSUL, GNUSUR, GNDSDL, GNDSDR )
C--
      DO 10 IX = 1, 4
         CX(1,1,IX) = GNESEL(1,IX)
         CX(2,1,IX) = GNESEL(2,IX)
         CX(1,2,IX) = GNESER(1,IX)
         CX(2,2,IX) = GNESER(2,IX)
10    CONTINUE
C--
C  Calculate sfermion-sfermion-gauge boson couplings.
C--
      PI    = ACOS(-1.D0)
      C2W   = 1 - S2W
      SW    = SQRT(S2W)
      CW    = SQRT(C2W)
C--
      C0      = SQRT(4*PI*ALF)
      CGE(1)  = C0
      CGE(2)  = C0
      CZE(1)  = -C0/(SW*CW)*(-0.5D0+S2W)
      CZE(2)  = -C0/(SW*CW)*(       S2W)
C--
      CGSE(1) = CGE(1)
      CGSE(2) = CGE(2)
      CZSE(1) = CZE(1)
      CZSE(2) = CZE(2)
C--
C  That's it.
C--
      RETURN
      END
