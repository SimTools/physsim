CC**********************************************************************
C*
C*====================================================
C* Subroutine AMPNNZZ(GZL,GZN, GWF, GWWZ, GWWH,GZZH,
C*                    AMW,GMW, AMZ,GMZ, AME, AMH,GMH,
C*                       SPEM,SPEP, SP1,SP2,VCT, AMP)
C*===============================================---==
C*
C* (Purpose)
C*    Calculate t-channel amplitude for e- e+ ---> nu nubar Z Z.
C*    Inputs are spinors and wave functions thus replaceable with
C*    spinors or wave functions made from final stable particles.
C* (Inputs)
C*       GZL(*)   : (R*8) : lepton Z coupling.
C*       GZN(*)   : (R*8) : neutrino Z coupling.
C*       GWF      : (R*8) : weak coupling.
C*       GWWZ     : (R*8) : WWZ coupling.
C*       GWWH     : (R*8) : WWH coupling.
C*       GZZH     : (R*8) : ZZH coupling.
C*       AMW      : (R*8) : W mass.
C*       GMW      : (R*8) : W width.
C*       AMZ      : (R*8) : Z mass.
C*       GMZ      : (R*8) : Z width.
C*       AME      : (R*8) : electron mass.
C*       AMH      : (R*8) : H mass.
C*       GMH      : (R*8) : H width.
C*       SPEM(*)  : (C*16) : beam e- spinor.
C*       SPEP(*)  : (C*16) : beam e+ spinor.
C*       SP1 (*)  : (C*16) : final-state nu spinor.
C*       SP2 (*)  : (C*16) : final-state nub spinor.
C*       VCT(*,1) : (C*16) : final-state Z1 wave function.
C*       VCT(*,2) : (C*16) : final-state Z2 wave function.
C* (Output)
C*       AMP(0)   : (C*16) : total t-channel amplitudes.
C*          (1)   : (C*16) : non-fusion diagrams.
C*          (2)   : (C*16) : fusion diagrams.
C*          (3)   : (C*16) : Higgs diagram.
C* (Relation)
C*    Invokes subroutines in HELAS.LOAD.
C* (Update Record)
C*    2005/03/04  K.Fujii	Original version.
C*
CC**********************************************************************
 

      SUBROUTINE AMPNNZZ(GZL,GZN,GWF,GWWZ,  GWWH, GZZH,
     .                     AMW,GMW, AMZ,GMZ, AME, AMH,GMH,   
     .                        SPEM,SPEP, SP1,SP2,VCT, AMP)
 
      IMPLICIT   REAL*8  ( A-H, O-Z )
      REAL   *8  GZL(2), GZN(2), GWF(2), GWWZ, GWWH, GZZH,
     .           AMW, GMW, AMZ, GMZ, AME, AMH, GMH
      COMPLEX*16 SPEM(6), SPEP(6), SP1(6), SP2(6), VCT(6,2), 
     .           AMP(0:3)
C--
      COMPLEX*16 WRK(6,18), TMP(6), ZMP(30)
C--
      DATA NCALL /0/
C
C========< Entry Point >================================================
C
C--
C  Initialization.
C--
      IF ( NCALL.EQ.0 ) THEN
         NCALL = 1
         AMN   = 0
         GMN   = 0
         GME   = 0
      ENDIF
C--
C  Calculate internal lines.
C--
C-- virtual electron with 1 Z
      CALL FVIXXX(SPEM     ,VCT(1,1),GZL,AME,GME, WRK(1, 1))
      CALL FVIXXX(SPEM     ,VCT(1,2),GZL,AME,GME, WRK(1, 2))
C-- virtual electron with 2 Z's
      CALL FVIXXX(WRK(1, 1),VCT(1,2),GZL,AME,GME, WRK(1, 3))
      CALL FVIXXX(WRK(1, 2),VCT(1,1),GZL,AME,GME, WRK(1, 4))
C-- virtual positron with 1 Z
      CALL FVOXXX(SPEP     ,VCT(1,1),GZL,AME,GME, WRK(1, 5))
      CALL FVOXXX(SPEP     ,VCT(1,2),GZL,AME,GME, WRK(1, 6))
C-- virtual positron with 2 Z's
      CALL FVOXXX(WRK(1, 5),VCT(1,2),GZL,AME,GME, WRK(1, 7))
      CALL FVOXXX(WRK(1, 6),VCT(1,1),GZL,AME,GME, WRK(1, 8))
C-- virtual nu with 1 Z
      CALL FVOXXX(SP1      ,VCT(1,1),GZN,AMN,GMN, WRK(1, 9))
      CALL FVOXXX(SP1      ,VCT(1,2),GZN,AMN,GMN, WRK(1,10))
C-- virtual nu with 2 Z's
      CALL FVOXXX(WRK(1, 9),VCT(1,2),GZN,AMN,GMN, WRK(1,11))
      CALL FVOXXX(WRK(1,10),VCT(1,1),GZN,AMN,GMN, WRK(1,12))
C-- virtual nubar with 1 Z
      CALL FVIXXX(SP2      ,VCT(1,1),GZN,AMN,GMN, WRK(1,13))
      CALL FVIXXX(SP2      ,VCT(1,2),GZN,AMN,GMN, WRK(1,14))
C-- virtual nubar with 2 Z's
      CALL FVIXXX(WRK(1,13),VCT(1,2),GZN,AMN,GMN, WRK(1,15))
      CALL FVIXXX(WRK(1,14),VCT(1,1),GZN,AMN,GMN, WRK(1,16))
C-- virtual W
      CALL JIOXXX(SPEM, SP1, GWF,AMW,GMW, WRK(1,17))
      CALL JIOXXX(SP2, SPEP, GWF,AMW,GMW, WRK(1,18))

C--
C  Now calculate amplitudes.
C--
C  Non-fusion (e+e-->nnbar + 2 Z's to fermion legs).
C--
      CALL AMPNN(GWF,AMW,GMW,WRK(1, 3),SPEP,SP1,SP2,ZMP( 1))
      CALL AMPNN(GWF,AMW,GMW,WRK(1, 4),SPEP,SP1,SP2,ZMP( 2))

      CALL AMPNN(GWF,AMW,GMW,WRK(1, 1),WRK(1, 6), SP1,SP2,ZMP( 3))
      CALL AMPNN(GWF,AMW,GMW,WRK(1, 2),WRK(1, 5), SP1,SP2,ZMP( 4))

      CALL AMPNN(GWF,AMW,GMW,WRK(1, 1),SPEP,WRK(1,10),SP2,ZMP( 5))
      CALL AMPNN(GWF,AMW,GMW,WRK(1, 2),SPEP,WRK(1, 9),SP2,ZMP( 6))

      CALL AMPNN(GWF,AMW,GMW,WRK(1, 1),SPEP,SP1,WRK(1,14),ZMP( 7))
      CALL AMPNN(GWF,AMW,GMW,WRK(1, 2),SPEP,SP1,WRK(1,13),ZMP( 8))

      CALL AMPNN(GWF,AMW,GMW,SPEM,WRK(1, 7),SP1,SP2,ZMP( 9))
      CALL AMPNN(GWF,AMW,GMW,SPEM,WRK(1, 8),SP1,SP2,ZMP(10))

      CALL AMPNN(GWF,AMW,GMW,SPEM,WRK(1, 5),WRK(1,10),SP2,ZMP(11))
      CALL AMPNN(GWF,AMW,GMW,SPEM,WRK(1, 6),WRK(1, 9),SP2,ZMP(12))

      CALL AMPNN(GWF,AMW,GMW,SPEM,WRK(1, 5),SP1,WRK(1,14),ZMP(13))
      CALL AMPNN(GWF,AMW,GMW,SPEM,WRK(1, 6),SP1,WRK(1,13),ZMP(14))

      CALL AMPNN(GWF,AMW,GMW,SPEM,SPEP,WRK(1,11),SP2,ZMP(15))
      CALL AMPNN(GWF,AMW,GMW,SPEM,SPEP,WRK(1,12),SP2,ZMP(16))

      CALL AMPNN(GWF,AMW,GMW,SPEM,SPEP,WRK(1, 9),WRK(1,14),ZMP(17))
      CALL AMPNN(GWF,AMW,GMW,SPEM,SPEP,WRK(1,10),WRK(1,13),ZMP(18))

      CALL AMPNN(GWF,AMW,GMW,SPEM,SPEP,SP1,WRK(1,15),ZMP(19))
      CALL AMPNN(GWF,AMW,GMW,SPEM,SPEP,SP1,WRK(1,16),ZMP(20))
C--
C  Non-fusion (1 Z from fermion legs and another from internal W)
C--
      CALL JIOXXX(WRK(1, 1),  SP1, GWF,AMW,GMW, TMP)
      CALL VVVXXX(WRK(1,18),TMP,VCT(1, 2),GWWZ, ZMP(21))
      CALL JIOXXX(WRK(1, 2),  SP1, GWF,AMW,GMW, TMP)
      CALL VVVXXX(WRK(1,18),TMP,VCT(1, 1),GWWZ, ZMP(22))

      CALL JIOXXX(SP2,  WRK(1, 5), GWF,AMW,GMW, TMP)
      CALL VVVXXX(TMP,WRK(1,17),VCT(1, 2),GWWZ, ZMP(23))
      CALL JIOXXX(SP2,  WRK(1, 6), GWF,AMW,GMW, TMP)
      CALL VVVXXX(TMP,WRK(1,17),VCT(1, 1),GWWZ, ZMP(24))

      CALL JIOXXX(SPEM, WRK(1, 9), GWF,AMW,GMW, TMP)
      CALL VVVXXX(WRK(1,18),TMP,VCT(1, 2),GWWZ, ZMP(25))
      CALL JIOXXX(SPEM, WRK(1,10), GWF,AMW,GMW, TMP)
      CALL VVVXXX(WRK(1,18),TMP,VCT(1, 1),GWWZ, ZMP(26))

      CALL JIOXXX(WRK(1,13), SPEP, GWF,AMW,GMW, TMP)
      CALL VVVXXX(TMP,WRK(1,17),VCT(1, 2),GWWZ, ZMP(27))
      CALL JIOXXX(WRK(1,14), SPEP, GWF,AMW,GMW, TMP)
      CALL VVVXXX(TMP,WRK(1,17),VCT(1, 1),GWWZ, ZMP(28))
C--
C  Fusion.
C--
      CALL W3W3XX(WRK(1,18),VCT(1, 1),WRK(1,17),VCT(1, 2),
     .                         GWWZ,GWWZ,AMW,GMW, ZMP(29))
C--
C  Higgs.
C--
      CALL HVVXXX(VCT(1, 1),VCT(1,  2),  GZZH,AMH,GMH, TMP)
      CALL VVSXXX(WRK(1,18),WRK(1, 17),TMP    ,GWWH,ZMP(30))
C--
C  Sum up amplitudes.
C--
      AMP(1) = 0
      DO 10 I = 1, 28 
         AMP(1) = AMP(1) + ZMP(I)
10    CONTINUE
      AMP(2) = ZMP(29)
      AMP(3) = ZMP(30)
C--
      AMP(0) = AMP(1) + AMP(2) + AMP(3)
C--
C  That's it.
C--
      RETURN
      END
