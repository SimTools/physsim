CC**********************************************************************
C*
C*============================================================
C* Subroutine AMPNNWW(GAL,GZL,GZN,GWF,GWWA,GWWZ,GW,GWWH,GZZH,
C*                            AMW,GMW, AMZ,GMZ, AME, AMH,GMH,
C*                               SPEM,SPEP, SP1,SP2,VCT, AMP)
C*=======================================================---==
C*
C* (Purpose)
C*    Calculate t-channel amplitude for e- e+ ---> nu nubar W- W+.
C*    Inputs are spinors and wave functions thus replaceable with
C*    spinors or wave functions made from final stable particles.
C* (Inputs)
C*       GAL(*)   : (R*8) : lepton photon coupling.
C*       GZL(*)   : (R*8) : lepton Z coupling.
C*       GZN(*)   : (R*8) : neutrino Z coupling.
C*       GWF      : (R*8) : weak coupling.
C*       GWWA     : (R*8) : WWA coupling.
C*       GWWZ     : (R*8) : WWZ coupling.
C*       GW       : (R*8) : gw.
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
C*       VCT(*,1) : (C*16) : final-state W- wave function.
C*       VCT(*,2) : (C*16) : final-state W+ wave function.
C* (Output)
C*       AMP(0)   : (C*16) : total t-channel amplitudes.
C*          (1)   : (C*16) : non-fusion diagrams.
C*          (2)   : (C*16) : fusion diagrams.
C*          (3)   : (C*16) : Higgs diagram.
C* (Relation)
C*    Invokes subroutines in HELAS.LOAD.
C* (Update Record)
C*    99/05/21  K.Fujii		Original version.
C*
CC**********************************************************************
 

      SUBROUTINE AMPNNWW(GAL,GZL,GZN,GWF,GWWA,GWWZ, GW, GWWH, GZZH,
     .                              AMW,GMW, AMZ,GMZ, AME, AMH,GMH,   
     .                                 SPEM,SPEP, SP1,SP2,VCT, AMP)
 
      IMPLICIT   REAL*8  ( A-H, O-Z )
      REAL   *8  GAL(2), GZL(2), GZN(2), GWF(2), GWWA, GWWZ, 
     .           GW, GWWH, GZZH,
     .           AMW, GMW, AMZ, GMZ, AME, AMH, GMH
      COMPLEX*16 SPEM(6), SPEP(6), SP1(6), SP2(6), VCT(6,2), 
     .		 AMP(0:3)
C--
      COMPLEX*16 WRK(6,40), TMP(6), ZMP(20)
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
         AMA   = 0
         GMN   = 0
         GME   = 0
         GMA   = 0
      ENDIF
C--
C  Calculate internal lines.
C--
      CALL JIOXXX(SPEM,SP1     ,GWF,AMW,GMW, WRK(1, 1))
      CALL FVOXXX(SPEP,WRK(1,1),GWF,AMN,GMN, WRK(1, 7))
      CALL FVIXXX(SP2 ,WRK(1,1),GWF,AME,GME, WRK(1, 9))
C--
      CALL JIOXXX(SP2 ,SPEP    ,GWF,AMW,GMW, WRK(1, 2))
      CALL FVIXXX(SPEM,WRK(1,2),GWF,AMN,GMN, WRK(1, 8))
      CALL FVOXXX(SP1 ,WRK(1,2),GWF,AME,GME, WRK(1,10))
C--
      CALL FVIXXX(SPEM,VCT(1,1),GWF,AMN,GMN, WRK(1, 3))
      CALL JIOXXX(WRK(1,3),SP1 ,GZN,AMZ,GMZ, WRK(1,11))
C--
      CALL FVOXXX(SPEP,VCT(1,2),GWF,AMN,GMN, WRK(1, 4))
      CALL JIOXXX(SP2 ,WRK(1,4),GZN,AMZ,GMZ, WRK(1,12))
C--
      CALL FVOXXX(SP1 ,VCT(1,1),GWF,AME,GME, WRK(1, 5))
C--
      CALL FVIXXX(SP2 ,VCT(1,2),GWF,AME,GME, WRK(1, 6))
C--
C  Now calculate amplitudes.
C    The amplitude numbers given here correspond to those on P86 of
C    H.Murayama's Ph.D thesis.
C--
C  Non-fusion.
C--
C( 1)
      CALL JIOXXX(SPEM     ,WRK(1, 5),    GAL,AMA,GMA, TMP)
      CALL IOVXXX(WRK(1, 6),SPEP     ,TMP      ,GAL, ZMP( 1))
      CALL JIOXXX(SPEM     ,WRK(1, 5),    GZL,AMZ,GMZ, TMP)
      CALL IOVXXX(WRK(1, 6),SPEP     ,TMP      ,GZL, ZMP( 2))
      ZMP( 1) = ZMP( 1) + ZMP( 2)
C( 2)
      CALL IOVXXX(WRK(1, 6),SPEP     ,WRK(1,11),GZL, ZMP( 2))
C( 3)
      CALL IOVXXX(SPEM     ,WRK(1, 5),WRK(1,12),GZL, ZMP( 3))
C( 4)
      CALL IOVXXX(SP2      ,WRK(1, 4),WRK(1,11),GZN, ZMP( 4))
C( 5)
      CALL JIOXXX(WRK(1, 8),SP1      ,    GZN,AMZ,GMZ, TMP)
      CALL VVVXXX(VCT(1, 1),VCT(1, 2),TMP      ,GWWZ,ZMP( 5))
C( 6)
      CALL J3XXXX(SPEM     ,WRK(1,10),GAL,GZL,AMZ,GMZ, TMP)
      CALL VVVXXX(VCT(1, 1),VCT(1, 2),TMP      ,GW,  ZMP( 6))
C( 7)
      CALL J3XXXX(WRK(1, 9),SPEP     ,GAL,GZL,AMZ,GMZ, TMP)
      CALL VVVXXX(VCT(1, 1),VCT(1, 2),TMP      ,GW,  ZMP( 7))
C( 8)
      CALL JIOXXX(SP2      ,WRK(1, 7),    GZN,AMZ,GMZ, TMP)
      CALL VVVXXX(VCT(1, 1),VCT(1, 2),TMP      ,GWWZ,ZMP( 8))
C( 9)
      CALL IOVXXX(WRK(1, 8),WRK(1, 5),VCT(1, 2),GWF, ZMP( 9))
C(10)
      CALL IOVXXX(WRK(1, 3),WRK(1,10),VCT(1, 2),GWF, ZMP(10))
C(11)
      CALL IOVXXX(WRK(1, 9),WRK(1, 4),VCT(1, 1),GWF, ZMP(11))
C(12)
      CALL IOVXXX(WRK(1, 6),WRK(1, 7),VCT(1, 1),GWF, ZMP(12))
C(13)
      CALL J3XXXX(SPEM     ,WRK(1, 5),GAL,GZL,AMZ,GMZ, TMP)
      CALL VVVXXX(WRK(1, 2),VCT(1, 2),TMP      ,GW,  ZMP(13))
C(14)
      CALL VVVXXX(WRK(1, 2),VCT(1, 2),WRK(1,11),GWWZ,ZMP(14))
C(15)
      CALL VVVXXX(VCT(1, 1),WRK(1, 1),WRK(1,12),GWWZ,ZMP(15))
C(16)
      CALL J3XXXX(WRK(1, 6),SPEP     ,GAL,GZL,AMZ,GMZ, TMP)
      CALL VVVXXX(VCT(1, 1),WRK(1, 1),TMP      ,GW,  ZMP(16))
C--
C  Fusion.
C--
C(17-19)
      CALL WWWWXX(VCT(1, 1),VCT(1, 2),WRK(1, 2),WRK(1, 1),
     .            		  GWWA,GWWZ,AMZ,GMZ, ZMP(17))
C--
C  Higgs.
C--
C(20)
      CALL HVVXXX(VCT(1, 1),VCT(1, 2),   GWWH,AMH,GMH, TMP)
      CALL VVSXXX(WRK(1, 1),WRK(1, 2),TMP      ,GWWH,ZMP(18))
C(21)
      CALL HVVXXX(VCT(1, 1),WRK(1, 1),   GWWH,AMH,GMH, TMP)
      CALL VVSXXX(WRK(1, 2),VCT(1, 2),TMP      ,GWWH,ZMP(19))
C--
C  Sum up amplitudes.
C--
      AMP(1) = 0
      DO 10 I = 1, 16
         AMP(1) = AMP(1) + ZMP(I)
10    CONTINUE
      AMP(2) = ZMP(17)
      AMP(3) = ZMP(18) + ZMP(19)
C--
      AMP(0) = AMP(1) + AMP(2) + AMP(3)
C--
C  That's it.
C--
      RETURN
      END
