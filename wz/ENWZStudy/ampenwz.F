CC**********************************************************************
C*
C*============================================================
C* Subroutine AMPENWZ(GAL,GZL,GZN,GWF,GWWA,GWWZ,GW,GWWH,GZZH,
C*                     AMW,GMW, AMZ,GMZ, AME, AMH,GMH,   EBM,
C*  		             E1, SHLF1,CHLF1,FI1, NHEB1,NHE1,
C*                               SPEM,SPEP, SP1,SP2,VCT, AMP)
C*=======================================================---==
C*
C* (Purpose)                                          __
C*    Calculate t-channel amplitude for e+ e- ---> e- nu_e W+ Z.
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
C*       EBM      : (R*8) : beam energy.
C*       E1       : (R*8) : final-state e- energy.
C*       SHLF1    : (R*8) : sin(theta_e-/2).
C*       CHLF1    : (R*8) : cos(theta_e-/2).
C*       FI1      : (R*8) : phi_e-
C*       NHEB1    : (I*4) : beam e- helicity.
C*       NHE1     : (I*4) : final-state e- helicity.
C*       SPEM(*)  : (C*16) : beam e- spinor.
C*       SPEP(*)  : (C*16) : beam e+ spinor.
C*       SP1 (*)  : (C*16) : final-state e- spinor.
C*       SP2 (*)  : (C*16) : final-state nu_e bar spinor.
C*       VCT(*,1) : (C*16) : final-state W+ wave function.
C*       VCT(*,2) : (C*16) : final-state Z  wave function.
C* (Output)
C*       AMP(0)   : (C*16) : total t-channel amplitudes.
C*          (1)   : (C*16) : non-fusion diagrams.
C*          (2)   : (C*16) : fusion diagrams.
C*          (3)   : (C*16) : Higgs diagram.
C* (Relation)
C*    Invokes subroutines in HELAS.LOAD.
C* (Update Record)
C*    2008/11/27  K.Fujii		Original version.
C*
CC**********************************************************************
 

      SUBROUTINE AMPENWZ(GAL,GZL,GZN,GWF,GWWA,GWWZ, GW, GWWH, GZZH,
     .                       AMW,GMW, AMZ,GMZ, AME, AMH,GMH,   EBM,
     .                       E1, SHLF1,CHLF1,FI1, NHEB1,NHE1,
     .                                 SPEM,SPEP, SP1,SP2,VCT, AMP)
 
      IMPLICIT   REAL*8  ( A-H, O-Z )
      REAL   *8  GAL(2), GZL(2), GZN(2), GWF(2), GWWA, GWWZ, 
     .           GW, GWWH, GZZH,
     .           AMW, GMW, AMZ, GMZ, AME, AMH, GMH, EBM,
     .            E1, SHLF1, CHLF1, FI1
      INTEGER*4  NHEB1, NHE1
      COMPLEX*16 SPEM(6), SPEP(6), SP1(6), SP2(6), VCT(6,2), 
     .           AMP(0:3)
C--
      COMPLEX*16 WRK(6,40), TMP(6), ZMP(25)
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
         PI    = ACOS(-1.D0)
         AMN   = 0
         GMN   = 0
         GME   = 0
         AMA   = 0
         GMA   = 0
      ENDIF
C--
C  Calculate internal lines.
C--
      CALL JEEXXX(EBM,E1,SHLF1,CHLF1,FI1,NHEB1,NHE1,+1, WRK(1, 1))
      CALL JIOXXX(SPEM,SP1,GZL,AMZ,GMZ,                 WRK(1, 2))
C--
      CALL JIOXXX(SP2,SPEP,GWF,AMW,GMW,                 WRK(1, 3))
C--
      CALL FVOXXX(SP1,VCT(1,2), GZL,AME,GME, WRK(1, 4))
      CALL JIOXXX(SPEM,WRK(1,4),GAL,AMA,GMA, WRK(1, 5))
      CALL JIOXXX(SPEM,WRK(1,4),GZL,AMZ,GMZ, WRK(1, 6))
C--
      CALL FVIXXX(SP2,VCT(1,2), GZN,AMN,GMN, WRK(1, 7))
      CALL FVIXXX(SP2,VCT(1,1), GWF,AME,GME, WRK(1, 8))
C--
      CALL FVOXXX(SP1 ,VCT(1,1),GWF,AMN,GMN, WRK(1, 9))
      CALL JIOXXX(SPEM,WRK(1,9),GWF,AMW,GMW, WRK(1,10))
C--
      CALL FVOXXX(SPEP,VCT(1,2),GZL,AME,GME, WRK(1,11))
      CALL JIOXXX(SP2,WRK(1,11),GWF,AMW,GMW, WRK(1,12))
C--
      CALL FVIXXX(WRK(1,7),VCT(1,1),GWF,AME,GME,WRK(1,13))
      CALL FVIXXX(WRK(1,8),VCT(1,2),GZL,AME,GME,WRK(1,14))
C--
      CALL JVVXXX(VCT(1,1),VCT(1,2),GWWZ,AMW,GMW,WRK(1,15))
      CALL FVIXXX(SP2,    WRK(1,15),GWF ,AME,GME,WRK(1,16))
C--
C  Calculate W3 combination.
C
      WRK(1,17) = WRK(1,1)*(GWWA/GW) + WRK(1,2)*(GWWZ/GW)
      WRK(2,17) = WRK(2,1)*(GWWA/GW) + WRK(2,2)*(GWWZ/GW)
      WRK(3,17) = WRK(3,1)*(GWWA/GW) + WRK(3,2)*(GWWZ/GW)
      WRK(4,17) = WRK(4,1)*(GWWA/GW) + WRK(4,2)*(GWWZ/GW)
      WRK(5,17) = WRK(5,1)
      WRK(6,17) = WRK(6,1)
C--
      CALL FVOXXX(SPEP,VCT(1,1),GWF,AMN,GMN, WRK(1,18))
C--
      CALL FVIXXX(SPEM,VCT(1,2),GZL,AME,GME, WRK(1,19))
C--
      CALL FVIXXX(SPEM,WRK(1,3),GWF,AMN,GMN, WRK(1,20))
C--
      CALL FVIXXX(SP2 ,WRK(1,2),GZN,AMN,GMN, WRK(1,21))
C--
      CALL JIOXXX(WRK(1,19),SP1,GAL,AMA,GMA, WRK(1,22))
      CALL JIOXXX(WRK(1,19),SP1,GZL,AMZ,GMZ, WRK(1,23))
C--
C  Now calculate amplitudes.
C    The amplitude numbers given here correspond to those on P86 of
C    H.Murayama's Ph.D thesis.
C--
C  Non-fusion.
C--
C( 1)
      CALL IOVXXX(WRK(1, 8),WRK(1,11),WRK(1, 1),GAL, TMP( 1))
      CALL IOVXXX(WRK(1, 8),WRK(1,11),WRK(1, 2),GZL, ZMP( 1))
      ZMP(1) = ZMP(1) + TMP(1)
C( 2)
      CALL IOVXXX(WRK(1,13),     SPEP,WRK(1, 1),GAL, TMP( 1))
      CALL IOVXXX(WRK(1,13),     SPEP,WRK(1, 2),GZL, ZMP( 2))
      ZMP(2) = ZMP(2) + TMP(1)
C( 3)
      CALL IOVXXX(WRK(1,14),     SPEP,WRK(1, 1),GAL, TMP( 1))
      CALL IOVXXX(WRK(1,14),     SPEP,WRK(1, 2),GZL, ZMP( 3))
      ZMP(3) = ZMP(3) + TMP(1)
C( 4)
      CALL IOVXXX(WRK(1,16),     SPEP,WRK(1, 1),GAL, TMP( 1))
      CALL IOVXXX(WRK(1,16),     SPEP,WRK(1, 2),GZL, ZMP( 4))
      ZMP(4) = ZMP(4) + TMP(1)
C( 5)
      CALL JIOXXX(WRK(1, 7), SPEP,GWF,AMW,GMW, TMP)
      CALL VVVXXX(TMP      , VCT(1,1),WRK(1,17),GW,  ZMP( 5))
C( 6)
      CALL VVVXXX(WRK(1,12), VCT(1,1),WRK(1,17),GW,  ZMP( 6))
C( 7)
      CALL IOVXXX(WRK(1, 7),WRK(1,18),WRK(1, 2),GZN, ZMP( 7))
C( 8)
      CALL IOVXXX(WRK(1,21),WRK(1,18),VCT(1, 2),GZN, ZMP( 8))
C( 9)
      CALL IOVXXX(WRK(1,21),WRK(1,11),VCT(1, 1),GWF, ZMP( 9))
C(10)
      CALL IOVXXX(WRK(1,21),     SPEP,WRK(1,15),GWF, ZMP(10))
C(11)
      CALL IOVXXX(SP2      ,WRK(1,18),WRK(1, 6),GZN, ZMP(11))
C(12)
      CALL IOVXXX(SP2      ,WRK(1,18),WRK(1,23),GZN, ZMP(12))
C(13)
      CALL IOVXXX(WRK(1,19),WRK(1, 9),WRK(1, 3),GWF, ZMP(13))
C(14)
      CALL IOVXXX(WRK(1,20),WRK(1, 4),VCT(1, 1),GWF, ZMP(14))
C(15)
      CALL IOVXXX(WRK(1,20),WRK(1, 9),VCT(1, 2),GZN, ZMP(15))
C(16)
      CALL IOVXXX(WRK(1,20),SP1      ,WRK(1,15),GWF, ZMP(16))
C(17)
      CALL J3XXXX(SPEM     ,WRK(1, 4),GAL,GZL,AMZ,GMZ, TMP)
      CALL VVVXXX(WRK(1, 3),VCT(1, 1),TMP      ,GW,  ZMP(17))
C(18)
      CALL J3XXXX(WRK(1,19),SP1      ,GAL,GZL,AMZ,GMZ, TMP)
      CALL VVVXXX(WRK(1, 3),VCT(1, 1),TMP      ,GW,  ZMP(18))
C(19)
      CALL IOVXXX(WRK(1, 8),SPEP     ,WRK(1, 5),GAL, TMP( 1))
      CALL IOVXXX(WRK(1, 8),SPEP     ,WRK(1, 6),GZL, ZMP(19))
      ZMP(19) = ZMP(19) + TMP(1)
C(20)
      CALL IOVXXX(WRK(1, 8),SPEP     ,WRK(1,22),GAL, TMP( 1))
      CALL IOVXXX(WRK(1, 8),SPEP     ,WRK(1,23),GZL, ZMP(20))
      ZMP(20) = ZMP(20) + TMP(1)
C(21)
      CALL IOVXXX(WRK(1, 7),SPEP     ,WRK(1,10),GWF, ZMP(21))
C(22)
      CALL IOVXXX(SP2      ,WRK(1,11),WRK(1,10),GWF, ZMP(22))
C(23)
      CALL VVVXXX(WRK(1, 3),WRK(1,10),VCT(1, 2),GWWZ,ZMP(23))
C--
C  Fusion.
C--
C(24-26)
      CALL W3W3XX(WRK(1, 3),WRK(1,17),VCT(1, 1),VCT(1, 2),
     .                              GW,GWWZ,AMW,GMW, ZMP(24))
C--
C  Higgs.
C--
C(27)
      CALL HVVXXX(WRK(1, 2),VCT(1, 2),GZZH,AMH,GMH, TMP)
      CALL VVSXXX(WRK(1, 3),VCT(1, 1),TMP,GWWH,      ZMP(25))
C--
C  Sum up amplitudes.
C--
      AMP(1) = 0
      DO 10 I = 1, 23
         AMP(1) = AMP(1) + ZMP(I)
10    CONTINUE
      AMP(2) = ZMP(24)
      AMP(3) = ZMP(25) 
C--
      AMP(0) = AMP(1) + AMP(2) + AMP(3)
C--
C  That's it.
C--
      RETURN
      END
