CC**********************************************************************
C*
C*============================================================
C* Subroutine AMPEEWW(GAL,GZL,GZN,GWF,GWWA,GWWZ,GW,GWWH,GZZH,
C*                     AMW,GMW, AMZ,GMZ, AME, AMH,GMH,   EBM,
C*  		             E1, SHLF1,CHLF1,FI1, NHEB1,NHE1,
C* 		             E2, SHLF2,CHLF2,FI2, NHEB2,NHE2,
C*                               SPEM,SPEP, SP1,SP2,VCT, AMP)
C*=======================================================---==
C*
C* (Purpose)
C*    Calculate t-channel amplitude for e+ e- ---> e+ e- W+ W-.
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
C*       E2       : (R*8) : final-state e+ energy.
C*       SHLF2    : (R*8) : sin(theta_e+/2).
C*       CHLF2    : (R*8) : cos(theta_e+/2).
C*       FI2      : (R*8) : phi_e+
C*       NHEB2    : (I*4) : beam e+ helicity.
C*       NHE2     : (I*4) : final-state e+ helicity.
C*       SPEM(*)  : (C*16) : beam e- spinor.
C*       SPEP(*)  : (C*16) : beam e+ spinor.
C*       SP1 (*)  : (C*16) : final-state e- spinor.
C*       SP2 (*)  : (C*16) : final-state e+ spinor.
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
C*    95/04/26  K.Fujii		Original version.
C*
CC**********************************************************************
 

      SUBROUTINE AMPEEWW(GAL,GZL,GZN,GWF,GWWA,GWWZ, GW, GWWH, GZZH,
     .                       AMW,GMW, AMZ,GMZ, AME, AMH,GMH,   EBM,
     . 		                   E1, SHLF1,CHLF1,FI1, NHEB1,NHE1,
     .		                   E2, SHLF2,CHLF2,FI2, NHEB2,NHE2,
     .                                 SPEM,SPEP, SP1,SP2,VCT, AMP)
 
      IMPLICIT   REAL*8  ( A-H, O-Z )
      REAL   *8  GAL(2), GZL(2), GZN(2), GWF(2), GWWA, GWWZ, 
     .           GW, GWWH, GZZH,
     .           AMW, GMW, AMZ, GMZ, AME, AMH, GMH, EBM,
     .		 E1, SHLF1, CHLF1, FI1,  E2, SHLF2, CHLF2, FI2    
      INTEGER*4  NHEB1, NHE1, NHEB2, NHE2
      COMPLEX*16 SPEM(6), SPEP(6), SP1(6), SP2(6), VCT(6,2), 
     .		 AMP(0:3)
C--
      COMPLEX*16 WRK(6,40), TMP(6), ZMP(18)
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
         FACT  = ABS(GAL(1))/SQRT(4*PI/128.D0)
         AMN   = 0
         GMN   = 0
         GME   = 0
      ENDIF
C--
C  Calculate internal lines.
C--
      CALL JEEXXX(EBM,E1,SHLF1,CHLF1,FI1,NHEB1,NHE1,+1, WRK(1, 1))
      CALL JIOXXX(SPEM,SP1,GZL,AMZ,GMZ,                 WRK(1, 2))
         WRK(1,1) = FACT*WRK(1,1)
         WRK(2,1) = FACT*WRK(2,1)
         WRK(3,1) = FACT*WRK(3,1)
         WRK(4,1) = FACT*WRK(4,1)
C--
      CALL JEEXXX(EBM,E2,SHLF2,CHLF2,FI2,NHEB2,NHE2,-1, WRK(1, 3))
      CALL JIOXXX(SP2,SPEP,GZL,AMZ,GMZ,                 WRK(1, 4))
         WRK(1,3) = FACT*WRK(1,3)
         WRK(2,3) = FACT*WRK(2,3)
         WRK(3,3) = FACT*WRK(3,3)
         WRK(4,3) = FACT*WRK(4,3)
C--
      CALL FVOXXX(SP1,VCT(1,2), GWF,AMN,GMN, WRK(1, 5))
      CALL JIOXXX(SPEM,WRK(1,5),GWF,AMW,GMW, WRK(1, 6))
C--
      CALL FVIXXX(SP2,VCT(1,1), GWF,AMN,GMN, WRK(1, 7))
      CALL JIOXXX(WRK(1,7),SPEP,GWF,AMW,GMW, WRK(1, 8))
C--
      CALL FVIXXX(SPEM,VCT(1,1),GWF,AMN,GMN, WRK(1, 9))
      CALL JIOXXX(WRK(1,9),SP1, GWF,AMW,GMW, WRK(1,10))
C--
      CALL FVOXXX(SPEP,VCT(1,2),GWF,AMN,GMN, WRK(1,11))
      CALL JIOXXX(SP2,WRK(1,11),GWF,AMW,GMW, WRK(1,12))
C--
      CALL FVIXXX(SPEM,WRK(1,3),GAL,AME,GME, WRK(1,13))
      CALL FVIXXX(SPEM,WRK(1,4),GZL,AME,GME, TMP      )
      WRK(1,13) = WRK(1,13) + TMP(1)
      WRK(2,13) = WRK(2,13) + TMP(2)
      WRK(3,13) = WRK(3,13) + TMP(3)
      WRK(4,13) = WRK(4,13) + TMP(4)
C--
      CALL FVOXXX(SPEP,WRK(1,1),GAL,AME,GME, WRK(1,14))
      CALL FVOXXX(SPEP,WRK(1,2),GZL,AME,GME, TMP      )
      WRK(1,14) = WRK(1,14) + TMP(1)
      WRK(2,14) = WRK(2,14) + TMP(2)
      WRK(3,14) = WRK(3,14) + TMP(3)
      WRK(4,14) = WRK(4,14) + TMP(4)
C--
      CALL FVOXXX(SP1,WRK(1,3), GAL,AME,GME, WRK(1,15))
      CALL FVOXXX(SP1,WRK(1,4), GZL,AME,GME, TMP      )
      WRK(1,15) = WRK(1,15) + TMP(1)
      WRK(2,15) = WRK(2,15) + TMP(2)
      WRK(3,15) = WRK(3,15) + TMP(3)
      WRK(4,15) = WRK(4,15) + TMP(4)
C--
      CALL FVIXXX(SP2,WRK(1,1), GAL,AME,GME, WRK(1,16))
      CALL FVIXXX(SP2,WRK(1,2), GZL,AME,GME, TMP      )
      WRK(1,16) = WRK(1,16) + TMP(1)
      WRK(2,16) = WRK(2,16) + TMP(2)
      WRK(3,16) = WRK(3,16) + TMP(3)
      WRK(4,16) = WRK(4,16) + TMP(4)
C--
      WRK(1,17) = WRK(1,1)*(GWWA/GW) + WRK(1,2)*(GWWZ/GW)
      WRK(2,17) = WRK(2,1)*(GWWA/GW) + WRK(2,2)*(GWWZ/GW)
      WRK(3,17) = WRK(3,1)*(GWWA/GW) + WRK(3,2)*(GWWZ/GW)
      WRK(4,17) = WRK(4,1)*(GWWA/GW) + WRK(4,2)*(GWWZ/GW)
      WRK(5,17) = WRK(5,1)
      WRK(6,17) = WRK(6,1)
C--
      WRK(1,18) = WRK(1,3)*(GWWA/GW) + WRK(1,4)*(GWWZ/GW)
      WRK(2,18) = WRK(2,3)*(GWWA/GW) + WRK(2,4)*(GWWZ/GW)
      WRK(3,18) = WRK(3,3)*(GWWA/GW) + WRK(3,4)*(GWWZ/GW)
      WRK(4,18) = WRK(4,3)*(GWWA/GW) + WRK(4,4)*(GWWZ/GW)
      WRK(5,18) = WRK(5,3)
      WRK(6,18) = WRK(6,3)
C--
C  Now calculate amplitudes.
C    The amplitude numbers given here correspond to those on P86 of
C    H.Murayama's Ph.D thesis.
C--
C  Non-fusion.
C--
C( 1)
      CALL IOVXXX(WRK(1, 7),SPEP,     WRK(1, 6),GWF, ZMP( 1))
C( 2)
      CALL IOVXXX(WRK(1, 7),WRK(1,11),WRK(1, 2),GZN, ZMP( 2))
C( 3)
      CALL IOVXXX(WRK(1, 9),WRK(1, 5),WRK(1, 4),GZN, ZMP( 3))
C( 4)
      CALL IOVXXX(SP2,      WRK(1,11),WRK(1,10),GWF, ZMP( 4))
C( 5)
      CALL J3XXXX(WRK(1,13),SP1, GAL,GZL,AMZ,GMZ, TMP)
      CALL VVVXXX(VCT(1, 1),VCT(1, 2),TMP      ,GW,  ZMP( 5))
C( 6)
      CALL J3XXXX(SPEM,WRK(1,15),GAL,GZL,AMZ,GMZ, TMP)
      CALL VVVXXX(VCT(1, 1),VCT(1, 2),TMP      ,GW,  ZMP( 6))
C( 7)
      CALL J3XXXX(WRK(1,16),SPEP,GAL,GZL,AMZ,GMZ, TMP)
      CALL VVVXXX(VCT(1, 1),VCT(1, 2),TMP      ,GW,  ZMP( 7))
C( 8)
      CALL J3XXXX(SP2, WRK(1,14),GAL,GZL,AMZ,GMZ, TMP)
      CALL VVVXXX(VCT(1, 1),VCT(1, 2),TMP      ,GW,  ZMP( 8))
C( 9)
      CALL IOVXXX(WRK(1,13),WRK(1, 5),VCT(1, 1),GWF, ZMP( 9))
C(10)
      CALL IOVXXX(WRK(1, 9),WRK(1,15),VCT(1, 2),GWF, ZMP(10))
C(11)
      CALL IOVXXX(WRK(1,16),WRK(1,11),VCT(1, 1),GWF, ZMP(11))
C(12)
      CALL IOVXXX(WRK(1, 7),WRK(1,14),VCT(1, 2),GWF, ZMP(12))
C(13)
      CALL VVVXXX(VCT(1, 1),WRK(1, 6),WRK(1,18),GW,  ZMP(13))
C(14)
      CALL VVVXXX(VCT(1, 1),WRK(1,12),WRK(1,17),GW,  ZMP(14))
C(15)
      CALL VVVXXX(WRK(1,10),VCT(1, 2),WRK(1,18),GW,  ZMP(15))
C(16)
      CALL VVVXXX(WRK(1, 8),VCT(1, 2),WRK(1,17),GW,  ZMP(16))
C--
C  Fusion.
C--
C(17-19)
      CALL W3W3XX(VCT(1, 1),WRK(1,17),VCT(1, 2),WRK(1,18),
     .            		      GW,GW,AMW,GMW, ZMP(17))
C--
C  Higgs.
C--
C(20)
      CALL HVVXXX(VCT(1, 1),VCT(1, 2),GWWH,AMH,GMH, TMP)
      CALL VVSXXX(WRK(1, 2),WRK(1, 4),TMP,GZZH,      ZMP(18))
C--
C  Sum up amplitudes.
C--
      AMP(1) = 0
      DO 10 I = 1, 16
         AMP(1) = AMP(1) + ZMP(I)
10    CONTINUE
      AMP(2) = ZMP(17)
      AMP(3) = ZMP(18) 
C--
      AMP(0) = AMP(1) + AMP(2) + AMP(3)
C--
C  That's it.
C--
      RETURN
      END
