CC**********************************************************************
C*
C*============================================================
C* Subroutine AMPEEZZ(GAL,GZL,GZZH, AMZ,GMZ,AME,AMH,GMH, EBM,
C*  		             E1, SHLF1,CHLF1,FI1, NHEB1,NHE1,
C* 		             E2, SHLF2,CHLF2,FI2, NHEB2,NHE2,
C*                               SPEM,SPEP, SP1,SP2,VCT, AMP)
C*=======================================================---==
C*
C* (Purpose)
C*    Calculate t-channel amplitude for e+ e- ---> e+ e- Z0 Z0.
C*    Inputs are spinors and wave functions thus replaceable with
C*    spinors or wave functions made from final stable particles.
C* (Inputs)
C*       GAL(*)   : (R*8) : lepton photon coupling.
C*       GZL(*)   : (R*8) : lepton Z coupling.
C*       GZZH     : (R*8) : ZZH coupling.
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
C*    05/03/06  K.Fujii		Original version.
C*
CC**********************************************************************
 

      SUBROUTINE AMPEEZZ(GAL,GZL,GZZH, AMZ,GMZ,AME,AMH,GMH,   EBM,
     .                             E1, SHLF1,CHLF1,FI1, NHEB1,NHE1,
     .                             E2, SHLF2,CHLF2,FI2, NHEB2,NHE2,
     .                                 SPEM,SPEP, SP1,SP2,VCT, AMP)
 
      IMPLICIT   REAL*8  ( A-H, O-Z )
      REAL   *8  GAL(2), GZL(2), GZZH, AMZ, GMZ, AME, AMH, GMH, EBM,
     .           E1, SHLF1, CHLF1, FI1,  E2, SHLF2, CHLF2, FI2    
      INTEGER*4  NHEB1, NHE1, NHEB2, NHE2
      COMPLEX*16 SPEM(6), SPEP(6), SP1(6), SP2(6), VCT(6,2), 
     .           AMP(0:3)
C--
      COMPLEX*16 WRK(6,20), TMP(6), ZMP(43)
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
         GME   = 0
         AMA   = 0
         GMA   = 0
      ENDIF
C--
C  Calculate internal lines.
C--
C-- virtual incoming e- with 1 Z
      CALL FVIXXX(SPEM     ,VCT(1,1),GZL,AME,GME, WRK(1, 1))
      CALL FVIXXX(SPEM     ,VCT(1,2),GZL,AME,GME, WRK(1, 2))
C-- virtual incoming e- with 2 Z's
      CALL FVIXXX(WRK(1, 1),VCT(1,2),GZL,AME,GME, WRK(1, 3))
      CALL FVIXXX(WRK(1, 2),VCT(1,1),GZL,AME,GME, WRK(1, 4))
C-- virtual outgoing e+ with 1 Z
      CALL FVOXXX(SPEP     ,VCT(1,1),GZL,AME,GME, WRK(1, 5))
      CALL FVOXXX(SPEP     ,VCT(1,2),GZL,AME,GME, WRK(1, 6))
C-- virtual outgoing e+ with 2 Z's
      CALL FVOXXX(WRK(1, 5),VCT(1,2),GZL,AME,GME, WRK(1, 7))
      CALL FVOXXX(WRK(1, 6),VCT(1,1),GZL,AME,GME, WRK(1, 8))
C-- virtual outgoing e- with 1 Z
      CALL FVOXXX(SP1      ,VCT(1,1),GZL,AME,GME, WRK(1, 9))
      CALL FVOXXX(SP1      ,VCT(1,2),GZL,AME,GME, WRK(1,10))
C-- virtual outgoing e- with 2 Z's
      CALL FVOXXX(WRK(1, 9),VCT(1,2),GZL,AME,GME, WRK(1,11))
      CALL FVOXXX(WRK(1,10),VCT(1,1),GZL,AME,GME, WRK(1,12))
C-- virtual incoming e+  with 1 Z
      CALL FVIXXX(SP2      ,VCT(1,1),GZL,AME,GME, WRK(1,13))
      CALL FVIXXX(SP2      ,VCT(1,2),GZL,AME,GME, WRK(1,14))
C-- virtual incoming e+ with 2 Z's
      CALL FVIXXX(WRK(1,13),VCT(1,2),GZL,AME,GME, WRK(1,15))
      CALL FVIXXX(WRK(1,14),VCT(1,1),GZL,AME,GME, WRK(1,16))
C-- virtual almost real A/Z
C     (e- leg)
      CALL JEEXXX(EBM,E1,SHLF1,CHLF1,FI1,NHEB1,NHE1,+1, WRK(1,17))
      CALL JIOXXX(SPEM,SP1,GZL,AMZ,GMZ,                 WRK(1,18))
      WRK(1,17) = FACT*WRK(1,17)
      WRK(2,17) = FACT*WRK(2,17)
      WRK(3,17) = FACT*WRK(3,17)
      WRK(4,17) = FACT*WRK(4,17)
C--   (e+ leg)
      CALL JEEXXX(EBM,E2,SHLF2,CHLF2,FI2,NHEB2,NHE2,-1, WRK(1,19))
      CALL JIOXXX(SP2,SPEP,GZL,AMZ,GMZ,                 WRK(1,20))
      WRK(1,19) = FACT*WRK(1,19)
      WRK(2,19) = FACT*WRK(2,19)
      WRK(3,19) = FACT*WRK(3,19)
      WRK(4,19) = FACT*WRK(4,19)
C--
C  Now calculate amplitudes.
C--
C  Non-fusion.
C--
C-- 2 Z's from initial e- leg
C     (A exchange)
      CALL IOVXXX(WRK(1, 3),SP1,  WRK(1,19), GAL, ZMP( 1)) 
      CALL IOVXXX(WRK(1, 4),SP1,  WRK(1,19), GAL, ZMP( 2)) 
C     (Z exchange)
      CALL IOVXXX(WRK(1, 3),SP1,  WRK(1,20), GZL, ZMP( 3)) 
      CALL IOVXXX(WRK(1, 4),SP1,  WRK(1,20), GZL, ZMP( 4)) 

C-- 1 Z from initial e-, 1 Z from initial e+
C     (A exchange)
      CALL JIOXXX(WRK(1, 1), SP1, GAL, AMA,GMA, TMP)
      CALL IOVXXX(SP2, WRK(1, 6),       TMP, GAL, ZMP( 5))
      CALL JIOXXX(WRK(1, 2), SP1, GAL, AMA,GMA, TMP)
      CALL IOVXXX(SP2, WRK(1, 5),       TMP, GAL, ZMP( 6))
C     (Z exchange)
      CALL JIOXXX(WRK(1, 1), SP1, GZL, AMZ,GMZ, TMP)
      CALL IOVXXX(SP2, WRK(1, 6),       TMP, GZL, ZMP( 7))
      CALL JIOXXX(WRK(1, 2), SP1, GZL, AMZ,GMZ, TMP)
      CALL IOVXXX(SP2, WRK(1, 5),       TMP, GZL, ZMP( 8))

C-- 1 Z from initial e-, 1 Z from final e-
C     (A exchange)
      CALL IOVXXX(WRK(1, 1),WRK(1,10),WRK(1,19), GAL, ZMP( 9))
      CALL IOVXXX(WRK(1, 2),WRK(1, 9),WRK(1,19), GAL, ZMP(10))
C     (Z exchange)
      CALL IOVXXX(WRK(1, 1),WRK(1,10),WRK(1,20), GZL, ZMP(11))
      CALL IOVXXX(WRK(1, 2),WRK(1, 9),WRK(1,20), GZL, ZMP(12))

C-- 1 Z from initial e-, 1 Z from final e+
C     (A exchange)
      CALL JIOXXX(WRK(1, 1), SP1, GAL, AMA,GMA, TMP)
      CALL IOVXXX(WRK(1,14), SPEP,      TMP, GAL, ZMP(13))
      CALL JIOXXX(WRK(1, 2), SP1, GAL, AMA,GMA, TMP)
      CALL IOVXXX(WRK(1,13), SPEP,      TMP, GAL, ZMP(14))
C     (Z exchange)
      CALL JIOXXX(WRK(1, 1), SP1, GZL, AMZ,GMZ, TMP)
      CALL IOVXXX(WRK(1,14), SPEP,      TMP, GZL, ZMP(15))
      CALL JIOXXX(WRK(1, 2), SP1, GZL, AMZ,GMZ, TMP)
      CALL IOVXXX(WRK(1,13), SPEP,      TMP, GZL, ZMP(16))

C-- 2 Z's from initial e+ leg
C     (A exchange)
      CALL IOVXXX(SP2,WRK(1, 7),  WRK(1,17), GAL, ZMP(17)) 
      CALL IOVXXX(SP2,WRK(1, 8),  WRK(1,17), GAL, ZMP(18)) 
C     (Z exchange)
      CALL IOVXXX(SP2,WRK(1, 7),  WRK(1,18), GZL, ZMP(19)) 
      CALL IOVXXX(SP2,WRK(1, 8),  WRK(1,18), GZL, ZMP(20)) 

C-- 1 Z from initial e+, 1 Z from final e-
C     (A exchange)
      CALL JIOXXX(SP2, WRK(1, 5), GAL, AMA,GMA, TMP)
      CALL IOVXXX(SPEM,WRK(1,10),       TMP, GAL, ZMP(21))
      CALL JIOXXX(SP2, WRK(1, 6), GAL, AMA,GMA, TMP)
      CALL IOVXXX(SPEM,WRK(1, 9),       TMP, GAL, ZMP(22))
C     (Z exchange)
      CALL JIOXXX(SP2, WRK(1, 5), GZL, AMZ,GMZ, TMP)
      CALL IOVXXX(SPEM,WRK(1,10),       TMP, GZL, ZMP(23))
      CALL JIOXXX(SP2, WRK(1, 6), GZL, AMZ,GMZ, TMP)
      CALL IOVXXX(SPEM,WRK(1, 9),       TMP, GZL, ZMP(24))

C-- 1 Z from initial e+, 1 Z from final e+
C     (A exchange)
      CALL IOVXXX(WRK(1,14),WRK(1, 5),WRK(1,17), GAL, ZMP(25))
      CALL IOVXXX(WRK(1,13),WRK(1, 6),WRK(1,17), GAL, ZMP(26))
C     (Z exchange)
      CALL IOVXXX(WRK(1,14),WRK(1, 5),WRK(1,18), GZL, ZMP(27))
      CALL IOVXXX(WRK(1,13),WRK(1, 6),WRK(1,18), GZL, ZMP(28))

C-- 2 Z's from final e- leg
C     (A exchange)
      CALL IOVXXX(SPEM,WRK(1,11), WRK(1,19), GAL, ZMP(29)) 
      CALL IOVXXX(SPEM,WRK(1,12), WRK(1,19), GAL, ZMP(30)) 
C     (Z exchange)
      CALL IOVXXX(SPEM,WRK(1,11), WRK(1,20), GZL, ZMP(31)) 
      CALL IOVXXX(SPEM,WRK(1,12), WRK(1,20), GZL, ZMP(32)) 

C-- 1 Z from final e-, 1 Z from final e+
C     (A exchange)
      CALL JIOXXX(SPEM, WRK(1, 9), GAL, AMA,GMA, TMP)
      CALL IOVXXX(WRK(1,14), SPEP,      TMP, GAL, ZMP(33))
      CALL JIOXXX(SPEM, WRK(1,10), GAL, AMA,GMA, TMP)
      CALL IOVXXX(WRK(1,13), SPEP,      TMP, GAL, ZMP(34))
C     (Z exchange)
      CALL JIOXXX(SPEM, WRK(1, 9), GZL, AMZ,GMZ, TMP)
      CALL IOVXXX(WRK(1,14), SPEP,      TMP, GZL, ZMP(35))
      CALL JIOXXX(SPEM, WRK(1,10), GZL, AMZ,GMZ, TMP)
      CALL IOVXXX(WRK(1,13), SPEP,      TMP, GZL, ZMP(36))

C-- 2 Z's from final e+ leg
C     (A exchange)
      CALL IOVXXX(WRK(1,15), SPEP,WRK(1,17), GAL, ZMP(37)) 
      CALL IOVXXX(WRK(1,16), SPEP,WRK(1,17), GAL, ZMP(38)) 
C     (Z exchange)
      CALL IOVXXX(WRK(1,15), SPEP,WRK(1,18), GZL, ZMP(39)) 
      CALL IOVXXX(WRK(1,16), SPEP,WRK(1,18), GZL, ZMP(40)) 
C--
C  Higgs.
C--
C-- t-channel H
      CALL HVVXXX(WRK(1,18),VCT(1, 1),GZZH,AMH,GMH, TMP)
      CALL VVSXXX(VCT(1, 2),WRK(1,20),TMP,GZZH,   ZMP(41))
C-- u-channel H
      CALL HVVXXX(WRK(1,18),VCT(1, 2),GZZH,AMH,GMH, TMP)
      CALL VVSXXX(VCT(1, 1),WRK(1,20),TMP,GZZH,   ZMP(42))
C-- s-channel H
      CALL HVVXXX(VCT(1, 1),VCT(1, 2),GZZH,AMH,GMH, TMP)
      CALL VVSXXX(WRK(1,18),WRK(1,20),TMP,GZZH,   ZMP(43))
C--
C  Sum up amplitudes.
C--
      AMP(1) = 0
      DO 10 I = 1, 40
         AMP(1) = AMP(1) + ZMP(I)
10    CONTINUE
      AMP(2) = ZMP(41) + ZMP(42)
      AMP(3) = ZMP(43) 
C--
      AMP(0) = AMP(1) + AMP(2) + AMP(3)
C--
C  That's it.
C--
      RETURN
      END
