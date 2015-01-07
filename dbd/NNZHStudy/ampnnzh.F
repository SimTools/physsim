CC**********************************************************************
C*
C*===================================================================
C* Subroutine AMPNNZH(GZL,GZN,GWF,GWWZ,GZZH,GWWH, AMW,GMW,AMZ,GMZ, AME,
C*                                SPEM,SPEP,SPNU,SPNB,VCTZ,SCH,AMP)
C*=============================================================---===
C*
C* (Purpose)
C*    Calculate t-channel amplitude for e+ + e- ---> nu+nubar + Z H.
C*    Inputs are spinors and wave functions thus replaceable with
C*    spinors or wave function made from final stable particles.
C* (Inputs)
C*       GZL(*)   : (R*8) : lepton Z coupling.
C*       GZN(*)   : (R*8) : neutrino Z coupling.
C*       GWF      : (R*8) : weak coupling.
C*       GWWZ     : (R*8) : WWZ coupling.
C*       GZZH     : (R*8) : ZZH coupling.
C*       GWWH     : (R*8) : WWH coupling.
C*       AMW      : (R*8) : W mass.
C*       GMW      : (R*8) : W width.
C*       AMZ      : (R*8) : Z mass.
C*       GMZ      : (R*8) : Z width.
C*       AME      : (R*8) : electron mass.
C*       SPEM(*)  : (C*16) : beam electron spinor.
C*       SPEP(*)  : (C*16) : beam positron spinor.
C*       SPNU (*) : (C*16) : final-state neutrino spinor.
C*       SPNB (*) : (C*16) : final-state neutrino-bar spinor.
C*       VCTZ (*) : (C*16) : final-state Z wave function.
C*       SCH  (*) : (C*16) : final-state H wave function.
C* (Output)
C*       AMP(*)   : (C*16) : amplitudes for t-channel diagram.
C* (Relation)
C*    Invokes AMPNNZ, AMPNNH and subroutines in HELAS.LOAD.
C* (Update Record)
C*    2010/06/19  J.Tian & K.Fujii		Original version.
C*
CC**********************************************************************
 
      SUBROUTINE AMPNNZH(GZL,GZN,GWF,GWWZ,GZZH,GWWH,
     .                   AMW,GMW,AMZ,GMZ, AME,
     .                   SPEM,SPEP,SPNU,SPNB,VCTZ,SCH, AMP)
 
      IMPLICIT   REAL*8  ( A-H, O-Z )
      REAL   *8  GZL(2), GZN(2), GWF(2), GWWZ, GZZH, GWWH, 
     .           AMW, GMW, AMZ, GMZ, AME
      COMPLEX*16 SPEM(6), SPEP(6), SPNU(6), SPNB(6), VCTZ(6), SCH(3),
     .           AMP(0:7)
C--
      COMPLEX*16 WRK(6,10)
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
C  Amplitude (1) : nnZ* followd by Z* -> ZH
C--
      CALL JVSXXX(VCTZ,SCH,GZZH,AMZ,GMZ, WRK(1,1))
      CALL AMPNNZ(GZL,GZN,GWF,GWWZ,AMW,GMW,AME,
     .            SPEM,SPEP,SPNU,SPNB,WRK(1,1),AMP(1))
C--
C  Amplitude (2-5) : e+(*) e-(*) -> nnH with e+ or e- being virtual (eZ)
C--
      CALL FVIXXX(SPEM,VCTZ,GZL,AME,GME,WRK(1,2))       
      CALL AMPNNH(GWF,GWWH,AMW,GMW,AME,AMH,GMH,
     .            WRK(1,2),SPEP, SPNU,SPNB, SCH, AMP(2))

      CALL FVOXXX(SPEP,VCTZ,GZL,AME,GME,WRK(1,3))       
      CALL AMPNNH(GWF,GWWH,AMW,GMW,AME,AMH,GMH,
     .            SPEM,WRK(1,3), SPNU,SPNB, SCH, AMP(3))

      CALL FVOXXX(SPNU,VCTZ,GZN,AMN,GMN,WRK(1,4))       
      CALL AMPNNH(GWF,GWWH,AMW,GMW,AME,AMH,GMH,
     .            SPEM,SPEP, WRK(1,4),SPNB, SCH, AMP(4))

      CALL FVIXXX(SPNB,VCTZ,GZN,AMN,GMN,WRK(1,5))       
      CALL AMPNNH(GWF,GWWH,AMW,GMW,AME,AMH,GMH,
     .            SPEM,SPEP, SPNU,WRK(1,5), SCH, AMP(5))
C--
C  Amplitude (6-7).
C--
      CALL JIOXXX(SPEM,SPNU,GWF,AMW,GMW, WRK(1,6))
      CALL JIOXXX(SPNB,SPEP,GWF,AMW,GMW, WRK(1,7))

      CALL JVVXXX(WRK(1,6),VCTZ, GWWZ, AMW,GMW, WRK(1,8))
      CALL VVSXXX(WRK(1,7),WRK(1,8),SCH,GWWH, AMP(6))

      CALL JVVXXX(VCTZ,WRK(1,7), GWWZ, AMW,GMW, WRK(1,9))
      CALL VVSXXX(WRK(1,9),WRK(1,6),SCH,GWWH, AMP(7))
C--
C  Sum up amplitudes.
C--
      AMP(0) = AMP(1) + AMP(2) + AMP(3) + AMP(4) + AMP(5)
     .                + AMP(6) + AMP(7)
c      AMP(0) = AMP(1)
C--
C  That's it.
C--
      RETURN
      END
