CC**********************************************************************
C*
C*===================================================================
C* Subroutine AMPEETT(GAL,GZL,GZN,GAU,GZU,GWF,GWWA,GWWZ,GW,GZZH,GTH,
C*                    AMW,GMW, AMZ,GMZ, AME,AMT,GMT, AMH,GMH,   EBM,
C*  		                    E1, SHLF1,CHLF1,FI1, NHEB1,NHE1,
C* 		                    E2, SHLF2,CHLF2,FI2, NHEB2,NHE2,
C*                                  SPEM,SPEP, SP1,SP2,ST1,ST2, AMP)
C*==============================================================---==
C*
C* (Purpose)
C*    Calculate t-channel amplitude for e+ e- ---> e+ e- t tbar.
C*    Inputs are spinors and wave functions thus replaceable with
C*    spinors or wave functions made from final stable particles.
C* (Inputs)
C*       GAL(*)   : (R*8)  : lepton photon coupling.
C*       GZL(*)   : (R*8)  : lepton Z coupling.
C*       GZN(*)   : (R*8)  : neutrino Z coupling.
C*       GAU(*)   : (R*8)  : up-type quark photon coupling.
C*       GZU(*)   : (R*8)  : up-type quark Z coupling.
C*       GWF      : (R*8)  : weak coupling.
C*       GWWA     : (R*8)  : WWA coupling.
C*       GWWZ     : (R*8)  : WWZ coupling.
C*       GW       : (R*8)  : gw.
C*       GZZH     : (R*8)  : ZZH coupling.
C*       GTH(*)   : (R*16) : ZZH coupling.
C*       AMW      : (R*8)  : W mass.
C*       GMW      : (R*8)  : W width.
C*       AMZ      : (R*8)  : Z mass.
C*       GMZ      : (R*8)  : Z width.
C*       AME      : (R*8)  : electron mass.
C*       AMT      : (R*8)  : t mass.
C*       GMT      : (R*8)  : t width.
C*       AMH      : (R*8)  : H mass.
C*       GMH      : (R*8)  : H width.
C*       EBM      : (R*8)  : beam energy.
C*       E1       : (R*8)  : final-state e- energy.
C*       SHLF1    : (R*8)  : sin(theta_e-/2).
C*       CHLF1    : (R*8)  : cos(theta_e-/2).
C*       FI1      : (R*8)  : phi_e-
C*       NHEB1    : (I*4)  : beam e- helicity.
C*       NHE1     : (I*4)  : final-state e- helicity.
C*       E2       : (R*8)  : final-state e+ energy.
C*       SHLF2    : (R*8)  : sin(theta_e+/2).
C*       CHLF2    : (R*8)  : cos(theta_e+/2).
C*       FI2      : (R*8)  : phi_e+
C*       NHEB2    : (I*4)  : beam e+ helicity.
C*       NHE2     : (I*4)  : final-state e+ helicity.
C*       SPEM(*)  : (C*16) : beam e- spinor.
C*       SPEP(*)  : (C*16) : beam e+ spinor.
C*       SP1 (*)  : (C*16) : final-state e- spinor.
C*       SP2 (*)  : (C*16) : final-state e+ spinor.
C*       ST1(*,1) : (C*16) : final-state t spinor.
C*       ST2(*,2) : (C*16) : final-state tbar spinor.
C* (Output)
C*       AMP(0)   : (C*16) : total t-channel amplitudes.
C*          (1)   : (C*16) : non-fusion diagrams.
C*          (2)   : (C*16) : fusion diagrams.
C*          (3)   : (C*16) : Higgs diagram.
C* (Relation)
C*    Invokes subroutines in HELAS.LOAD.
C* (Update Record)
C*    95/06/12  K.Fujii		Original version.
C*
CC**********************************************************************
 
      SUBROUTINE AMPEETT(GAL,GZL,GZN,GAU,GZU,GWF,GWWA,GWWZ,GW,GZZH,GTH,
     .			 AMW,GMW, AMZ,GMZ, AME,AMT,GMT, AMH,GMH,   EBM,
     .				       E1, SHLF1,CHLF1,FI1, NHEB1,NHE1,
     .				       E2, SHLF2,CHLF2,FI2, NHEB2,NHE2,
     .				       SPEM,SPEP, SP1,SP2,ST1,ST2, AMP)
 
      IMPLICIT   REAL*8  ( A-H, O-Z )
      REAL   *8  GAL(2), GZL(2), GZN(2), GAU(2), GZU(2), GWF(2), 
     .           GWWA, GWWZ, GW, GZZH,
     .           AMW, GMW, AMZ, GMZ, AME, AMT, GMT, AMH, GMH, EBM,
     .		 E1, SHLF1, CHLF1, FI1,  E2, SHLF2, CHLF2, FI2    
      INTEGER*4  NHEB1, NHE1, NHEB2, NHE2
      COMPLEX*16 SPEM(6), SPEP(6), SP1(6), SP2(6), ST1(6), ST2(6),
     .		 GTH(2), AMP(0:3)
C--
      COMPLEX*16 WRK(6,6), TMP1(6), TMP2(6), ZMP(7), ZTMP
C--
      DATA AMHMX /10000.D0/
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
         AMG   = 0
         GMG   = 0
         GME   = 0
      ENDIF
C--
C  Calculate internal lines.
C--
      CALL JEEXXX(EBM,E1,SHLF1,CHLF1,FI1,NHEB1,NHE1,+1, WRK(1, 1))
      CALL JIOXXX(SPEM,SP1,GZL,AMZ,GMZ,			WRK(1, 2))
         WRK(1,1) = FACT*WRK(1,1)
         WRK(2,1) = FACT*WRK(2,1)
         WRK(3,1) = FACT*WRK(3,1)
         WRK(4,1) = FACT*WRK(4,1)
C--
      CALL JEEXXX(EBM,E2,SHLF2,CHLF2,FI2,NHEB2,NHE2,-1, WRK(1, 3))
      CALL JIOXXX(SP2,SPEP,GZL,AMZ,GMZ,			WRK(1, 4))
         WRK(1,3) = FACT*WRK(1,3)
         WRK(2,3) = FACT*WRK(2,3)
         WRK(3,3) = FACT*WRK(3,3)
         WRK(4,3) = FACT*WRK(4,3)
C--
      CALL JIOXXX(ST2,ST1,GAU,AMG,GMG, WRK(1,5))
      CALL JIOXXX(ST2,ST1,GZU,AMZ,GMZ, WRK(1,6))
C--
C  Now calculate amplitudes.
C    The amplitude numbers given here correspond to those on P87 of
C    H.Murayama's Ph.D thesis.
C--
C  Non-fusion.
C--
C (1)
      CALL FVOXXX(SP1,WRK(1,5),GAL,AME,GME,TMP1)
      CALL FVOXXX(SP1,WRK(1,6),GZL,AME,GME,TMP2)
         TMP1(1) = TMP1(1) + TMP2(1)
         TMP1(2) = TMP1(2) + TMP2(2)
         TMP1(3) = TMP1(3) + TMP2(3)
         TMP1(4) = TMP1(4) + TMP2(4)
      CALL IOVXXX(SPEM,TMP1,WRK(1,3),GAL,ZTMP)
      CALL IOVXXX(SPEM,TMP1,WRK(1,4),GZL,ZMP(1))
      ZMP(1) = ZMP(1) + ZTMP
C (2)
      CALL FVIXXX(SPEM,WRK(1,5),GAL,AME,GME,TMP1)
      CALL FVIXXX(SPEM,WRK(1,6),GZL,AME,GME,TMP2)
         TMP1(1) = TMP1(1) + TMP2(1)
         TMP1(2) = TMP1(2) + TMP2(2)
         TMP1(3) = TMP1(3) + TMP2(3)
         TMP1(4) = TMP1(4) + TMP2(4)
      CALL IOVXXX(TMP1,SP1,WRK(1,3),GAL,ZTMP)
      CALL IOVXXX(TMP1,SP1,WRK(1,4),GZL,ZMP(2))
      ZMP(2) = ZMP(2) + ZTMP
C (3)
      CALL FVIXXX(SP2,WRK(1,5),GAL,AME,GME,TMP1)
      CALL FVIXXX(SP2,WRK(1,6),GZL,AME,GME,TMP2)
         TMP1(1) = TMP1(1) + TMP2(1)
         TMP1(2) = TMP1(2) + TMP2(2)
         TMP1(3) = TMP1(3) + TMP2(3)
         TMP1(4) = TMP1(4) + TMP2(4)
      CALL IOVXXX(TMP1,SPEP,WRK(1,1),GAL,ZTMP)
      CALL IOVXXX(TMP1,SPEP,WRK(1,2),GZL,ZMP(3))
      ZMP(3) = ZMP(3) + ZTMP
C (4)
      CALL FVOXXX(SPEP,WRK(1,5),GAL,AME,GME,TMP1)
      CALL FVOXXX(SPEP,WRK(1,6),GZL,AME,GME,TMP2)
         TMP1(1) = TMP1(1) + TMP2(1)
         TMP1(2) = TMP1(2) + TMP2(2)
         TMP1(3) = TMP1(3) + TMP2(3)
         TMP1(4) = TMP1(4) + TMP2(4)
      CALL IOVXXX(SP2,TMP1,WRK(1,1),GAL,ZTMP)
      CALL IOVXXX(SP2,TMP1,WRK(1,2),GZL,ZMP(4))
      ZMP(4) = ZMP(4) + ZTMP
C--
C  Fusion.
C--
C (5) 
      CALL FVIXXX(ST2,WRK(1,3),GAU,AMT,GMT,TMP1)
      CALL FVIXXX(ST2,WRK(1,4),GZU,AMT,GMT,TMP2)
         TMP1(1) = TMP1(1) + TMP2(1)
         TMP1(2) = TMP1(2) + TMP2(2)
         TMP1(3) = TMP1(3) + TMP2(3)
         TMP1(4) = TMP1(4) + TMP2(4)
      CALL IOVXXX(TMP1, ST1,WRK(1,1),GAU,ZTMP)
      CALL IOVXXX(TMP1, ST1,WRK(1,2),GZU,ZMP(5))
      ZMP(5) = ZMP(5) + ZTMP
C (6) 
      CALL FVIXXX(ST2,WRK(1,1),GAU,AMT,GMT,TMP1)
      CALL FVIXXX(ST2,WRK(1,2),GZU,AMT,GMT,TMP2)
         TMP1(1) = TMP1(1) + TMP2(1)
         TMP1(2) = TMP1(2) + TMP2(2)
         TMP1(3) = TMP1(3) + TMP2(3)
         TMP1(4) = TMP1(4) + TMP2(4)
      CALL IOVXXX(TMP1, ST1,WRK(1,3),GAU,ZTMP)
      CALL IOVXXX(TMP1, ST1,WRK(1,4),GZU,ZMP(6))
      ZMP(6) = ZMP(6) + ZTMP
C--
C  Higgs.
C--
C (7)
      IF ( AMH.LT.AMHMX ) THEN
         CALL HIOXXX(ST2,ST1,GTH,AMH,GMH,TMP1)
         CALL VVSXXX(WRK(1,2),WRK(1,4),TMP1,GZZH, ZMP(7))
      ELSE
         ZMP(7) = 0
      ENDIF
C--
C  Sum up amplitudes.
C--
      AMP(1) = 0
      DO 10 I = 1, 4
         AMP(1) = AMP(1) + ZMP(I)
10    CONTINUE
      AMP(2) = ZMP(5) + ZMP(6)
      AMP(3) = ZMP(7) 
C--
      AMP(0) = AMP(1) + AMP(2) + AMP(3)
C--
C  That's it.
C--
      RETURN
      END
