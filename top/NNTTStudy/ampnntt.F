CC**********************************************************************
C*
C*============================================================
C* Subroutine AMPNNTT(GAL,GZL,GZN,GAU,GZU,GWF, GW, GWWH, GTH,
C*                    AMW,GMW, AMZ,GMZ, AME,AMB,GMB, AMH,GMH,
C*                          SPEM,SPEP, SP1,SP2, ST1,ST2, AMP)
C*=======================================================---==
C*
C* (Purpose)
C*    Calculate t-channel amplitude for e+ e- ---> nu nub t tb.
C*    Inputs are spinors and thus replaceable with virtual
C*    spinors made from final stable particles.
C* (Inputs)
C*       GAL(*)   : (R*8)  : lepton photon coupling.
C*       GZL(*)   : (R*8)  : lepton Z coupling.
C*       GZN(*)   : (R*8)  : neutrino Z coupling.
C*       GAU(*)   : (R*8)  : up-type quark photon coupling.
C*       GZU(*)   : (R*8)  : up-type quark Z coupling.
C*       GWF      : (R*8)  : weak coupling.
C*       GW       : (R*8)  : gw.
C*       GWWH     : (R*8)  : WWH coupling.
C*       GTH(*)   : (R*16) : ZZH coupling.
C*       AMW      : (R*8)  : W mass.
C*       GMW      : (R*8)  : W width.
C*       AMZ      : (R*8)  : Z mass.
C*       GMZ      : (R*8)  : Z width.
C*       AME      : (R*8)  : electron mass.
C*       AMB      : (R*8)  : b mass.
C*       GMB      : (R*8)  : b width.
C*       AMH      : (R*8)  : H mass.
C*       GMH      : (R*8)  : H width.
C*       SPEM(*)  : (C*16) : beam e- spinor.
C*       SPEP(*)  : (C*16) : beam e+ spinor.
C*       SP1 (*)  : (C*16) : final-state nu spinor.
C*       SP2 (*)  : (C*16) : final-state nub spinor.
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
C*    95/06/09  K.Fujii		Original version.
C*
CC**********************************************************************
 

      SUBROUTINE AMPNNTT(GAL,GZL,GZN,GAU,GZU,GWF, GW, GWWH, GTH,
     .                   AMW,GMW, AMZ,GMZ, AME,AMB,GMB, AMH,GMH,
     .                         SPEM,SPEP, SP1,SP2, ST1,ST2, AMP)
 
      IMPLICIT   REAL*8  ( A-H, O-Z )
C--
      REAL   *8  GAL(2), GZL(2), GZN(2), GAU(2), GZU(2), GWF(2), GW,
     .           GWWH, 
     .           AMW, GMW, AMZ, GMZ, AME, AMB, GMB, AMH, GMH
      COMPLEX*16 SPEM(6), SPEP(6), SP1(6), SP2(6), ST1(6), ST2(6),
     .		 GTH(2), AMP(0:3)
C--
      COMPLEX*16 WRK(6,4), TMP(6), ZMP(7), ZTMP
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
         AMG   = 0
         GMG   = 0
         AMN   = 0
         GMN   = 0
         GME   = 0
       ENDIF
C--
C  Calculate internal lines.
C--
      CALL JIOXXX(SPEM,SP1,GWF,AMW,GMW, WRK(1, 1))
      CALL JIOXXX(SP2,SPEP,GWF,AMW,GMW, WRK(1, 2))
      CALL JIOXXX(ST2, ST1,GAU,AMG,GMG, WRK(1, 3))
      CALL JIOXXX(ST2, ST1,GZU,AMZ,GMZ, WRK(1, 4))
C--
C  Now calculate amplitudes.
C    The amplitude numbers given here correspond to those on P85 of
C    H.Murayama's Ph.D thesis.
C--
C  Non-fusion.
C--
C (1)
      CALL FVOXXX(SP1,WRK(1,4),GZN,AMN,GMN,TMP)
      CALL IOVXXX(SPEM,TMP,WRK(1,2),GWF, ZMP(1))
C (2)
      CALL FVOXXX(SP1,WRK(1,2),GWF,AME,GME,TMP)
      CALL IOVXXX(SPEM,TMP,WRK(1,3),GAL, ZTMP)
      CALL IOVXXX(SPEM,TMP,WRK(1,4),GZL, ZMP(2))
      ZMP(2) = ZMP(2) + ZTMP
C (3)
      CALL FVIXXX(SP2,WRK(1,4),GZN,AMN,GMN,TMP)
      CALL IOVXXX(TMP,SPEP,WRK(1,1),GWF, ZMP(3))
C (4)
      CALL FVIXXX(SP2,WRK(1,1),GWF,AME,GME,TMP)
      CALL IOVXXX(TMP,SPEP,WRK(1,3),GAL, ZTMP)
      CALL IOVXXX(TMP,SPEP,WRK(1,4),GZL, ZMP(4))
      ZMP(4) = ZMP(4) + ZTMP
C--
C  Fusion.
C--
C (5) 
      CALL FVIXXX(ST2,WRK(1,1),GWF,AMB,GMB,TMP)
      CALL IOVXXX(TMP, ST1,WRK(1,2),GWF, ZMP(5))
C (6)
      CALL J3XXXX(ST2,ST1,GAU,GZU,AMZ,GMZ,TMP)
      CALL VVVXXX(WRK(1,2),WRK(1,1),TMP,GW, ZMP(6))
C--
C  Higgs.
C--
C (7)
      IF ( AMH.LT.AMHMX ) THEN
         CALL HIOXXX(ST2,ST1,GTH,AMH,GMH,TMP)
         CALL VVSXXX(WRK(1,2),WRK(1,1),TMP,GWWH, ZMP(7))
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
