CC**********************************************************************
C*
C*====================================================
C* Subroutine AMPWWZ(GAL,GZN,GZL,GWF,GZH,GWH,GW,GWWZ,
C*                   AME,AMW,GMW,AMZ,GMZ,AMH,GMH,
C*                   VCT,EIN,EOT, AMPSUM)
C*================================------==============
C*
C* (Purpose)
C*    Calculate amplitude for e+ + e- ---> W+ + W- + Z.
C*    Inputs are spinors and wave function thus replaceable with
C*    spinors or wave function made from final stable particle
C*    4-momenta.
C* (Inputs)
C*       GAL(*)  : (R*4) : lepton photon coupling.
C*       GZN(*)  : (R*4) : neutrino Z coupling.
C*       GZL(*)  : (R*4) : lepton Z coupling.
C*       GWF(*)  : (R*4) : fermion W coupling.
C*       GZH     : (R*4) : ZZ-H coupling.
C*       GWH     : (R*4) : WW-H coupling.
C*       GW      : (R*4) : weak coupling.
C*       GWWZ    : (R*4) : WW-Z coupling.
C*       AME     : (R*4) : electron mass.
C*       AMW     : (R*4) : W mass.
C*       GMW     : (R*4) : W width.
C*       AMZ     : (R*4) : Z mass.
C*       GMZ     : (R*4) : Z width.
C*       AMH     : (R*4) : Higss mass.
C*       GMH     : (R*4) : Higss width.
C*       VCT(*,i): (C*8) : i = (1,2,3) = (W-,W+,Z).
C*       EIN(*)  : (C*8) : incoming electron spinor.
C*       EOT(*)  : (C*8) : outgoing electron spinor.
C* (Output)
C*    AMPSUM(0)  : (C*8) : amplitude sum.
C*          (1)  : (C*8) : sum of non-Higgs diagrams.
C*          (2)  : (C*8) : sum of Higgs diagram.
C* (Relation)
C*    Invokes subroutines in HELAS.LOAD.
C* (Update Record)
C*    90/08/07  K.Fujii         Original version.
C*    95/03/12  K.Fujii		Modiifed to HELAS V204.
C*
CC**********************************************************************
 
      SUBROUTINE AMPWWZ(GAL,GZN,GZL,GWF,GZH,GWH,GW,GWWZ,
#ifdef ANOM_HVV
     .                  ALAM,A,B,BTILDE,
#endif
     .                  AME,AMW,GMW,AMZ,GMZ,AMH,GMH,
     .                  VCT,EIN,EOT, AMPSUM)
 
      IMPLICIT   REAL*4  ( A-H, O-Z )
      REAL   *4  GAL(2), GZN(2), GZL(2), GWF(2), GZH, GWH,
     .           GW, GWWZ,
#ifdef ANOM_HVV
     .           ALAM, A, B, BTILDE
#endif
     .           AME, AMW, GMW, AMZ, GMZ, AMH, GMH
      COMPLEX*8  EIN(6), EOT(6), VCT(6,3), AMPSUM(0:2)
C--
      COMPLEX*8  WRK(6,20), AMP(9)
C
C========< Entry Point >================================================
C
C--
C  Zero electron width.
C--
      GME = 0
      AMN = 0
      GMN = 0
C--
C  Calculate diagrm (1).
C--
      CALL J3XXXX(EIN,EOT,GAL,GZL,AMZ,GMZ,WRK(1,1))
      CALL W3W3XX(VCT(1,1),WRK(1,1),VCT(1,2),VCT(1,3),
     .            GW,GWWZ,AMW,GMW,AMP(1))
C--
C  Calculate diagrm (2).
C--
      CALL JVVXXX(VCT(1,3),VCT(1,1),GWWZ,AMW,GMW,WRK(1,2))
      CALL FVOXXX(EOT,VCT(1,2),GWF,AMN,GMN,WRK(1,3))
      CALL IOVXXX(EIN,WRK(1,3),WRK(1,2),GWF,AMP(2))
C--
C  Calculate diagrm (3).
C--
      CALL FVIXXX(EIN,VCT(1,1),GWF,AMN,GMN,WRK(1,4))
      CALL JVVXXX(VCT(1,2),VCT(1,3),GWWZ,AMW,GMW,WRK(1,5))
      CALL IOVXXX(WRK(1,4),EOT,WRK(1,5),GWF,AMP(3))
C--
C  Calculate diagrm (4).
C--
      CALL FVIXXX(EIN,VCT(1,3),GZL,AME,GME,WRK(1,6))
      CALL J3XXXX(WRK(1,6),EOT,GAL,GZL,AMZ,GMZ,WRK(1,7))
      CALL VVVXXX(VCT(1,1),VCT(1,2),WRK(1,7),GW,AMP(4))
C--
C  Calculate diagrm (5).
C--
      CALL FVOXXX(EOT,VCT(1,3),GZL,AME,GME,WRK(1,8))
      CALL J3XXXX(EIN,WRK(1,8),GAL,GZL,AMZ,GMZ,WRK(1,9))
      CALL VVVXXX(VCT(1,1),VCT(1,2),WRK(1,9),GW,AMP(5))
C--
C  Calculate diagrm (6).
C--
      CALL IOVXXX(WRK(1,4),WRK(1,8),VCT(1,2),GWF,AMP(6))
C--
C  Calculate diagrm (7).
C--
      CALL IOVXXX(WRK(1,4),WRK(1,3),VCT(1,3),GZN,AMP(7))
C--
C  Calculate diagrm (8).
C--
      CALL IOVXXX(WRK(1,6),WRK(1,3),VCT(1,1),GWF,AMP(8))
C--
C  Calculate diagrm (9).
C--
      CALL JIOXXX(EIN,EOT,GZL,AMZ,GMZ,WRK(1,10))
#ifndef ANOM_HVV
      CALL HVVXXX(VCT(1,1),VCT(1,2),GWH,AMH,GMH,WRK(1,11))
#else
      G1 = GWH + 2*AMW**2*A/ALAM
      G2 = -2*B/ALAM
      G3 = -4*BTILDE/ALAM
      CALL HVVXXA(VCT(1,1),VCT(1,2),G1,G2,G3,AMH,GMH,WRK(1,11))
#endif
      CALL VVSXXX(VCT(1,3),WRK(1,10),WRK(1,11),GZH,AMP(9))
C--
C  Sum non-Higgs and Higgs diagrams.
C--
      AMPSUM(1) = 0
      DO 10 I = 1, 8
         AMPSUM(1) = AMPSUM(1) + AMP(I)
10    CONTINUE
      AMPSUM(2) = AMP(9)
      AMPSUM(0) = AMPSUM(1) + AMPSUM(2)
C--
C  That's it.
C--
      RETURN
      END
