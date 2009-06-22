CC**********************************************************************
C*
C*====================================================
C* Subroutine AMPZZZ(GAL,GZN,GZL,GWF,GZH,GWH,GW,GZZZ,
C*                   AME,AMW,GMW,AMZ,GMZ,AMH,GMH,
C*                   VCT,EIN,EOT, AMPSUM)
C*================================------==============
C*
C* (Purpose)
C*    Calculate amplitude for e+ + e- ---> Z2 + Z1 + Z.
C*    Inputs are spinors and wave function thus replaceable with
C*    spinors or wave function made from final stable particle
C*    4-momenta.
C* (Inputs)
C*       GAL(*)  : (R*4) : lepton photon coupling.
C*       GZN(*)  : (R*4) : neutrino Z coupling.
C*       GZL(*)  : (R*4) : lepton Z coupling.
C*       GWF(*)  : (R*4) : fermion W coupling.
C*       GZH     : (R*4) : ZZ-H coupling.
C*       GWH     : (R*4) : WZ1H coupling.
C*       GW      : (R*4) : weak coupling.
C*       GZZZ    : (R*4) : WZ1Z coupling.
C*       AME     : (R*4) : electron mass.
C*       AMW     : (R*4) : W mass.
C*       GMW     : (R*4) : W width.
C*       AMZ     : (R*4) : Z mass.
C*       GMZ     : (R*4) : Z width.
C*       AMH     : (R*4) : Higss mass.
C*       GMH     : (R*4) : Higss width.
C*       VCT(*,i): (C*8) : i = (1,2,3) = (Z1,Z2,Z).
C*       EIN(*)  : (C*8) : incoming electron spinor.
C*       EOT(*)  : (C*8) : outgoing electron spinor.
C* (Output)
C*    AMPSUM(0)  : (C*8) : amplitude sum.
C*          (1)  : (C*8) : sum of non-Higgs diagrams.
C*          (2)  : (C*8) : sum of Higgs diagram.
C* (Relation)
C*    Invokes subroutines in HELAS.LOAD.
C* (Update Record)
C*    09/06/18  K.Fujii         Original version.
C*    09/06/22  J.Tian          another method used
C*
CC**********************************************************************
 
      SUBROUTINE AMPZZZ(GAL,GZN,GZL,GWF,GZH,GWH,GW,GZZZ,
     .                  AME,AMW,GMW,AMZ,GMZ,AMH,GMH,
     .                  VCT,EIN,EOT, AMPSUM)
 
      IMPLICIT   REAL*4  ( A-H, O-Z )
      REAL   *4  GAL(2), GZN(2), GZL(2), GWF(2), GZH, GWH,
     .           GW, GZZZ,
     .           AME, AMW, GMW, AMZ, GMZ, AMH, GMH
      COMPLEX*8  EIN(6), EOT(6), VCT(6,3), AMPSUM(0:2)
C--
      COMPLEX*8  WRK(6,20), AMP(9), TMP(0:2), ZVCT(6,2,3)
C
C========< Entry Point >================================================
C
C--
C  Zero electron width.
C--
      GME = 0

C--
C  Calculate diagrm (1).
C--
      CALL FVIXXX(EIN,VCT(1,1),GZL,AME,GME,WRK(1,1))
      CALL FVOXXX(EOT,VCT(1,2),GZL,AME,GME,WRK(1,2))
      CALL IOVXXX(WRK(1,1),WRK(1,2),VCT(1,3),GZL,AMP(1))

C--
C  Calculate diagrm (2).
C--
      CALL FVIXXX(EIN,VCT(1,1),GZL,AME,GME,WRK(1,3))
      CALL FVOXXX(EOT,VCT(1,3),GZL,AME,GME,WRK(1,4))
      CALL IOVXXX(WRK(1,3),WRK(1,4),VCT(1,2),GZL,AMP(2))

C--
C  Calculate diagrm (3).
C--
      CALL FVIXXX(EIN,VCT(1,2),GZL,AME,GME,WRK(1,5))
      CALL FVOXXX(EOT,VCT(1,1),GZL,AME,GME,WRK(1,6))
      CALL IOVXXX(WRK(1,5),WRK(1,6),VCT(1,3),GZL,AMP(3))

C--
C  Calculate diagrm (4).
C--
      CALL FVIXXX(EIN,VCT(1,2),GZL,AME,GME,WRK(1,7))
      CALL FVOXXX(EOT,VCT(1,3),GZL,AME,GME,WRK(1,8))
      CALL IOVXXX(WRK(1,7),WRK(1,8),VCT(1,1),GZL,AMP(4))

C--
C  Calculate diagrm (5).
C--
      CALL FVIXXX(EIN,VCT(1,3),GZL,AME,GME,WRK(1,9))
      CALL FVOXXX(EOT,VCT(1,1),GZL,AME,GME,WRK(1,10))
      CALL IOVXXX(WRK(1,9),WRK(1,10),VCT(1,2),GZL,AMP(5))

C--
C  Calculate diagrm (6).
C--
      CALL FVIXXX(EIN,VCT(1,3),GZL,AME,GME,WRK(1,11))
      CALL FVOXXX(EOT,VCT(1,2),GZL,AME,GME,WRK(1,12))
      CALL IOVXXX(WRK(1,11),WRK(1,12),VCT(1,1),GZL,AMP(6))

C--
C  Calculate diagrm (7).
C--
      CALL HVVXXX(VCT(1,1),VCT(1,2),GZH,AMH,GMH,WRK(1,13))
      CALL AMPZH (GZH,GZL,AMZ,GMZ,VCT(1,3),WRK(1,13),EIN,EOT, AMP(7))

C--
C  Calculate diagrm (8).
C--
      CALL HVVXXX(VCT(1,1),VCT(1,3),GZH,AMH,GMH,WRK(1,14))
      CALL AMPZH (GZH,GZL,AMZ,GMZ,VCT(1,2),WRK(1,14),EIN,EOT, AMP(8))

C--
C  Calculate diagrm (9).
C--
      CALL HVVXXX(VCT(1,2),VCT(1,3),GZH,AMH,GMH,WRK(1,15))
      CALL AMPZH (GZH,GZL,AMZ,GMZ,VCT(1,1),WRK(1,15),EIN,EOT, AMP(9))

C--
C  Sum non-Higgs and Higgs diagrams.
C--
      AMPSUM(1) = 0
      DO 10 I = 1, 6
         AMPSUM(1) = AMPSUM(1) + AMP(I)
10    CONTINUE
      AMPSUM(2) = AMP(7)+AMP(8)+AMP(9)
      AMPSUM(0) = AMPSUM(1) + AMPSUM(2)

C--
C  That's it.
C--
      RETURN
      END
