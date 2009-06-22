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
C  Calculate virtual e-.
C--
      CALL FVIXXX(EIN,VCT(1,1),GZL,AME,GME,WRK(1,1))
      CALL FVIXXX(EIN,VCT(1,2),GZL,AME,GME,WRK(1,2))
      CALL FVIXXX(EIN,VCT(1,3),GZL,AME,GME,WRK(1,3))
C--
C  Calculate virtual e+.
C--
      CALL FVOXXX(EOT,VCT(1,1),GZL,AME,GME,WRK(1,4))
      CALL FVOXXX(EOT,VCT(1,2),GZL,AME,GME,WRK(1,5))
      CALL FVOXXX(EOT,VCT(1,3),GZL,AME,GME,WRK(1,6))
C--
C  Set ZVCT.
C--
      DO 100 IZ = 1, 3
         IF (IZ .EQ. 1) THEN
            IZ1 = 2
            IZ2 = 3
         ELSE IF (IZ .EQ. 2) THEN
            IZ1 = 3
            IZ2 = 1
         ELSE
            IZ1 = 1
            IZ2 = 2
         ENDIF
         DO 10 I = 1, 6
            ZVCT(I,1,IZ) = VCT(I,IZ1)
            ZVCT(I,2,IZ) = VCT(I,IZ2)
10       CONTINUE
100   CONTINUE
C--
C  Calculate diagrams.
C--
      DO 200 IZ = 1, 3
C--      Non-H
         CALL AMPZZ (GZL, AME,GME, ZVCT(1,1,IZ), WRK(1,IZ),EOT, TMP)
         AMP(IZ  ) = TMP(0)
         CALL AMPZZ (GZL, AME,GME, ZVCT(1,1,IZ), EIN,WRK(1,IZ), TMP)
         AMP(IZ+3) = TMP(0)
C--      Higgs
         CALL HVVXXX(ZVCT(1,1,IZ),ZVCT(1,2,IZ),GZH,AMH,GMH,WRK(1,7))
         CALL AMPZH (GZH,GZL,AMZ,GMZ,VCT(1,IZ),WRK(1,7),EIN,EOT, TMP)
         AMP(IZ+6) = TMP(0)
200   CONTINUE
C--
C  Sum non-Higgs and Higgs diagrams.
C--
      AMPSUM(1) = 0
      DO 30 I = 1, 3
         AMPSUM(1) = AMPSUM(1) + AMP(I) + AMP(I+3)
         AMPSUM(2) = AMPSUM(2) + AMP(I+6)
30    CONTINUE
      AMPSUM(0) = AMPSUM(1) + AMPSUM(2)
C--
C  That's it.
C--
      RETURN
      END
