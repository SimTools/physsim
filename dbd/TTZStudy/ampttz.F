CC**********************************************************************
C*
C*=================================================================
C* Subroutine AMPTTZ(GAL,GAU,GZL,GZU,GZH,GTH,
C*                   AMZ,GMZ,AME,AMT,GMT,AMH,GMH, ZVCT,
C*                   EIN,EOT,TIN,TOT, AMPSUM)
C*====================================------=======================
C*
C* (Purpose)
C*    Calculate amplitude for e+ + e- ---> t + t-bar + Z.
C*    Inputs are spinors and wave function thus replaceable with
C*    spinors or wave function made from final stable particle
C*    4-momenta.
C* (Inputs)
C*       GAL(*)  : (R*4) : lepton photon coupling.
C*       GAU(*)  : (R*4) : up-type quark photon coupling.
C*       GZL(*)  : (R*4) : lepton Z coupling.
C*       GZU(*)  : (R*4) : up-type quark Z coupling.
C*       GZH     : (R*4) : ZZ-H coupling.
C*       GTH(*)  : (C*8) : TT-H COUPLING.
C*       AMZ     : (R*4) : Z mass.
C*       GMZ     : (R*4) : Z width.
C*       AME     : (R*4) : electron mass.
C*       AMT     : (R*4) : top mass.
C*       GMT     : (R*4) : top width.
C*       AMH     : (R*4) : Higss mass.
C*       GMH     : (R*4) : Higss width.
C*       ZVCT(*) : (C*8) : Z wave function.
C*       EIN(*)  : (C*8) : incoming electron spinor.
C*       EOT(*)  : (C*8) : outgoing electron spinor.
C*       TIN(*)  : (C*8) : incoming top spinor.
C*       TOT(*)  : (C*8) : outgoing top spinor.
C* (Output)
C*    AMPSUM(0)  : (C*8) : amplitude sum.
C*          (1)  : (C*8) : sum of non-Higgs diagrams.
C*          (2)  : (C*8) : sum of Higgs diagram.
C* (Relation)
C*    Invokes AMPTT, AMPZH, and subroutines in HELAS.LOAD.
C* (Update Record)
C*    90/06/13  K.Fujii         Original version.
C*    90/08/02  K.Fujii         Modified to accomodate with HELAS 1.06.
C*    95/03/16  K.Fujii         Modified to accomodate with HELAS 2.04.
C*    95/03/28  K.Fujii         Modified the calling sequence of AMPZH 
C*				to be consistent with the one for 
C*				ffh_gen.
C*
CC**********************************************************************
 
      SUBROUTINE AMPTTZ(GAL,GAU,GZL,GZU,GZH,GTH,
     .                  AMZ,GMZ,AME,AMT,GMT,AMH,GMH, ZVCT,
     .                  EIN,EOT,TIN,TOT, AMPSUM)
 
      IMPLICIT   REAL*4  ( A-H, O-Z )
      REAL   *4  GAL(2), GAU(2), GZL(2), GZU(2), GZH,
     .           AMZ, GMZ, AME, AMT, GMT, AMH, GMH
      COMPLEX*8  GTH(2)
      COMPLEX*8  EIN(6), EOT(6), TIN(6), TOT(6), ZVCT(0:5), AMPSUM(0:2)
C--
      COMPLEX*8  VEIN(6), VEOT(6), VTIN(6), VTOT(6), SCLR(3), AMP(0:2)
C
C========< Entry Point >================================================
C
C--
C  Prepare spinors for virtual e, e-bar, t, and t-bar.
C--
      GME = 0
      CALL FVIXXX(EIN,ZVCT,GZL,AME,GME,VEIN)
      CALL FVIXXX(TIN,ZVCT,GZU,AMT,GMT,VTIN)
      CALL FVOXXX(EOT,ZVCT,GZL,AME,GME,VEOT)
      CALL FVOXXX(TOT,ZVCT,GZU,AMT,GMT,VTOT)
C--
C  Calculate non-Higgs diagrams.
C--
      AMPSUM(1) = 0
C--
      CALL AMPTT (GAL,GAU,GZL,GZU,AMZ,GMZ,
     .            VEIN, EOT,  TIN,  TOT, AMP)
      AMPSUM(1) = AMPSUM(1) + AMP(0)
C--
      CALL AMPTT (GAL,GAU,GZL,GZU,AMZ,GMZ,
     .            EIN, VEOT,  TIN,  TOT, AMP)
      AMPSUM(1) = AMPSUM(1) + AMP(0)
C--
      CALL AMPTT (GAL,GAU,GZL,GZU,AMZ,GMZ,
     .            EIN,  EOT, VTIN,  TOT, AMP)
      AMPSUM(1) = AMPSUM(1) + AMP(0)
C--
      CALL AMPTT (GAL,GAU,GZL,GZU,AMZ,GMZ,
     .            EIN,  EOT,  TIN, VTOT, AMP)
      AMPSUM(1) = AMPSUM(1) + AMP(0)
C--
C  Calculate Higgs diagram.
C--
      CALL HIOXXX(TIN,TOT,GTH,AMH,GMH,SCLR)
      CALL AMPZH (GZL,GZH, AMZ,GMZ, ZVCT, SCLR, EIN,EOT, AMP)
      AMPSUM(2) = AMP(0)
C--
C  Sum non-Higgs and Higgs diagrams.
C--
      AMPSUM(0) = AMPSUM(1) + AMPSUM(2)
C--
C  That's it.
C--
      RETURN
      END
