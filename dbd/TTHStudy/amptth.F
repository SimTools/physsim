CC**********************************************************************
C*
C*==============================================================
C* Subroutine AMPTTH(GAL,GAF,GZL,GZF,GVH,GCH, AMZ,GMZ,AMF,GMF,
C*                   EIN,EOT, FIN,FOT,HOT, AMP)
C*=========================================---==================
C*
C* (Purpose)
C*    Calculate amplitude for e+ + e- ---> f + f-bar H.
C* (Inputs)
C*       GAL(2)   : (R*4) : eeA coupling.
C*       GAF(2)   : (R*4) : ffA coupling.
C*       GZL(2)   : (R*4) : eeZ coupling.
C*       GZF(2)   : (R*4) : ffZ coupling.
C*       GVH      : (R*4) : VVH coupling.
C*       GCH(2)   : (C*8) : ffH coupling.
C*       GCA(2)   : (C*8) : ffA coupling.
C*       AMZ      : (R*4) : Z mass.
C*       GMZ      : (R*4) : Z width.
C*       AMF      : (R*4) : f mass.
C*       GMF      : (R*4) : f width.
C*       EIN(6)   : (C*8) : e- spinor.
C*       EOT(6)   : (C*8) : e+ spinor.
C*       FIN(6)   : (C*8) : f-bar spinor.
C*       FOT(6)   : (C*8) : f spinor.
C*       HOT(3)   : (C*8) : scalar.
C* (Outputs)
C*    AMP (0)     : (C*8) : amplitude sum.
C*        (1)     : (C*8) : bremsstrahlung.
C*        (2)     : (C*8) : ZH.
C* (Relation)
C*    Invokes AMPFF and subroutines in HELAS.LOAD.
C* (Update Record)
C*   95/03/28  K.Fujii	Derived from AMPFFH in ffh_gen.
C*
CC**********************************************************************
 
      SUBROUTINE AMPTTH(GAL,GAF,GZL,GZF,GVH,GCH, AMZ,GMZ,AMF,GMF,
     .                  EIN,EOT, FIN,FOT,HOT, AMP)
 
      IMPLICIT   REAL*4  ( A-H, O-Z )
      REAL   *4  GAL(2), GAF(2), GZL(2), GZF(2), GVH,
     .           AMZ, GMZ, AMF, GMF
      COMPLEX*8  GCH(2), 
     .           EIN(6), EOT(6), FIN(6), FOT(6), HOT(3), 
     .           AMP(0:2)
C--
      COMPLEX*8  SPWRK(6), SCWRK(3), VCWRK(6), TMP(0:2)
C
C========< Entry Point >================================================
C
C--
C  Bremsstrahlung off f-bar.
C--
      CALL FSIXXX(FIN,HOT,GCH,AMF,GMF,SPWRK)
      CALL AMPTT (GAL,GAF,GZL,GZF,AMZ,GMZ,
     .            EIN,EOT,SPWRK,FOT,TMP)
      AMP(1) = TMP(0)
C--
C  Bremsstrahlung off f.
C--
      CALL FSOXXX(FOT,HOT,GCH,AMF,GMF,SPWRK)
      CALL AMPTT (GAL,GAF,GZL,GZF,AMZ,GMZ,
     .            EIN,EOT,FIN,SPWRK,TMP)
      AMP(1) = AMP(1) + TMP(0)
C--
C  ZH and Z to ff-bar.
C--
      CALL JIOXXX(FIN,FOT,GZF,AMZ,GMZ,VCWRK)
      CALL AMPZH (GZL,GVH,AMZ,GMZ,VCWRK,HOT,EIN,EOT,AMP(2))
C--
C  Sum up amplitudes.
C--
      AMP(0) = AMP(1) + AMP(2)
C--
C  That's it.
C--
      RETURN
      END
