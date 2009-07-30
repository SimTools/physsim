CC**********************************************************************
C*
C*==============================================================
C* Subroutine AMPTTG(GAL,GAF,GZL,GZF,GCH, AMZ,GMZ,AMF,GMF,
C*                   EIN,EOT, FIN,FOT,HOT, AMP)
C*=========================================---==================
C*
C* (Purpose)
C*    Calculate amplitude for e+ + e- ---> f + f-bar gluon.
C* (Inputs)
C*       GAL(2)   : (R*4) : eeA coupling.
C*       GAF(2)   : (R*4) : ffA coupling.
C*       GZL(2)   : (R*4) : eeZ coupling.
C*       GZF(2)   : (R*4) : ffZ coupling.
C*       GCH(2)   : (C*8) : ffG coupling.
C*       GCA(2)   : (C*8) : ffA coupling.
C*       AMZ      : (R*4) : Z mass.
C*       GMZ      : (R*4) : Z width.
C*       AMF      : (R*4) : f mass.
C*       GMF      : (R*4) : f width.
C*       EIN(6)   : (C*8) : e- spinor.
C*       EOT(6)   : (C*8) : e+ spinor.
C*       FIN(6)   : (C*8) : f-bar spinor.
C*       FOT(6)   : (C*8) : f spinor.
C*       HOT(6)   : (C*8) : vector.
C* (Outputs)
C*    AMP         : (C*8) : amplitude.
C* (Relation)
C*    Invokes AMPFF and subroutines in HELAS.LOAD.
C* (Update Record)
C*   2009/07/30  T.Tanabe	Inital version based on TTHStudy by K.Fujii
C*
CC**********************************************************************
 
      SUBROUTINE AMPTTG(GAL,GAF,GZL,GZF,GCH, AMZ,GMZ,AMF,GMF,
     .                  EIN,EOT, FIN,FOT,HOT, AMP)
 
      IMPLICIT   REAL*4  ( A-H, O-Z )
      REAL   *4  GAL(2), GAF(2), GZL(2), GZF(2),
     .           AMZ, GMZ, AMF, GMF
      REAL   *4  GCH(2)
      COMPLEX*8  EIN(6), EOT(6), FIN(6), FOT(6), HOT(6), 
     .           AMP
C--
      COMPLEX*8  SPWRK(6), SCWRK(3), VCWRK(6), TMP(0:2)
C
C========< Entry Point >================================================
C
C--
C  Bremsstrahlung off f-bar.
C--
      CALL FVIXXX(FIN,HOT,GCH,AMF,GMF,SPWRK)
      CALL AMPTT (GAL,GAF,GZL,GZF,AMZ,GMZ,
     .            EIN,EOT,SPWRK,FOT,TMP)
      AMP = TMP(0)
C--
C  Bremsstrahlung off f.
C--
      CALL FVOXXX(FOT,HOT,GCH,AMF,GMF,SPWRK)
      CALL AMPTT (GAL,GAF,GZL,GZF,AMZ,GMZ,
     .            EIN,EOT,FIN,SPWRK,TMP)
      AMP = AMP + TMP(0)
C--
C  That's it.
C--
      RETURN
      END
