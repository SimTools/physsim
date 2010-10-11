C***********************************************************************
C*
C*==========================================================
C* Subrouine SGXX00(GAL,GAX,GZL,GZX,GSNX,AMSN,GMSN,AMZ,GMZ,
C*                                AMX, GMX ,RS, POL, SGEFF)
C*===================================================-----==
C*
C* (Purpose)
C*   Calculates the tree level total cross section for
C*   chargino pair production.
C* (Inputs)
C*      GAL(2)    : (R*4) : L-A-L chiral couplings.
C*      GAX(2)    : (R*4) : XC-A-XC chiral couplings.
C*      GZL(2)    : (R*4) : L-Z-L chiral couplings.
C*      GZX(2)    : (R*4) : XC-Z-XC chiral couplings.
C*      GSNX(2,2) : (R*4) : E-SNUE-XC couplings.
C*      AMSN      : (R*4) : sneutrino mass.
C*      GMSN      : (R*4) : sneutrino width (not used).
C*      AMZ       : (R*4) : Z mass.
C*      GMZ       : (R*4) : Z width.
C*      AMX(2)    : (R*4) : XC masses.
C*      GMX(2)    : (R*4) : XC widths (not used).
C*      RS        : (R*4) : sqrt(s).
C*      POL       : (R*4) : electron beam polarization.
C* (Output)
C*      SGEFF     : (R*4) : effective cross section.
C* (Relation)
C*   Invokes 
C*		SGCXXA
C* (Update Record)
C*   92/06/26  K.Fujii       Original version. R*4 interface to
C*			     SGCXXA.
C*
C***********************************************************************
 
      SUBROUTINE SGXX00(GAL,GAX,GZL,GZX,GSNX,AMSN,GMSN,AMZ,GMZ,
     .                  AMX, GMX ,RS, POL, SGEFF)
 
      IMPLICIT     REAL*8 ( A-H, O-Z )
      REAL   *4    GAL(2), GAX(2), GZL(2), GZX(2), GSNX(2,2), POL,
     .             AMSN, GMSN, AMZ, GMZ, AMX(2), GMX(2), RS, SGEFF
C--
      REAL   *8    GALD(2), GAXD(2), GZLD(2), GZXD(2), GSNXD(2,2),
     .             AMXD(2), GMXD(2), SG(0:3)
C
C========< Entry Point >================================================
C
C--
C  Convert the inputs to real*8 variables.
C--
         GALD(1)    = GAL(1)
         GALD(2)    = GAL(2)
         GAXD(1)    = GAX(1)
         GAXD(2)    = GAX(2)
         GZLD(1)    = GZL(1)
         GZLD(2)    = GZL(2)
         GZXD(1)    = GZX(1)
         GZXD(2)    = GZX(2)
         GSNXD(1,1) = GSNX(1,1)
         GSNXD(2,1) = GSNX(2,1)
         GSNXD(1,2) = GSNX(1,2)
         GSNXD(2,2) = GSNX(2,2)
         AMSND      = AMSN
         GMSND      = GMSN
         AMZD       = AMZ
         GMZD       = GMZ
         AMXD(1)    = AMX(1)
         AMXD(2)    = AMX(2)
         GMXD(1)    = GMX(1)
         GMXD(2)    = GMX(2)
         RSD        = RS
         POLD       = POL
         CALL SGCXXA(GALD,GAXD,GZLD,GZXD,GSNXD, AMSND,GMSND,
     .               AMZD,GMZD,AMXD,GMXD, RSD, POLD, SG)
         SGEFF = SG(0)
C--
C  That's it.
C--
      RETURN
      END
