CC**********************************************************************
C*
C*===================================================================
C* Subroutine AMPENW(GAL,GZL,GZN,GWF,GWWA,GWWZ,AMZ,GMZ,AMW,GMW,AME,
C*    				    EBM,EE,SHLF,CHLF,FIE,
C*   				    NHEB,NHEF,
C*                                  SPEMB,SPEPB,SPEM,SPNB,VCWP,AMP)
C*=============================================================---===
C*
C* (Purpose)
C*    Calculate t-channel amplitude for e+ + e- ---> e- + nuebar + W+.
C*    Inputs are spinors and wave functions thus replaceable with
C*    spinors or wave function made from final stable particles.
C* (Inputs)
C*       GAL(*)   : (R*8) : lepton photon coupling.
C*       GZL(*)   : (R*8) : lepton Z coupling.
C*       GZN(*)   : (R*8) : neutrino Z coupling.
C*       GWF      : (R*8) : weak coupling.
C*       GWWA     : (R*8) : WWA coupling.
C*       GWWZ     : (R*8) : WWZ coupling.
C*       AMZ      : (R*8) : Z mass.
C*       GMZ      : (R*8) : Z width.
C*       AMW      : (R*8) : W mass.
C*       GMW      : (R*8) : W width.
C*       AME      : (R*8) : electron mass.
C*       EBM      : (R*8) : beam energy.
C*       EE       : (R*8) : final-state electron energy.
C*       SHLF     : (R*8) : sin(theta_e/2).
C*       CHLF     : (R*8) : cos(theta_e/2).
C*       FIEF     : (R*8) : phi_e
C*       NHEB     : (I*4) : beam electron helicity.
C*       NHEF     : (I*4) : final-state electron helicity.
C*       SPEMB(*) : (C*16) : beam electron spinor.
C*       SPEPB(*) : (C*16) : beam positron spinor.
C*       SPEM (*) : (C*16) : final-state electron spinor.
C*       SPNB (*) : (C*16) : final-state neutrino-bar spinor.
C*       VCWP (*) : (C*16) : final-state W+ wave function.
C* (Output)
C*       AMP     : (C*16) : amplitude for t-channel diagram.
C* (Relation)
C*    Invokes subroutines in HELAS.LOAD.
C* (Update Record)
C*    95/04/20  K.Fujii		Original version.
C*
CC**********************************************************************
 
      SUBROUTINE AMPENW(GAL,GZL,GZN,GWF,GWWA,GWWZ,AMZ,GMZ,AMW,GMW,AME,
     . 				       EBM,EE,SHLF,CHLF,FIE,
     .				       NHEB,NHEF,
     .                                 SPEMB,SPEPB,SPEM,SPNB,VCWP,AMP)
 
      IMPLICIT   REAL*8  ( A-H, O-Z )
      REAL   *8  GAL(2), GZL(2), GZN(2), GWF(2), GWWA, GWWZ,
     .           AMZ, GMZ, AMW, GMW, AME,
     .		 EBM, EE, SHLF, CHLF, FIE
      INTEGER*4  NHEB, NHEF
      COMPLEX*16 SPEMB(6), SPEPB(6), SPEM(6), SPNB(6), VCWP(6), 
     .		 AMP(0:6)
C--
      COMPLEX*16 WRK(6,6)
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
         PI    = ACOS(-1.D0)
         FACT  = ABS(GAL(1))/SQRT(4*PI/128.D0)
         AMN   = 0
         GMN   = 0
         GME   = 0
      ENDIF
C--
C  Calculate t-channel currents.
C--
      CALL JEEXXX(EBM,EE,SHLF,CHLF,FIE,NHEB,NHEF,+1, WRK(1,1))
      CALL JIOXXX(SPEMB,SPEM,GZL,AMZ,GMZ, WRK(1,2))
      CALL JIOXXX(SPNB,SPEPB,GWF,AMW,GMW, WRK(1,3))
C--
C  Calculate virtual fermion spinors.
C--
      CALL FVOXXX(SPEPB,VCWP,GWF,AMN,GMN, WRK(1,4))
      CALL FVIXXX(SPNB ,VCWP,GWF,AME,GME, WRK(1,5))
      CALL FVOXXX(SPEM ,VCWP,GWF,AMN,GMN, WRK(1,6))
C--
C  Diagram (a).
C--
      CALL IOVXXX(WRK(1,5),SPEPB,WRK(1,1),GAL, AMP(1))
      CALL IOVXXX(WRK(1,5),SPEPB,WRK(1,2),GZL, AMP(2))
C--
C  Diagram (b).
C--
      CALL VVVXXX(WRK(1,3),VCWP,WRK(1,1),GWWA, AMP(3))
      CALL VVVXXX(WRK(1,3),VCWP,WRK(1,2),GWWZ, AMP(4))
C--
C  Diagram (c).
C--
      CALL IOVXXX(SPNB,WRK(1,4),WRK(1,2),GZN,  AMP(5))
C--
C  Diagram (d).
C--
      CALL IOVXXX(SPEMB,WRK(1,6),WRK(1,3),GWF, AMP(6))
C--
C  Corrections for alpha modification.
C--
      AMP(1) = FACT*AMP(1)
      AMP(3) = FACT*AMP(3)
C--
C  Sum up amplitudes.
C--
      AMP(0) = AMP(1) + AMP(2) + AMP(3) + AMP(4) + AMP(5) + AMP(6)
C--
C  That's it.
C--
      RETURN
      END
