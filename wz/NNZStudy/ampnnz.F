CC**********************************************************************
C*
C*===================================================================
C* Subroutine AMPNNZ(GZL,GZN,GWF,GWWZ,AMW,GMW,AME,
C*                                  SPEM,SPEP,SPNU,SPNB,VCTZ,AMP)
C*=============================================================---===
C*
C* (Purpose)
C*    Calculate t-channel amplitude for e+ + e- ---> nu + nuebar + Z.
C*    Inputs are spinors and wave functions thus replaceable with
C*    spinors or wave function made from final stable particles.
C* (Inputs)
C*       GZL(*)   : (R*8) : lepton Z coupling.
C*       GZN(*)   : (R*8) : neutrino Z coupling.
C*       GWF      : (R*8) : weak coupling.
C*       GWWZ     : (R*8) : WWZ coupling.
C*       AMW      : (R*8) : W mass.
C*       GMW      : (R*8) : W width.
C*       AME      : (R*8) : electron mass.
C*       SPEM(*)  : (C*16) : beam electron spinor.
C*       SPEP(*)  : (C*16) : beam positron spinor.
C*       SPNU (*) : (C*16) : final-state neutrino spinor.
C*       SPNB (*) : (C*16) : final-state neutrino-bar spinor.
C*       VCTZ (*) : (C*16) : final-state Z wave function.
C* (Output)
C*       AMP(*)   : (C*16) : amplitudes for t-channel diagram.
C* (Relation)
C*    Invokes subroutines in HELAS.LOAD.
C* (Update Record)
C*    95/04/26  K.Fujii		Original version.
C*
CC**********************************************************************
 
      SUBROUTINE AMPNNZ(GZL,GZN,GWF,GWWZ,AMW,GMW,AME,
     .                                 SPEM,SPEP,SPNU,SPNB,VCTZ,AMP)
 
      IMPLICIT   REAL*8  ( A-H, O-Z )
      REAL   *8  GZL(2), GZN(2), GWF(2), GWWZ, AMW, GMW, AME
      COMPLEX*16 SPEM(6), SPEP(6), SPNU(6), SPNB(6), VCTZ(6), 
     .		 AMP(0:5)
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
         AMN   = 0
         GMN   = 0
         GME   = 0
      ENDIF
C--
C  Amplitude (1).
C--
      CALL JIOXXX(SPEM,SPNU,GWF,AMW,GMW,WRK(1,1))                 
      CALL FVIXXX(SPNB,WRK(1,1),GWF,AME,GME,WRK(1,2))       
      CALL IOVXXX(WRK(1,2),SPEP,VCTZ,GZL, AMP(1))                     
C--
C  Amplitude (2).
C--
      CALL JIOXXX(SPNB,SPEP,GWF,AMW,GMW,WRK(1,3))                 
      CALL FVIXXX(SPEM,VCTZ,GZL,AME,GME,WRK(1,4))       
      CALL IOVXXX(WRK(1,4),SPNU,WRK(1,3),GWF, AMP(2))    
C--
C  Amplitude (3).
C--
      CALL FVOXXX(SPNU,VCTZ,GZN,AMN,GMN,WRK(1,5))       
      CALL IOVXXX(SPEM,WRK(1,5),WRK(1,3),GWF, AMP(3))                     
C--
C  Amplitude (4).
C--
      CALL FVIXXX(SPNB,VCTZ,GZN,AMN,GMN,WRK(1,6))
      CALL IOVXXX(WRK(1,6),SPEP,WRK(1,1),GWF, AMP(4))                     
C--
C  Amplitude (5).
      CALL VVVXXX(WRK(1,3),WRK(1,1),VCTZ,GWWZ, AMP(5))                       
C--
C  Sum up amplitudes.
C--
      AMP(0) = AMP(1) + AMP(2) + AMP(3) + AMP(4) + AMP(5)
C--
C  That's it.
C--
      RETURN
      END
