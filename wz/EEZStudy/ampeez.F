CC**********************************************************************
C*
C*====================================================
C* Subroutine AMPEEZ(GAL,GZL, AMZ,GMZ, AME,      EBM,
C*  		     E1, SHLF1,CHLF1,FI1, NHEB1,NHE1,
C* 		     E2, SHLF2,CHLF2,FI2, NHEB2,NHE2,
C*                       SPEM,SPEP, SP1,SP2,VCT, AMP)
C*===============================================---==
C*
C* (Purpose)
C*    Calculate t-channel amplitude for e+ e- ---> e+ e- W+ W-.
C*    Inputs are spinors and wave functions thus replaceable with
C*    spinors or wave functions made from final stable particles.
C* (Inputs)
C*       GAL(*)   : (R*8) : lepton photon coupling.
C*       GZL(*)   : (R*8) : lepton Z coupling.
C*       AMZ      : (R*8) : Z mass.
C*       GMZ      : (R*8) : Z width.
C*       AME      : (R*8) : electron mass.
C*       EBM      : (R*8) : beam energy.
C*       E1       : (R*8) : final-state e- energy.
C*       SHLF1    : (R*8) : sin(theta_e-/2).
C*       CHLF1    : (R*8) : cos(theta_e-/2).
C*       FI1      : (R*8) : phi_e-
C*       NHEB1    : (I*4) : beam e- helicity.
C*       NHE1     : (I*4) : final-state e- helicity.
C*       E2       : (R*8) : final-state e+ energy.
C*       SHLF2    : (R*8) : sin(theta_e+/2).
C*       CHLF2    : (R*8) : cos(theta_e+/2).
C*       FI2      : (R*8) : phi_e+
C*       NHEB2    : (I*4) : beam e+ helicity.
C*       NHE2     : (I*4) : final-state e+ helicity.
C*       SPEM(*)  : (C*16) : beam e- spinor.
C*       SPEP(*)  : (C*16) : beam e+ spinor.
C*       SP1 (*)  : (C*16) : final-state e- spinor.
C*       SP2 (*)  : (C*16) : final-state e+ spinor.
C*       VCT (*)  : (C*16) : final-state Z wave function.
C* (Output)
C*       AMP      : (C*16) : total t-channel amplitudes.
C* (Relation)
C*    Invokes subroutines in HELAS.LOAD.
C* (Update Record)
C*    95/04/26  K.Fujii		Original version.
C*
CC**********************************************************************
 

      SUBROUTINE AMPEEZ(GAL,GZL, AMZ,GMZ, AME,      EBM,
     . 		        E1, SHLF1,CHLF1,FI1, NHEB1,NHE1,
     .		        E2, SHLF2,CHLF2,FI2, NHEB2,NHE2,
     .                      SPEM,SPEP, SP1,SP2,VCT, AMP)
 
      IMPLICIT   REAL*8  ( A-H, O-Z )
      REAL   *8  GAL(2), GZL(2), AMZ, GMZ, AME, EBM,
     .		 E1, SHLF1, CHLF1, FI1,  E2, SHLF2, CHLF2, FI2    
      INTEGER*4  NHEB1, NHE1, NHEB2, NHE2
      COMPLEX*16 SPEM(6), SPEP(6), SP1(6), SP2(6), VCT(6), AMP 
C--
      COMPLEX*16 WRK(6,4), TMP(6), ZMP(8)
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
         GME   = 0
      ENDIF
C--
C  Calculate internal lines.
C--
      CALL JEEXXX(EBM,E1,SHLF1,CHLF1,FI1,NHEB1,NHE1,+1, WRK(1, 1))
      CALL JIOXXX(SPEM,SP1,GZL,AMZ,GMZ,                 WRK(1, 2))
         WRK(1,1) = FACT*WRK(1,1)
         WRK(2,1) = FACT*WRK(2,1)
         WRK(3,1) = FACT*WRK(3,1)
         WRK(4,1) = FACT*WRK(4,1)
C--
      CALL JEEXXX(EBM,E2,SHLF2,CHLF2,FI2,NHEB2,NHE2,-1, WRK(1, 3))
      CALL JIOXXX(SP2,SPEP,GZL,AMZ,GMZ,                 WRK(1, 4))
         WRK(1,3) = FACT*WRK(1,3)
         WRK(2,3) = FACT*WRK(2,3)
         WRK(3,3) = FACT*WRK(3,3)
         WRK(4,3) = FACT*WRK(4,3)
C--
C  Now calculate amplitudes.
C--
      CALL FVIXXX(SP2, VCT,GZL,AME,GME,TMP)
      CALL IOVXXX(TMP,SPEP,WRK(1, 1),GAL,ZMP(1))
      CALL IOVXXX(TMP,SPEP,WRK(1, 2),GZL,ZMP(2))
C--
      CALL FVOXXX(SPEP,VCT,GZL,AME,GME,TMP)      
      CALL IOVXXX(SP2, TMP,WRK(1, 1),GAL,ZMP(3))
      CALL IOVXXX(SP2, TMP,WRK(1, 2),GZL,ZMP(4))
C--
      CALL FVOXXX(SP1, VCT,GZL,AME,GME,TMP)
      CALL IOVXXX(SPEM,TMP,WRK(1, 3),GAL,ZMP(5))
      CALL IOVXXX(SPEM,TMP,WRK(1, 4),GZL,ZMP(6))
C--
      CALL FVIXXX(SPEM,VCT,GZL,AME,GME,TMP)
      CALL IOVXXX(TMP, SP1,WRK(1, 3),GAL,ZMP(7))
      CALL IOVXXX(TMP, SP1,WRK(1, 4),GZL,ZMP(8))
C--
C  Sum up amplitudes.
C--
      AMP = ZMP(1) + ZMP(2) + ZMP(3) + ZMP(4) 
     .    + ZMP(5) + ZMP(6) + ZMP(7) + ZMP(8) 
C--
C  That's it.
C--
      RETURN
      END
