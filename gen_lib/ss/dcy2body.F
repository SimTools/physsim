C==============
C subrouines to calculate decay widthes of 2-body decay channel
C
C      Fermion_1 ---> Fermion_2 + Vector_boson
C      Scaler    ---> Fermion_1 + Fermion_2
C      Fermion_1 ---> Fermion_2 + Scaler
C
C  << Note >>
C  If both oF fermion_1 & 2 are Majorana particles,
C  one should give appropreate GL & GR or keep in mind that
C  output is not correct because these routines
C  don't concern the case.
C
C==============
C----------------------------------------------------------------------C
      SUBROUTINE F1F2VD( F1MAS, F2MAS, VMAS, GL, GR, DCYWDT )
C----------------------------------------------------------------------C
C Purpose
C     calculate the decay widthes of
C       Fermion_1 --> Fermion_2 + Vector_boson
C
C Inputs ( note that GL,GR are complex*8 )
C     F1MAS (Real*4) :   mass of fermion_1 (GeV)
C     F2MAS (Real*4) :   mass of fermion_2 (GeV)
C     VMAS  (Real*4) :   mass of Vector_boson  (GeV)
C     GL, GR (complex*8) :  left/right-handed coupling of
C                           f_1 - f_2 - V
C Output
C     DCYWDT : decay width (GeV)
C----------------------------------------------------------------------C
      COMPLEX*8 GL, GR
C----------------------------------------------------------------------C
      IF( F1MAS .LE. F2MAS+VMAS ) THEN
        DCYWDT = 0.
        RETURN
      ENDIF
 
      PI = ACOS(0.)*2.
      F1MAS2 = F1MAS**2
      F2MAS2 = F2MAS**2
      VMAS2 = VMAS**2
 
      PMOM2 = (F1MAS2-(F2MAS+VMAS)**2)*(F1MAS2-(F2MAS-VMAS)**2)
     &         /4./F1MAS2
      PMOM = SQRT(PMOM2)
 
      TMAT2 = 0.5* ( ABS(GL)**2 + ABS(GR)**2 )
     &        * ( F1MAS2+F2MAS2-VMAS2
     &             + ((F1MAS2-F2MAS2)**2-VMAS2**2)/VMAS2      )
     &       - 6.*F1MAS*F2MAS*REAL( GL*CONJG(GR) )
 
      DCYWDT = PMOM*TMAT2/8./PI/F1MAS2
 
 
      RETURN
      END
C--------------
      SUBROUTINE SC1SC2V(SCMAS1,SCMAS2,VMAS, GV, DCYWDT )
C Purpose
C     calculate the decay widthes of
C       Scaler_1 ---> Scaler_2 + V
C Inputs ( note that GL,GR are complex*8 )
C     SCMAS1  (Real*4) :    mass of scaler1 (GeV)
C     SCMAS2  (Real*4) :    mass of scaler_2  (GeV)
C     VMAS  (Real*4) :    mass of vector  (GeV)
C     GV (complex*8) :   coupling of  scaler1 - scaler2 - V
C                       
C Output
C     DCYWDT : decay width (GeV)
C         
C----------------------------------------------------------------------C
      COMPLEX*8 GV
C----------------------------------------------------------------------C
      IF( SCMAS1 .LE. SCMAS2+VMAs) THEN
        DCYWDT = 0.
        RETURN
      ENDIF
 
      PI = ACOS(0.)*2.
      SCMS12 = SCMAS1**2
      SCMS22 = SCMAS2**2
      VMAS2 = VMAS**2
      BETA2 = ( 1. - (SCMAS2 + VMAS)**2/SCMS12)*
     &         ( 1. - (SCMAS2 - VMAS)**2/SCMS12)
      BETA =SQRT(BETA2)
      TMAT2= (SCMS12 - SCMS22)**2/VMAS2 -1.5*(SCMS12+SCMS22) 
     & + 0.5*VMAS2
      DCYWDT = BETA * TMAT2 * ABS(GV)**2/16./PI/SCMAS1
 
 
      RETURN
      END

C----------------------------------------------------------------------C
      SUBROUTINE SCF1F2( SCMAS, F1MAS, F2MAS, GL, GR, DCYWDT )
C----------------------------------------------------------------------C
C Purpose
C     calculate the decay widthes of
C       Scaler ---> fermion_1 + fermion_2
C Inputs ( note that GL,GR are complex*8 )
C     SCMAS  (Real*4) :    mass of scaler (GeV)
C     F1MAS  (Real*4) :    mass of fermion_1  (GeV)
C     F2MAS  (Real*4) :    mass of fermion_2  (GeV)
C     GL, GR (complex*8) :  left/right-handed coupling of
C                           scaler - f_1 - f_2
C Output
C     DCYWDT : decay width (GeV)
C----------------------------------------------------------------------C
      COMPLEX*8 GL, GR
C----------------------------------------------------------------------C
      IF( SCMAS .LE. F1MAS+F2MAS ) THEN
        DCYWDT = 0.
        RETURN
      ENDIF
 
      PI = ACOS(0.)*2.
      SCMAS2 = SCMAS**2
      F1MAS2 = F1MAS**2
      F2MAS2 = F2MAS**2
 
      PMOM2 = (SCMAS2-(F1MAS+F2MAS)**2)*(SCMAS2-(F1MAS-F2MAS)**2)
     &         /4./SCMAS2
      PMOM = SQRT(PMOM2)
 
      TMAT2 = ( ABS(GL)**2 + ABS(GR)**2 ) * ( SCMAS2-F1MAS2-F2MAS2 )
     &       - 4.*F2MAS*F1MAS*REAL( GL*CONJG(GR) )
 
      DCYWDT = PMOM*TMAT2/8./PI/SCMAS2
 
 
      RETURN
      END
C----------------------------------------------------------------------C
      SUBROUTINE F1F2SC( F1MAS, F2MAS, SCMAS, GL, GR, DCYWDT )
C----------------------------------------------------------------------C
C Purpose
C     calculate the decay widthes of
C       fermion_1 --> fermion_2 + scaler
C Inputs ( note that GL,GR are complex*8 )
C     F1MAS  (Real*4) :    mass of fermion_1  (GeV)
C     F2MAS  (Real*4) :    mass of fermion_2  (GeV)
C     SCMAS  (Real*4) :    mass of scaler (GeV)
C     GL, GR (complex*8) :  left/right-handed coupling of
C                           f_1 - f_2 - scaler
C Output
C     DCYWDT : decay width (GeV)
C----------------------------------------------------------------------C
      COMPLEX*8 GL, GR
C----------------------------------------------------------------------C
      IF( F1MAS .LE. F2MAS+SCMAS ) THEN
        DCYWDT = 0.
        RETURN
      ENDIF
 
      PI = ACOS(0.)*2.
      F1MAS2 = F1MAS**2
      F2MAS2 = F2MAS**2
      SCMAS2 = SCMAS**2
 
      PMOM2 = (F1MAS2-(F2MAS+SCMAS)**2)*(F1MAS2-(F2MAS-SCMAS)**2)
     &         /4./F1MAS2
      PMOM = SQRT(PMOM2)
 
      TMAT2 = ( ABS(GL)**2 + ABS(GR)**2 ) * ( F1MAS2+F2MAS2-SCMAS2 )
     &       + 4.*F1MAS*F2MAS*REAL( GL*CONJG(GR) )
C      spin average for initial fermion f1
      TMAT2 = TMAT2/2.
 
      DCYWDT = PMOM*TMAT2/8./PI/F1MAS2
 
 
      RETURN
      END
