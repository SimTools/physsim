CC**********************************************************************
C*
C*===================================---===
C* Subroutine AMPWWT(GWF,VCT,EIN,EOT,AMP)
C*===================================---===
C*
C* (Purpose)
C*    Calculate t-channel amplitude for e+ + e- ---> W+ + W-.
C*    Inputs are spinors and wave functions thus replaceable with
C*    spinors or wave function made from final stable particle
C*    4-momenta.
C* (Inputs)
C*       GWF(*)  : (R*4) : fermion W coupling.
C*       VCT(*,i): (C*8) : i = (1,2) = (W-,W+).
C*       EIN(*)  : (C*8) : incoming electron spinor.
C*       EOT(*)  : (C*8) : outgoing electron spinor.
C* (Output)
C*       AMP     : (C*8) : amplitude for t-channel diagram.
C* (Relation)
C*    Invokes subroutines in HELAS.LOAD.
C* (Update Record)
C*    90/08/08  K.Fujii         Original version.
C*    95/03/12  K.Fujii		Modiifed to HELAS V204.
C*
CC**********************************************************************
 
      SUBROUTINE AMPWWT(GWF,VCT,EIN,EOT,AMP)
 
      IMPLICIT   REAL*4  ( A-H, O-Z )
      REAL   *4  GWF(2)
      COMPLEX*8  EIN(6), EOT(6), VCT(6,2), AMP
C--
      COMPLEX*8  SPINOR(6)
C
C========< Entry Point >================================================
C
C--
C  Zero neutrino mass and width.
C--
      AMN = 0
      GMN = 0
C--
C  Calculate the t-channel diagrams.
C--
      CALL FVIXXX(EIN,VCT(1,1),GWF,AMN,GMN,SPINOR)
      CALL IOVXXX(SPINOR,EOT,VCT(1,2),GWF,AMP)
C--
C  That's it.
C--
      RETURN
      END
