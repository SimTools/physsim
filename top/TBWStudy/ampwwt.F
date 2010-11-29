CC**********************************************************************
C*
C*=====================================---===
C* Subroutine AMPWWT(GWF,WM,WP,EIN,EOT,AMP)
C*=====================================---===
C*
C* (Purpose)
C*    Calculate t-channel amplitude for e+ + e- ---> W+ + W-.
C*    Inputs are spinors and wave functions thus replaceable with
C*    spinors or wave function made from final stable particle
C*    4-momenta.
C* (Inputs)
C*       GWF(*)  : (R*4) : fermion W coupling.
C*       WM (*)  : (C*8) : W- wave function.
C*       WP (*)  : (C*8) : W+ wave function.
C*       EIN(*)  : (C*8) : incoming electron spinor.
C*       EOT(*)  : (C*8) : outgoing electron spinor.
C* (Output)
C*       AMP     : (C*8) : amplitude for t-channel diagram.
C* (Relation)
C*    Invokes subroutines in HELAS.LOAD.
C* (Update Record)
C*    90/08/08  K.Fujii         Original version.
C*    95/03/12  K.Fujii		Modiifed to HELAS V204.
C*  2010/11/29  K.Fujii		Changed arguments for W+ and W-.
C*
CC**********************************************************************
 
      SUBROUTINE AMPWWT(GWF,WM,WP,EIN,EOT,AMP)
 
      IMPLICIT    REAL*8  ( A-H, O-Z )
      REAL   *8   GWF(2)
      COMPLEX*16  EIN(6), EOT(6), WM(6), WP(6), AMP
C--
      COMPLEX*16  SPINOR(6)
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
      CALL FVIXXX(EIN,WM,GWF,AMN,GMN,SPINOR)
      CALL IOVXXX(SPINOR,EOT,WP,GWF,AMP)
C--
C  That's it.
C--
      RETURN
      END
