CC**********************************************************************
C*
C*======================================================---===
C* Subroutine AMPWWS(GAL,GZL,GWF,GW,AMZ,GMZ,VCT,EIN,EOT,AMP)
C*======================================================---===
C*
C* (Purpose)
C*    Calculate s-channel amplitude for e+ + e- ---> W+ + W-.
C*    Inputs are spinors and wave functions thus replaceable with
C*    spinors or wave function made from final stable particle
C*    4-momenta.
C* (Inputs)
C*       GAL(*)  : (R*4) : lepton photon coupling.
C*       GZL(*)  : (R*4) : lepton Z coupling.
C*       GW      : (R*4) : weak coupling.
C*       AMZ     : (R*4) : Z mass.
C*       GMZ     : (R*4) : Z width.
C*       WM (*)  : (C*8) : W- wave function.
C*       WP (*)  : (C*8) : W+ wave function.
C*       EIN(*)  : (C*8) : incoming electron spinor.
C*       EOT(*)  : (C*8) : outgoing electron spinor.
C* (Output)
C*       AMP     : (C*8) : amplitude for s-channel diagram.
C* (Relation)
C*    Invokes subroutines in HELAS.LOAD.
C* (Update Record)
C*    90/08/08  K.Fujii		Original version.
C*    95/03/12  K.Fujii		Modiifed to HELAS V204.
C*  2010/11/29  K.Fujii		Changed arguments for W+ and W-.
C*
CC**********************************************************************
 
      SUBROUTINE AMPWWS(GAL,GZL,GW,AMZ,GMZ,WM,WP,EIN,EOT,AMP)
 
      IMPLICIT    REAL*8  ( A-H, O-Z )
      REAL   *8   GAL(2), GZL(2), GW, AMZ, GMZ
      COMPLEX*16  EIN(6), EOT(6), WM(6), WP(6), AMP
C--
      COMPLEX*16  CURR(6)
C
C========< Entry Point >================================================
C
C--
C  Calculate the s-channel diagrams.
C--
      CALL J3XXXX(EIN,EOT,GAL,GZL,AMZ,GMZ,CURR)
      CALL VVVXXX(WM,WP,CURR,GW,AMP)
C--
C  That's it.
C--
      RETURN
      END
