CC**********************************************************************
C*
C*=======================================================---===
C* Subroutine AMPWW(GAL,GZL,GWF,GW,AMZ,GMZ,WM,WP,EIN,EOT,AMP)
C*=======================================================---===
C*
C* (Purpose)
C*    Calculate amplitude for e+ + e- ---> W+ + W-.
C*    Inputs are spinors and wave functions thus replaceable with
C*    spinors or wave function made from final stable particle
C*    4-momenta.
C* (Inputs)
C*       GAL(*)  : (R*4) : lepton photon coupling.
C*       GZL(*)  : (R*4) : lepton Z coupling.
C*       GWF(*)  : (R*4) : fermion W coupling.
C*       GW      : (R*4) : weak coupling.
C*       AMZ     : (R*4) : Z mass.
C*       GMZ     : (R*4) : Z width.
C*       WM (*)  : (C*8) : W- wave function.
C*       WP (*)  : (C*8) : W+ wave function.
C*       EIN(*)  : (C*8) : incoming electron spinor.
C*       EOT(*)  : (C*8) : outgoing electron spinor.
C* (Output)
C*       AMP(0)  : (C*8) : amplitude sum.
C*          (1)  : (C*8) : sum of s-channel diagrams.
C*          (2)  : (C*8) : sum of t-channel diagrams.
C* (Relation)
C*    Invokes subroutines in HELAS.LOAD.
C* (Update Record)
C*    90/08/08  K.Fujii         Original version.
C*    95/03/12  K.Fujii		Modiifed to HELAS V204.
C*  2010/11/29  K.Fujii		Changed arguments for W+ and W-.
C*
CC**********************************************************************
 
      SUBROUTINE AMPWW(GAL,GZL,GWF,GW,AMZ,GMZ,WM,WP,EIN,EOT,AMP)
 
      IMPLICIT    REAL*8  ( A-H, O-Z )
      REAL   *8   GAL(2), GZL(2), GWF(2), GW, AMZ, GMZ
      COMPLEX*16  EIN(6), EOT(6), WM(6), WP(6), AMP(0:2)
C
C========< Entry Point >================================================
C
C--
C  Calculate the s-channel diagrams.
C--
      CALL AMPWWS(GAL,GZL,GW,AMZ,GMZ,WM,WP,EIN,EOT,AMP(1))
C--
C  Calculate the t-channel diagrams.
C--
      CALL AMPWWT(GWF,WM,WP,EIN,EOT,AMP(2))
C--
C  Sum s- and t-channel diagrams.
C--
      AMP(0) = AMP(1) + AMP(2)
C--
C  That's it.
C--
      RETURN
      END
