CC**********************************************************************
C*
C*===========================================---===
C* Subroutine AMPZZ(GZL,AME,GME,ZVCT,EIN,EOT,AMP)
C*===========================================---===
C*
C* (Purpose)
C*    Calculate amplitudes for e+ + e- ---> Z0 + Z0.
C*    Inputs are spinors and wave functions thus replaceable with
C*    spinors or wave function made from final stable particle
C*    4-momenta.
C* (Inputs)
C*       GZL(*)   :(C*8): electron Z coupling.
C*       AME      :(R*4): electron mass.
C*       GME      :(R*4): electron width.
C*       ZVCT(*,i):(C*8): i = (1,2) = (Z0,Z0).
C*       EIN(*)   :(C*8): incoming electron spinor.
C*       EOT(*)   :(C*8): outgoing electron spinor.
C* (Output)
C*       AMP(0)   :(C*8): total amplitude.
C*          (1)   :(C*8): t-channel amplitude.
C*          (2)   :(C*8): u-channel amplitude.
C* (Relation)
C*    Invokes subroutines in HELAS.LOAD.
C* (Update Record)
C*    95/03/12  K.Fujii		Modiifed to HELAS V204.
C*
CC**********************************************************************

      SUBROUTINE AMPZZ (GZL, AME,GME, ZVCT, EIN,EOT, AMP)
 
      IMPLICIT   REAL*4  ( A-H, O-Z )
      REAL   *4  GZL(2), AME, GME
      COMPLEX*8  EIN(6), EOT(6), ZVCT(6,2), AMP(0:2)
C--
      COMPLEX*8  SPWRK(6,2)
C
C========< Entry Point >================================================
C
C--
C  Compute t-channel exchange diagram.
C--
      CALL FVIXXX(EIN,ZVCT(1,1),GZL,AME,GME,SPWRK(1,1))
      CALL IOVXXX(SPWRK(1,1),EOT,ZVCT(1,2),GZL,AMP(1))
C--
C  Compute u-channel exchange diagram.
C--
      CALL FVIXXX(EIN,ZVCT(1,2),GZL,AME,GME,SPWRK(1,1))
      CALL IOVXXX(SPWRK(1,1),EOT,ZVCT(1,1),GZL,AMP(2))
C--
C  Sum up t- and u-channel exchange diagrams.
C--
      AMP(0) = AMP(1) + AMP(2)
C--
C  That's it.
C--
      RETURN
      END
