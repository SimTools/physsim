C*(Update Record)
C*  95/03/16  K.Fujii	New version with HELAS V204.
C*			New argument HSCLR(*) for Higgs wave function.
C*  95/03/28  K.Fujii	Modified the calling sequence of AMPZH 
C*			to be consistent with the one for ffh_gen.
C*
      SUBROUTINE AMPZH (GZL,GZH, AMZ,GMZ, ZVCT, HSCLR, EIN,EOT, AMP)
 
      IMPLICIT   REAL*4  ( A-H, O-Z )
      REAL   *4  GZL(2), GZH, AMZ, GMZ
      COMPLEX*8  EIN(6), EOT(6), ZVCT(6), HSCLR(3), AMP
C--
      COMPLEX*8  CURR(6)
C
C========< Entry Point >================================================
C
C--
C  Compute Z s-channel exchange diagram.
C--
      CALL JIOXXX(EIN,EOT,GZL,AMZ,GMZ,CURR)
      CALL VVSXXX(ZVCT,CURR,HSCLR,GZH,AMP)
C--
C  That's it.
C--
      RETURN
      END
