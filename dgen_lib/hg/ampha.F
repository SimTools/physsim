C*(Notice)
C*  GCHA are defined by out-going particle charges.
C*(Update Record)
C*  95/03/30  K.Fujii	GHA was made C*8 and renamed GCHA.
C*			VSSXXX was replaced by VSSXXC, accordingly.
C*
 
      SUBROUTINE AMPHA(GZL,GCHA,AMZ,GMZ,HIN,AOT,EIN,EOT,AMP)
 
      IMPLICIT   REAL*8  ( A-H, O-Z )
      REAL   *8  GZL(2), AMZ, GMZ
      COMPLEX*16 GCHA, EIN(6), EOT(6), HIN(3), AOT(3), AMP
C--
      COMPLEX*16 CURR(6)
C
C========< Entry Point >================================================
C
C--
C  Compute Z s-channel exchange diagram.
C--
      CALL JIOXXX(EIN,EOT,GZL,AMZ,GMZ, CURR)
      CALL VSSXXC(CURR,AOT,HIN,GCHA,AMP)
C--
C  That's it.
C--
      RETURN
      END
