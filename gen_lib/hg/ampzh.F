      SUBROUTINE AMPZH (GZL,GVH,AMZ,GMZ,ZVCT,HOT,EIN,EOT,AMP)
 
      IMPLICIT   REAL*4  ( A-H, O-Z )
      REAL   *4  GZL(2), GVH, AMZ, GMZ
      COMPLEX*8  EIN(6), EOT(6), ZVCT(6), HOT(3), AMP
C--
      COMPLEX*8  CURR(6)
C
C========< Entry Point >================================================
C
C--
C  Compute photon s-channel exchange diagram.
C--
      CALL JIOXXX(EIN,EOT,GZL,AMZ,GMZ,CURR)
      CALL VVSXXX(CURR,ZVCT,HOT,GVH,AMP)
C--
C  That's it.
C--
      RETURN
      END
