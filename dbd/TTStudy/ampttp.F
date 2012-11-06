C*(Update Record)
C*  95/03/06  K.Fujii		New version with HELAS V204.
C*

      SUBROUTINE AMPTTP(GAL,GAU,GZL,GZU,AMZ,GMZ,
     .                  EIN,EOT,TIN,TOT,AMP)
 
      IMPLICIT   REAL*4  ( A-H, O-Z )
      REAL   *4  GAL(2), GAU(2), GZL(2), GZU(2), AMZ, GMZ
      COMPLEX*8  EIN(6), EOT(6), TIN(6), TOT(6), AMP(0:2)
C--
      REAL   *4  VZT(2), AZT(2)
      COMPLEX*8  CURR(6)
C
C========< Entry Point >================================================
C
C--
C  Compute photon s-channel exchange diagram.
C--
      CALL JIOXXX(EIN,EOT,GAL,0.,0., CURR)
      CALL IOVXXX(TIN,TOT,CURR,GAU,AMP(0))
C--
C  Compute Z s-channel exchange diagram.
C--
      VZT(1) = ( GZU(1) + GZU(2) )/2
      VZT(2) = VZT(1)
      CALL JIOXXX(EIN,EOT,GZL,AMZ,GMZ, CURR)
      CALL IOVXXX(TIN,TOT,CURR,VZT,AMP(1))
      AMP(1) = AMP(1) + AMP(0)
C--
      AZT(1) = ( GZU(1) - GZU(2) )/2
      AZT(2) = ( GZU(2) - GZU(1) )/2
      CALL JIOXXX(EIN,EOT,GZL,AMZ,GMZ, CURR)
      CALL IOVXXX(TIN,TOT,CURR,AZT,AMP(2))
C--
C  Sum two amplitudes.
C     AMP(0) = sum.
C        (1) = v_t part.
C        (2) = a_t part.
C--
      AMP(0) = AMP(1) + AMP(2)
C--
C  That's it.
C--
      RETURN
      END
