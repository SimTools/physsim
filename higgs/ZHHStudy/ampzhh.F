C*(Update Record)
C*  08/11/30  K.Fujii		Original version.
C*
      SUBROUTINE AMPZHH (GZZH,GHHH,GZZHH, GZL, 
     .                   AMZ,GMZ,AMH,GMH, ZVCT, H1, H2, EIN, EOT, AMP)
 
      IMPLICIT   REAL*4  ( A-H, O-Z )
      REAL   *4  GZZH, GHHH, GZZHH, GZL(2), AMZ, GMZ, AMH, GMH
      COMPLEX*8  EIN(6), EOT(6), ZVCT(6), H1(3), H2(3), AMP(0:3)
C--
      COMPLEX*8  CURR(6), SCURR(6), VCURR(6), TMP
C
C========< Entry Point >================================================
C
C--
C  Compute Z s-channel exchange diagram.
C--
      CALL JIOXXX(EIN,EOT,GZL,AMZ,GMZ,CURR)
C
C (1) self-coupling
C
      CALL HSSXXX(H1,H2,GHHH,AMH,GMH,SCURR)
      CALL VVSXXX(ZVCT,CURR,SCURR,GZZH,AMP(1))
C
C (2) 2 H's from Z
C
      CALL JVSXXX(ZVCT,H1,GZZH,AMZ,GMZ,VCURR)
      CALL VVSXXX(VCURR,CURR,H2,GZZH,TMP)

      CALL JVSXXX(ZVCT,H2,GZZH,AMZ,GMZ,VCURR)
      CALL VVSXXX(VCURR,CURR,H1,GZZH,AMP(2))

      AMP(2) = AMP(2) + TMP
C
C (3) 4-point HHZZ
C
      CALL VVSSXX(ZVCT,CURR,H1,H2,GZZHH,AMP(3))
C
      AMP(0) = AMP(1) + AMP(2) + AMP(3)
C--
C  That's it.
C--
      RETURN
      END
