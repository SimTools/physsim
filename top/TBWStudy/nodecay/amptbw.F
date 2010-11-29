C*(Update Record)
C*  2010/11/29  K.Fujii		original version.
C*

      SUBROUTINE AMPTBW(GAL,GAU,GAD,GZL,GZU,GZD,GWF,AMW,GMW,AMZ,GMZ,
     .                  AMT,GMT,AMB,GMB,EIN,EOT,BIN,TOT,WM,AMP)
 
      IMPLICIT   REAL*8  ( A-H, O-Z )
      REAL   *8  GAL(2), GAU(2), GAD(2), GZL(2), GZU(2), GZD(2), GWF(2),
     .           AMW, GMW, AMZ, GMZ, AMT, GMT, AMB, GMB
      COMPLEX*16 EIN(6), EOT(6), BIN(6), WM(6), TOT(6), AMP(0:3)
C--
      COMPLEX*16 TIN(6), BOT(6), WP(6), TMP(0:2,3)
C
C========< Entry Point >================================================
C
C--
C  Compute tt diagram.
C--
      CALL FVIXXX(BIN,WM,GWF,AMT,GMT,TIN)
      CALL AMPTT(GAL,GAU,GZL,GZU,AMZ,GMZ,EIN,EOT,TIN,TOT,TMP(0,1))
C--
C  Compute bb diagram.
C--
      CALL FVOXXX(TOT,WM,GWF,AMT,GMT,BOT)
      CALL AMPBB(GAL,GAD,GZL,GZD,AMZ,GMZ,EIN,EOT,BIN,BOT,TMP(0,2))
C--
C  Compute WW diagram.
C--
      CALL JIOXXX(BIN,TOT,GWF,AMW,GMW,WP)
      CALL AMPWW(GAL,GZL,GWF,GW,AMZ,GMZ,WM,WP,EIN,EOT,TMP(0,3))
C--
C  Sum two amplitudes.
C     AMP(0) = sum.
C        (1) = tt
C        (2) = bb
C        (3) = WW
C--
      AMP(1) = TMP(0,1)
      AMP(2) = TMP(0,2)
      AMP(3) = TMP(0,3)
      AMP(0) = AMP(1) + AMP(2) + AMP(3)
C--
C  That's it.
C--
      RETURN
      END
