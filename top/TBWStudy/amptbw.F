C*(Update Record)
C*  2010/11/29  K.Fujii		original version.
C*

      SUBROUTINE AMPTBW(GAL,GAU,GAD,GZL,GZU,GZD,GWF,GW,AMW,GMW,AMZ,GMZ,
     .                  AMT,GMT,AMB,GMB,
     .                  EIN,EOT,BIN,TOT,WM,BOT,TIN,WP, AMP)
 
      IMPLICIT   REAL*8  ( A-H, O-Z )
      REAL   *8  GAL(2), GAU(2), GAD(2), GZL(2), GZU(2), GZD(2), GWF(2),
     .           GW, AMW, GMW, AMZ, GMZ, AMT, GMT, AMB, GMB
      COMPLEX*16 EIN(6), EOT(6), BIN(6), WM(6), TOT(6), 
     .           BOT(6), TIN(6), WP(6), AMP(0:3)
C--
      COMPLEX*16 TINP(6), BOTP(6), TOTP(6), BINP(6), 
     .           WMP(6), WPP(6), TMP(0:2,5)
C
C========< Entry Point >================================================
C
C--
C  Compute tt diagram.
C--
      CALL AMPTT(GAL,GAU,GZL,GZU,AMZ,GMZ,EIN,EOT,TIN,TOT,TMP(0,1))
C--
C  Compute bb diagram.
C--
      CALL FVOXXX(TOT,WM,GWF,AMB,GMB,BOTP)
      CALL AMPBB(GAL,GAD,GZL,GZD,AMZ,GMZ,EIN,EOT,BIN,BOTP,TMP(0,2))
      CALL FVIXXX(TIN,WP,GWF,AMB,GMB,BINP)
      CALL AMPBB(GAL,GAD,GZL,GZD,AMZ,GMZ,EIN,EOT,BINP,BOT,TMP(0,3))
C--
C  Compute WW diagram.
C--
      CALL JIOXXX(BIN,TOT,GWF,AMW,GMW,WPP)
      CALL AMPWW(GAL,GZL,GWF,GW,AMZ,GMZ,WM,WPP,EIN,EOT,TMP(0,4))
      CALL JIOXXX(TIN,BOT,GWF,AMW,GMW,WMP)
      CALL AMPWW(GAL,GZL,GWF,GW,AMZ,GMZ,WMP,WP,EIN,EOT,TMP(0,5))
C--
C  Sum two amplitudes.
C     AMP(0) = sum.
C        (1) = tt
C        (2) = bb
C        (3) = WW
C--
      AMP(1) = TMP(0,1)
      AMP(2) = TMP(0,2) + TMP(0,3)
      AMP(3) = TMP(0,4) + TMP(0,5)
      AMP(0) = AMP(1) + AMP(2) + AMP(3)
C--
C  That's it.
C--
      RETURN
      END
