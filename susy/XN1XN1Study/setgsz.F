      SUBROUTINE SETGSZ
 
      IMPLICIT   REAL*8  ( A-H, O-Z )
#include "SMPTAB.inc"
#include "SSPTAB.inc"
#include "SSCUPL.inc"
      COMPLEX*16 GSTX(2)
C
C========< Entry Point >================================================
C
#if 0
C--
C  This part is for usual MSUGRA.
C--
C  Select Ino: JX = 2 means chi^0_2.
C--
      JX = 2
C--
C  Calculate decay width for
C        SZJ -> SZ1 + Z
C--
      CALL F1F2VD(SZM(JX),SZM(1),AMZ,GNNZ(1,1,JX),GNNZ(2,1,JX),GMXSZ1)
C--
      PRINT *, '>>>> SETGSZ modified GMXTOT'
      PRINT *, ' GM(X->SZ1+Z)     = ', GMXSZ1, ' GeV'
C--
C  Calculate decay widths for
C        SZJ -> SeR+ + eR- & eR+ + SeR-
C            -> SmR+ + mR- & mR+ + SmR-
C--
      GSTX(1) = CONJG(GNESER(2,JX))
      GSTX(2) = CONJG(GNESER(1,JX))
C--
      CALL GMXFSF(GSTX,SZM(JX),AMF(1,2,1),SFM(3),GMXSER)
C--
      GMXSER = 2*GMXSER
      PRINT *, ' GM(X->seR+eR)    = ', GMXSER, ' GeV'
C--
      CALL GMXFSF(GSTX,SZM(JX),AMF(2,2,1),SFM(3),GMXSMR)
C--
      GMXSMR = 2*GMXSMR
      PRINT *, ' GM(X->smR+mR)    = ', GMXSMR, ' GeV'
#else
C--
C  This part is for GMSB study with stau LSP.
C--
C  Select Ino: JX = 1 means chi^0_1.
C--
      JX = 1
C--
C  No decay to itself.
C--
      GMXSZ1 = 0
C--
C  Override the 1st and 2nd generation slepton masses
C     SFM  (I): (SNL,SEL,SER,SUL,SUR,SDL,SDR)
C  while the 3rd generation masses 
C     SFM3M(I): (SNL,SE1,SE2,SU1,SU2,SD1,SD2)
C  are untouched.  Adjust AMSLEPL and AMSLEPR for your purpose.
C--
C      AMSLEPL = SFM(2)
      AMSLEPL = 400.
      AMSLEPR = 180.
      PRINT *, ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      PRINT *, ' CAUTION !!!! SETGSZ modified SFM(*) !!!!!'
      PRINT *, ' SFM  (2) = ', SFM(2), ' --> ', AMSLEPL
      PRINT *, ' SFM  (3) = ', SFM(3), ' --> ', AMSLEPR
      PRINT *, ' Stau1 mass stays the same:'
      PRINT *, ' SFM3M(2) = ', SFM3M(2)
      PRINT *, ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      SFM(2)  = AMSLEPL
      SFM(3)  = AMSLEPR
C--
C  Switch off the following #if block to make SZ1 -> tau stau 100%.
C--
#if 1
C--
C  Calculate decay width for
C        SZ1 -> SER+ + e- & e+ + SER-
C--
      GSTX(1) = CONJG(GNESER(2,JX))
      GSTX(2) = CONJG(GNESER(1,JX))
C--
      CALL GMXFSF(GSTX,SZM(JX),AMF(1,2,1),SFM(3),GMXSER)
      GMXSER = 2*GMXSER
      PRINT *, ' GM(X->seR+eR)    = ', GMXSER, ' GeV'
C--
C  Calculate decay width for
C        SZ1 -> SMR+ + mu- & mu+ + SMR-
C--
      CALL GMXFSF(GSTX,SZM(JX),AMF(2,2,1),SFM(3),GMXSMR)
      GMXSMR = 2*GMXSMR
      PRINT *, ' GM(X->smR+mR)    = ', GMXSMR, ' GeV'
C--
C  Calculate decay width for
C        SZ1 -> SEL+ + e- & e+ + SEL-
C--
      GSTX(1) = CONJG(GNESEL(2,JX))
      GSTX(2) = CONJG(GNESEL(1,JX))
C--
      CALL GMXFSF(GSTX,SZM(JX),AMF(1,2,1),SFM(2),GMXSEL)
      GMXSEL = 2*GMXSEL
      PRINT *, ' GM(X->seL+eL)    = ', GMXSEL, ' GeV'
C--
C  Calculate decay width for
C        SZ1 -> SML+ + mu- & mu+ + SML-
C--
      CALL GMXFSF(GSTX,SZM(JX),AMF(2,2,1),SFM(2),GMXSML)
      GMXSML = 2*GMXSML
      PRINT *, ' GM(X->smL+mL)    = ', GMXSML, ' GeV'
#else
C--
C  No cdecay to e + SER and mu + SMU in this "#else" block.
C--
      GMXSER = 0
      GMXSMR = 0
      GMXSEL = 0
      GMXSML = 0
#endif
#endif
C--
C  Calculate decay width for
C        SZ1 -> Stau1+ + tau- & tau+ + Stau1-
C--
      GSTX(1) = CONJG(GNESE1(2,JX))
      GSTX(2) = CONJG(GNESE1(1,JX))
C--
      CALL GMXFSF(GSTX,SZM(JX),AMF(3,2,1),SFM3M(2),GMXSTU)
C--
      GMXSTU = 2*GMXSTU
      PRINT *, ' GM(X->stau1+tau) = ', GMXSTU, ' GeV'
C-- 
C  Sum up partial widths.
C-- 
      GMSZ(JX) = GMXSZ1 + GMXSER + GMXSMR + GMXSTU
     .                  + GMXSEL + GMXSML
      PRINT *, ' ---'
      PRINT *, ' GMSZ(', JX, ') = ', GMSZ(JX), ' GeV'
C--
C  That's it.
C--
      RETURN
      END
