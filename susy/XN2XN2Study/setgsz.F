      SUBROUTINE SETGSZ
 
      IMPLICIT   REAL*8  ( A-H, O-Z )
#include "SMPTAB.inc"
#include "SSPTAB.inc"
#include "SSCUPL.inc"
      COMPLEX*16 GSTX(2)
C
C========< Entry Point >================================================
C
C--
C  Select Ino.
C--
      JX = 2
C--
C  Prepare coupling constants.
C        SZ2 -> SZ1 + Z
C--
      CALL F1F2VD(SZM(JX),SZM(1),AMZ,GNNZ(1,1,JX),GNNZ(2,1,JX),GMXSZ1)
C--
      PRINT *, '>>>> SETGSZ modified GMXTOT'
      PRINT *, ' GM(X->SZ1+Z)     = ', GMXSZ1, ' GeV'
C--
C  Prepare coupling constants.
C        SZ2 -> SeR+ + eR- & eR+ + SeR-
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
C--
C  Prepare coupling constants.
C        SZ2 -> Stau1+ + tau- & tau+ + Stau1-
C--
      GSTX(1) = CONJG(GNESE1(2,JX))
      GSTX(2) = CONJG(GNESE1(1,JX))
C--
      CALL GMXFSF(GSTX,SZM(JX),AMF(3,2,1),SFM3M(2),GMXSTU)
C--
      GMXSTU = 2*GMXSTU
      PRINT *, ' GM(X->stau1+tau) = ', GMXSTU, ' GeV'
C-- 
      GMSZ(JX) = GMXSZ1 + GMXSER + GMXSMR + GMXSTU
      PRINT *, ' GMSZ(', JX, ') = ', GMSZ(JX), ' GeV'
C--
C  That's it.
C--
      RETURN
      END
