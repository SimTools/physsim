      SUBROUTINE USROUT

#include "USRPRM.inc"
#include "USMPRM.inc"
#include "BS51PR.inc"
      COMMON /PLOTLU/ LU
      INTEGER*4       LU
C--
C  Set logical units and open the histogram file.
C--
      LOU = LUHIST
      LU  = LUINFO
      OPEN(LUHIST, FILE=BSHSTF)
C--
C  Output integration results.
C-- 
      LURSLT = 30
      WRITE(LURSLT,'(5E12.5)') SQRTS ,ESTIM, SIGMA, AMST, AMSH
C--
C  Output histograms.
C-- 
      DO 100 IH = 1, 19
         CALL XHSAVE( LOU,  IH )
100   CONTINUE       
C--
C  That's it.
C-- 
      RETURN
      END
