      SUBROUTINE USROUT

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
C  Output histograms.
C-- 
       CALL XHSAVE(LOU, 1 )
       CALL XHSAVE(LOU, 2 )
       CALL XHSAVE(LOU, 3 )
       CALL XHSAVE(LOU, 4 )
       CALL XHSAVE(LOU, 5 )
       CALL XHSAVE(LOU, 6 )
       CALL XHSAVE(LOU, 7 )
       CALL XHSAVE(LOU, 8 )
       CALL XHSAVE(LOU, 9 )
C--
C  That's it.
C-- 
      RETURN
      END
