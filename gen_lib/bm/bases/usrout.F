      SUBROUTINE USROUT
 
      INCLUDE 'USRPRM.inc'
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
      CALL XHSAVE( 24, 1 )
      CALL XHSAVE( 24, 2 )
      CALL XHSAVE( 24, 3 )
      CALL XHSAVE( 24, 4 )
      CALL XHSAVE( 24, 5 )
C--
C  That's it.
C-- 
      RETURN
      END
