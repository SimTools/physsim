C*
C* (Inputs)
C*    AMT    : (R*4) : top mass.
C*    AMD    : (R*4) : down type quark mass.
C*    AMW    : (R*4) : W mass.
C*    VFF    : (R*4) : KM matrix element.
C* (Output)
C*    GM     : (R*4) : partial width for t --> D + W.
C*
 
      SUBROUTINE GMT2QW(AMT,AMD,AMW,VFF,GM)
 
      IMPLICIT REAL*4 ( A-H, O-Z )
#include "../sm/include/SMCONS.inc"
      REAL*4 AMT, AMD, AMW, VFF,  GM
C--
C  Statement function.
C--
      BETA(X1,X2) = SQRT( 1 - 2*(X1+X2) + (X1-X2)**2 )
C
C========< Entry Point >================================================
C
C--
C  Calculate width.
C--
      IF ( AMD+AMW.GE.AMT ) THEN
         GM = 0
         RETURN
      ENDIF
      X1  = (AMD/AMT)**2
      X2  = (AMW/AMT)**2
      FAC = (xALF/xSIN2W)/16
      GM  = FAC*(AMT**3/AMW**2)*(VFF**2)
     .         *BETA(X1,X2)*( (1-X1)**2 + X2*(1+X1) - 2*X2*X2 )
C--
C  That's it.
C--
      RETURN
      END
