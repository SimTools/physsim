C*
C* (Inputs)
C*    AMT    : (R*4) : top mass.
C*    AMD    : (R*4) : down type quark mass.
C*    AMHM   : (R*4) : H- mass.
C*    AMW    : (R*4) : W mass.
C* (Output)
C*    GM     : (R*4) : partial width for t --> D + W.
C*
 
      SUBROUTINE GMT2HB(AMT,AMD,AMHM,AMW,TNB,GM)
 
      IMPLICIT REAL*8 ( A-H, O-Z    )
#include "../sm/include/SMCONS.inc"
      REAL*8 AMT, AMD, AMHM, AMW, TNB, GM
C--
      BETA(X1,X2) = SQRT( 1 - 2*(X1+X2) + (X1-X2)**2 )
C
C========< Entry Point >================================================
C
C--
C  Initialize constants.
C--
      FACT   = -xALF/xSIN2W/AMW**2/32
C--
C  Calculate width.
C--
      IF ( AMT.LE.AMHM+AMD) THEN
         GM = 0
         RETURN
      ELSE
         BT = BETA((AMD/AMT)**2,(AMHM/AMT)**2)
      ENDIF
      A2 = ( AMT/TNB + AMD*TNB )**2
      B2 = ( AMT/TNB - AMD*TNB )**2
      GM = FACT*( A2*(AMHM-AMT-AMD)*(AMHM+AMT+AMD)
     .          + B2*(AMHM-AMT+AMD)*(AMHM+AMT-AMD) )*BT/AMT
C--
C  That's it.
C--
      RETURN
      END
