C* (Inputs)
C*    AMW    : (R*8) : W mass.
C*    AMU    : (R*8) : up type fermion mass.
C*    AMD    : (R*8) : down type fermion mass.
C*    VFF    : (R*8) : KM matrix element.
C*    C      : (R*8) : color factor.
C* (Output)
C*    GM     : (R*8) : partial width for W --> fu + fd.
C*
 
      SUBROUTINE GMW2FF(AMW,AMU,AMD,VFF,C,GM)
 
      IMPLICIT REAL*8 ( A-H, O-Z )
#include "include/SMCONS.inc"
      REAL*8 AMW, AMU, AMD, VFF, C, GM
C--
C  Statement function.
C--
      PSTR(A,B,C) = SQRT( (A-(B+C))*(A+(B+C))*(A-(B-C))*(A+(B-C)) )/2/A
C
C========< Entry Point >================================================
C
C--
C  Calculate width.
C--
      IF ( AMU+AMD.GE.AMW ) THEN
         GM = 0
         RETURN
      ELSE
         AMU2 = AMU*AMU
         AMD2 = AMD*AMD
         AMW2 = AMW*AMW
         P1P2 = ( AMW2 - AMU2 - AMD2 )/2
         P1   = PSTR(AMW,AMU,AMD)
      ENDIF
      TTA = ( P1P2 + 2*(AMU2+P1P2)*(AMD2+P1P2)/AMW2 )/3
      FAC = (xALF/xSIN2W)/2
      GM  = FAC*TTA*C*VFF*VFF*P1/AMW2
C--
C  That's it.
C--
      RETURN
      END
