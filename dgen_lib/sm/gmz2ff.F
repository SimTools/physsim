C* (Inputs)
C*    AMZ    : (R*8) : Z mass.
C*    T3F    : (R*8) : T3 of fermion.
C*    QF     : (R*8) : fermion charge.
C*    C      : (R*8) : color factor.
C*    AM     : (R*8) : fermion mass.
C* (Output)
C*    GM     : (R*8) : partial width for Z --> f + f.
C*
 
      SUBROUTINE GMZ2FF(AMZ,T3F,QF,C,AM,GM)
 
      IMPLICIT REAL*8 ( A-H, O-Z )
#include "include/SMCONS.inc"
      REAL*8 AMZ, T3F, QF, C, AM, GM
C
C========< Entry Point >================================================
C
C--
C  Calculate width.
C--
      AMZ2 = AMZ*AMZ
      P1   = AMZ2/4 - AM*AM
      IF ( P1.LE.0.D0 ) THEN
         GM = 0
         RETURN
      ELSE
         P1   = SQRT(P1)
      ENDIF
      GV  =  T3F/2 - xSIN2W*QF
      GA  = -T3F/2
      TTA = 2*( (GV**2+GA**2)*( AMZ2 - 4*P1*P1/3 )
     .         - 4*GA*GA*AM*AM )
      FAC =  (xALF/xSIN2W/xCOS2W)/2
      GM  = FAC*TTA*C*P1/AMZ2
C--
C  That's it.
C--
      RETURN
      END
