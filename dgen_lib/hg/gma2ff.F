C* (Inputs)
C*    AMH    : (R*4) : pseudo-scalar Higss mass.
C*    AMF    : (R*4) : fermion mass.
C*    AMZ    : (R*4) : Z mass.
C* (Output)
C*    GM     : (R*4) : (1/1)*( Higss partical width ) for F = L.
C*                   : (1/3)*( Higss partical width ) for F = Q.
C*                   : In both cases, one should multiply GM by
C*                   : tan(beta)**2 for f=D, or cot(beta)**2
C*                   : for f=U.
C*
 
      SUBROUTINE GMA2FF(AMH,AMF,AMZ,GM)
 
      IMPLICIT REAL*8 ( A-H, O-Z    )
#include "../sm/include/SMCONS.inc"
      REAL*8 AMH, AMF, AMZ, GM
C
C========< Entry Point >================================================
C
C--
C  Initialize constants.
C--
      FACT   = xALF/xSIN2W/xCOS2W/AMZ**2/8
C--
C  Calculate width.
C--
      EF = AMH/2
      BT = (EF-AMF)*(EF+AMF)
      IF ( BT.LE.0.D0 ) THEN
         GM = 0
         RETURN
      ELSE
         BT   = SQRT(BT)/EF
      ENDIF
      GM = FACT*AMF*AMF*AMH*BT
C--
C  That's it.
C--
      RETURN
      END
